//! Functionality directly related to demultiplexing a [`RecordSet`] of FASTQ records.

use std::{borrow::Cow, iter::IntoIterator, vec::Vec};

use anyhow::{bail, ensure, Context, Result};
use bstr::ByteSlice;
use itertools::Itertools;
use log::debug;
use read_structure::{ReadSegment, ReadStructure, SegmentType};
use seq_io::{
    fastq::{OwnedRecord, Record, RecordSet, RefRecord},
    BaseRecord,
};

use crate::{
    fastq_header::FastqHeader,
    matcher::{MatchResult, Matcher},
    metrics::{BaseQualCounter, DemuxedGroupMetrics, SampleBarcodeHopChecker},
    sample_metadata::SampleMetadata,
    utils::MultiZip,
};

#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
struct Barcode<'a>(Cow<'a, [u8]>);

/// A helper struct to build up the extracted barcode sequences in both concatenated and delimited form for sample and UMI barcodes.
#[derive(Debug, Default)]
struct ExtractedBarcode {
    /// The extracted barcode in a concatenated format (no delimiters).
    concatenated: Vec<u8>,
    /// The extracted barcode in a delimited format, with `+` as the delimiter.
    delimited: Vec<u8>,
}

impl ExtractedBarcode {
    /// Create a new empty [`ExtractedBarcode`].
    fn new() -> Self {
        Self::default()
    }

    /// Create a new [`ExtractedBarcode`] from the delimited barcode, with `+` as the delimiter.
    fn from_delimited(delimited: &[u8]) -> Self {
        let concatenated = delimited.splitn(2, |c| *c == b'+').collect_vec().concat();
        ExtractedBarcode { concatenated, delimited: delimited.to_vec() }
    }

    /// Add bases to the [`ExtractedBarcode`], this will add to both the concatenated and delimited forms.
    fn add_bases(&mut self, bases: &[u8]) {
        if self.non_empty() {
            self.delimited.push(b'+');
        }
        self.delimited.extend(bases);
        self.concatenated.extend(bases);
    }

    /// Check if the [`ExtractedBarcode`] is empty or not.
    #[inline]
    fn is_empty(&self) -> bool {
        self.concatenated.is_empty() && self.delimited.is_empty()
    }

    /// Check if the [`ExtractedBarcode`] is **not** empty.
    #[inline]
    fn non_empty(&self) -> bool {
        !self.is_empty()
    }

    /// Count the number of `N`s in the barcode
    // TODO: verify this is faster than `naive_count_32` or just filter -> count
    fn num_no_calls(&self) -> usize {
        bytecount::count(&self.concatenated, b'N')
    }

    /// Take ownership of the two inner vecs: `(concatenated, delimited)`
    fn into_inner(self) -> (Vec<u8>, Vec<u8>) {
        (self.concatenated, self.delimited)
    }
}

/// A [`PerFastqRecordSet`] is a a collection of [`RecordSet`], one for each input FASTQ be demultiplexed.
#[derive(Debug)]
pub struct PerFastqRecordSet {
    /// The record sets, one per FASTQ
    per_raw_fastq_reads: Vec<RecordSet>,
}

impl PerFastqRecordSet {
    /// Create a new [`PerFastqRecordSet`] from a [`Vec`] of [`RecordSet`].
    ///
    /// Each record set must contain the same number of reads.
    pub fn new(per_raw_fastq_reads: Vec<RecordSet>) -> Result<Self> {
        ensure!(
            per_raw_fastq_reads.iter().all(|s| s.len() == per_raw_fastq_reads[0].len()),
            "Unequal number of reads in each record set (likely uneven input FASTQs)"
        );
        Ok(Self { per_raw_fastq_reads })
    }
}

/// Helper struct to hold onto and reuse the filter parameters that will used bye the [`Demultiplexer`].
#[derive(Debug)]
pub struct DemuxReadFilterConfig {
    /// If true, filter out reads that have the control field set in the FASTQ header
    pub filter_control_reads: bool,
    /// If true, filter out reads that have the failing quality filter field set in the FASTQ heder
    pub filter_failing_quality: bool,
    /// The quality thresholds at which to mask template bases to N. Sample barcode/index
    /// and UMI bases are never masked. The vec must contain one entry per input FASTQ file being
    /// masked.  If a base quality in a read is less than the corresponding threshold value for the
    /// FASTQ then the base will be masked. A `quality_mask_threshold` of 0 indicates that no
    /// masking checks shall occur.
    pub quality_mask_thresholds: Vec<u8>,
    /// Max no-calls (N's) in a barcode before it is considered unmatchable.
    ///
    /// A barcode with total N's greater than `max_no_call` will be considered unmatchable.
    pub max_no_calls: Option<usize>,
}

impl DemuxReadFilterConfig {
    /// Create a new [`DemuxReadFilterConfig`].
    pub fn new(
        filter_control_reads: bool,
        filter_failing_quality: bool,
        quality_mask_thresholds: Vec<u8>,
        max_no_calls: Option<usize>,
    ) -> Self {
        Self { filter_control_reads, filter_failing_quality, quality_mask_thresholds, max_no_calls }
    }
}

impl Default for DemuxReadFilterConfig {
    /// The defaults for [`DemuxReadFilterConfig`] are all set to perform no filtering or masking.
    fn default() -> Self {
        Self {
            filter_control_reads: false,
            filter_failing_quality: false,
            quality_mask_thresholds: vec![],
            max_no_calls: None,
        }
    }
}

/// [`ReadSegment`]s grouped together by [`SegmentType`].
struct GroupedReadSegments {
    segments_by_kind: Vec<(SegmentType, Vec<ReadSegment>)>,
}

/// The reason that a read has failed the filter step
#[non_exhaustive]
enum ReadFailedFilterResult {
    QualityFilter,
    ControlFilter,
}

/// A trait that defines what is necessary to demultiplex a [`PerFastqRecordSet`].
pub trait Demultiplex: Send + Sync {
    /// Demultiplex a [`PerFastqRecordSet`].
    fn demultiplex(&self, record_set: &PerFastqRecordSet) -> Result<DemuxedGroup>;
    /// A function that should be available to the [`Demultiplex::demultiplex`] function to check unmatched reads for index hopping.
    fn demux_hop_checker(&self) -> Option<&SampleBarcodeHopChecker>;
    /// Get the index of the undetermined sample.
    fn get_undetermined_index(&self) -> usize;
    /// Get the list of samples being demultiplexed into, including the undetermined sample.
    fn samples(&self) -> &[SampleMetadata];
}

/// An implementation of [`Demultiplex`]  that supports checking for index hopping.
///
/// This implementation is generic over the [`Matcher`] that is used.
pub struct Demultiplexer<'a, M: Matcher> {
    /// The ordered list of samples.
    ///
    /// **Note**: the last sample in this list is expected to be the undetermined sample.
    samples: &'a [SampleMetadata],
    /// The index of the undetermined sample in `samples`.
    undetermined_sample_index: usize,
    /// The [`ReadSegment`]s grouped by [`SegmentType`].
    grouped_read_segments: Vec<GroupedReadSegments>,
    /// The [`SegmentType`]s to write to output FASTQs.
    output_segment_types: &'a [SegmentType],
    /// The number of separate output files, derived from the `output_segment_type` and `read_structure`.
    num_outputs: usize,
    /// The criteria for filtering out reads or masking bases.
    read_filter_config: &'a DemuxReadFilterConfig,
    /// The matcher to use for determining which sample the found barcode best matches to.
    matcher: M,
    /// An optional index hop checker that will check unmatched barcodes for index hopping.
    ///
    /// If this is present then an [`crate::metrics::SampleBarcodeHopTracker`] will be created on each [`DemuxedGroupMetrics`].
    sample_barcode_hop_checker: Option<SampleBarcodeHopChecker>,
    /// If this is true, the unmatched barcodes will be collected on a vec in [`DemuxedGroup`], and then sent
    /// to the unmatched barcode counter thread later in demultiplexing.
    ///
    /// If false, unmatched barcodes will not be collected, since there will be no thread to send to if unmatched are not
    /// being collected.
    collect_unmatched: bool,
    /// If this is true, then the read names across FASTQs will not be enforced to be the same.  This may be useful when
    /// the read names are known to be the same and performance matters.
    skip_read_name_check: bool,
    /// If this is true, then the sample barcode is expected to be in the FASTQ read header.  For
    /// dual indexed data, the barcodes must be `+` (plus) delimited.  Additionally, if true, then
    /// neither index FASTQ files nor sample barcode segments in the read structure may be
    /// specified.
    sample_barcode_in_fastq_header: bool,
}

impl<'a, M> Demultiplexer<'a, M>
where
    M: Matcher,
{
    /// Create a new [`Demultiplexer`].
    ///
    /// **Note** the last sample in `samples` _must_ be the undetermined sample.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        samples: &'a [SampleMetadata],
        read_structures: &'a [ReadStructure],
        output_segment_types: &'a [SegmentType],
        read_filter_config: &'a DemuxReadFilterConfig,
        matcher: M,
        collect_unmatched: bool,
        skip_read_name_check: bool,
        sample_barcode_in_fastq_header: bool,
    ) -> Result<Self> {
        // TODO: use display instead of debug formatting
        ensure!(
            read_structures.iter().any(|structure| structure.templates().count() > 0),
            "No template sequence found in read structures: {:?}",
            read_structures
        );

        let num_outputs = Self::determine_num_outputs(read_structures, output_segment_types);
        let sample_barcode_hop_checker = SampleBarcodeHopChecker::try_new(samples, read_structures)
            .context("Failed to create Index Hop Checker")?;
        let read_structures: Vec<GroupedReadSegments> = read_structures
            .iter()
            .cloned()
            .map(|s| GroupedReadSegments {
                segments_by_kind: s
                    .iter()
                    // Sort by SegmentType for group_by
                    .sorted_by_key(|segment| &segment.kind)
                    .group_by(|segment| segment.kind)
                    .into_iter()
                    .map(|(kind, group)| (kind, group.copied().collect()))
                    .collect(),
            })
            .collect();
        Ok(Self {
            samples,
            undetermined_sample_index: samples.len() - 1,
            grouped_read_segments: read_structures,
            output_segment_types,
            num_outputs,
            read_filter_config,
            matcher,
            sample_barcode_hop_checker,
            collect_unmatched,
            skip_read_name_check,
            sample_barcode_in_fastq_header,
        })
    }

    /// Determine how many outputs there will be given the read structures and output types.
    ///
    /// There will be one output per [`ReadStructure`] per [`SegmentType`] that is present in `output_types`.
    ///
    /// # Example
    ///
    /// A `read_structure` of [10B100T20S100T 150T] with `output_types` of [`SegmentType::Template`] and
    /// [`SegmentType::SampleBarcode`] will yield 3 outputs:
    ///
    /// - 1 for the Barcode sequence in structure 1
    /// - 1 for the two Template sequences in structure 1
    /// - 1 for the Template sequence in structure 2
    ///
    /// A `read_structure` of [10B100T20S100T 150T 10B] with `output_types` of [`SegmentType::Template`] and
    /// [`SegmentType::SampleBarcode`] will yield 4 outputs:
    ///
    /// - 1 for the Barcode sequence in structure 1
    /// - 1 for the two Template sequences in structure 1
    /// - 1 for the Template sequence in structure 2
    /// - 1 for the Barcode sequence in structure 3
    ///
    /// A `read_structure` of [10M5S+T, 8B, 8B, 10M5S+T] with `output_types` of [`SegementType::Template`] will
    /// yield 2 outputs:
    ///
    /// - 1 for the Template sequence in structure 1
    /// - 1 for the Template sequence in structure 4
    fn determine_num_outputs(
        read_structures: &[ReadStructure],
        output_types: &[SegmentType],
    ) -> usize {
        read_structures
            .iter()
            .map(|r| {
                r.iter()
                    .map(|s| s.kind)
                    .sorted()
                    .dedup()
                    .filter(|k| output_types.contains(k))
                    .count()
            })
            .sum()
    }

    /// Apply any filters that have been configured.
    ///
    /// If `Some(ReadFailedFilterResult)` is returned, the grouped reads should be skipped.
    ///
    /// This will increment the filter counts.
    fn apply_filters(&self, zipped_reads: &[RefRecord]) -> Result<Option<ReadFailedFilterResult>> {
        if self.read_filter_config.filter_control_reads
            || self.read_filter_config.filter_failing_quality
        {
            for rec in zipped_reads.iter() {
                let header = FastqHeader::try_from(rec.head())?;
                if self.read_filter_config.filter_failing_quality
                    && header.is_failed_quality_filter()
                {
                    return Ok(Some(ReadFailedFilterResult::QualityFilter));
                }
                if self.read_filter_config.filter_control_reads && header.is_control() {
                    return Ok(Some(ReadFailedFilterResult::ControlFilter));
                }
            }
        }
        Ok(None)
    }

    /// Process all bases and qualities for segments of the given `kind` from a `read`.
    ///
    /// Processing here entails masking low quality bases, extracting the segments that match the given `kind`
    /// and generating an [`OwnedRecord`] with the processed and extracted bases and qualities.
    ///
    /// # Arguments
    /// - `kind` - the [`SegmentType`] being extracted.
    /// - `read` - the [`RefRecord`] the segments are being extracted from.
    /// - `segments` - all segments of the same [`SegmentType`] based on the [`ReadStructure`] for the given read.
    /// - `sample_barcode` - the [`ExtractedBarcode`] that being built up as the sample barcode across reads.
    /// - `umi_barcode` - the [`ExtractedBarcode`] that being built up as the umi barcode across reads.
    ///
    /// # Returns
    ///
    /// Assuming no errors were encountered, a tuple is returned with the first value being an optional [`OwnedRecord`].
    /// If the `kind` of segments being processed is one of the [`Self::output_segment_types`], then an [`OwnedRecord`]
    /// will be created. Otherwise it will be none.
    ///
    /// The second value in the tuple is the [`BaseQualCounter`] which aggregates metrics about the base qualities
    /// in the segments being processed.
    #[inline]
    fn process_segments(
        &self,
        kind: SegmentType,
        read: &RefRecord,
        segments: &[ReadSegment],
        quality_mask_threshold: u8,
        sample_barcode: &mut ExtractedBarcode,
        umi_barcode: &mut ExtractedBarcode,
    ) -> Result<(Option<OwnedRecord>, BaseQualCounter)> {
        let mut qual_counter = BaseQualCounter::new();
        let mut builder_read = if self.output_segment_types.contains(&kind) {
            let mut read = read.to_owned_record();
            read.seq.clear();
            read.qual.clear();
            Some(read)
        } else {
            None
        };

        for segment in segments {
            let (bases, quals) =
                segment.extract_bases_and_quals(read.seq(), read.qual()).with_context(|| {
                    format!(
                        "Failed to extract bases and quality scores for {:?} from {:?}",
                        segment,
                        String::from_utf8_lossy(read.head())
                    )
                })?;

            // Mask low quality bases
            let do_masking = quality_mask_threshold > 0 && kind == SegmentType::Template;
            let sub_read_bases = if do_masking {
                let mut new_bases = Vec::with_capacity(bases.len());
                for (base, qual) in
                    bases.iter().zip(qual_counter.update_with_iter(quals.iter(), &kind))
                {
                    new_bases.push(if qual - 33 < quality_mask_threshold { b'N' } else { *base });
                }
                Cow::from(new_bases)
            } else {
                qual_counter.update(quals, &kind);
                Cow::from(bases)
            };

            if kind == SegmentType::SampleBarcode {
                sample_barcode.add_bases(bases);
            } else if kind == SegmentType::MolecularBarcode {
                umi_barcode.add_bases(bases);
            }

            if let Some(builder) = builder_read.as_mut() {
                builder.seq.extend_from_slice(&sub_read_bases);
                builder.qual.extend_from_slice(quals);
            }
        }

        Ok((builder_read, qual_counter))
    }
}

impl<'a, M> Demultiplex for Demultiplexer<'a, M>
where
    M: Matcher + Send + Sync,
{
    /// Demultiplex the [`PerFastqRecordSet`].
    #[allow(clippy::too_many_lines)]
    fn demultiplex(&self, record_set: &PerFastqRecordSet) -> Result<DemuxedGroup> {
        let mut demuxed_group = DemuxedGroup::new(
            self.samples.len(),
            self.num_outputs,
            self.demux_hop_checker(),
            self.collect_unmatched,
        );

        ensure!(
            record_set.per_raw_fastq_reads.len() == self.grouped_read_segments.len(),
            "Number of input FASTQs ({}) does not match number expected from read structures ({})",
            record_set.per_raw_fastq_reads.len(),
            self.grouped_read_segments.len()
        );

        let iterators =
            record_set.per_raw_fastq_reads.iter().map(IntoIterator::into_iter).collect();

        // Iterate over the reads within each chunk, zipped together
        for zipped_reads in MultiZip::new(iterators) {
            // Iterate over each read + read structure combo, collecting barcodes
            let mut sample_barcode = ExtractedBarcode::new();
            let mut umi_barcode = ExtractedBarcode::new();
            let mut output_reads = vec![];
            let mut qual_counter = BaseQualCounter::new();

            if !self.skip_read_name_check {
                let first_head = zipped_reads[0].head();
                let end_index = zipped_reads[0].head().find_byte(b' ').unwrap_or(first_head.len());
                for read in zipped_reads.iter().dropping(1) {
                    let cur_head = read.head();
                    let ok = cur_head.len() == end_index
                        || (cur_head.len() > end_index && cur_head[end_index] == b' ');
                    let ok = ok && first_head[0..end_index] == cur_head[0..end_index];
                    ensure!(
                        ok,
                        "Read names did not match: {:?} != {:?}",
                        String::from_utf8_lossy(first_head),
                        String::from_utf8_lossy(cur_head)
                    );
                }
            }

            // If any filtering was specified, check it here by peeking at the header of the first read in the set of reads.
            match self.apply_filters(&zipped_reads)? {
                Some(ReadFailedFilterResult::ControlFilter) => {
                    demuxed_group.group_metrics.num_reads_filtered_as_control += 1;
                    continue;
                }
                Some(ReadFailedFilterResult::QualityFilter) => {
                    demuxed_group.group_metrics.num_reads_failing_filters += 1;
                    continue;
                }
                None => (),
            }

            // If desired, get the sample barcode bases from the FASTQ header (of the first read)
            if self.sample_barcode_in_fastq_header {
                // get the FASTQ header from the first read
                let read = &zipped_reads[0];
                let header =
                    FastqHeader::try_from(read.head()).with_context(|| {
                        format!(
                            "Unable to parse read header: {}",
                            String::from_utf8_lossy(read.head()),
                        )
                    })?;
                // get the sample barcode from the FASTQ header
                match header.sample_barcode() {
                    None => bail!(
                        "Could not read the sample barcode from the read header: {}",
                        String::from_utf8_lossy(read.head())
                    ),
                    Some(barcode) => sample_barcode = ExtractedBarcode::from_delimited(barcode),
                };
            }

            // Extract the barcodes from the read structure and build up the record to be output
            for (idx, (read, structure)) in
                zipped_reads.iter().zip(self.grouped_read_segments.iter()).enumerate()
            {
                for (kind, segments) in &structure.segments_by_kind {
                    let (read, counter) = self.process_segments(
                        *kind,
                        read,
                        segments,
                        self.read_filter_config.quality_mask_thresholds[idx],
                        &mut sample_barcode,
                        &mut umi_barcode,
                    )?;
                    if let Some(read) = read {
                        output_reads.push(read);
                    }

                    qual_counter.update_with_self(&counter);
                }
            }

            let num_no_calls = sample_barcode.num_no_calls();
            let (concat, delim) = sample_barcode.into_inner();

            // Determine which sample (if any) this barcode belongs to
            let match_result = match self.read_filter_config.max_no_calls {
                Some(max_no_calls) if num_no_calls > max_no_calls => {
                    MatchResult::NoMatch { barcode: concat }
                }
                _ => self.matcher.find(concat),
            };

            let sample_idx = match match_result {
                MatchResult::Match { sample_index, .. } => sample_index,
                MatchResult::NoMatch { .. } => self.undetermined_sample_index,
            };

            let sample_metrics = &mut demuxed_group.group_metrics.per_sample_metrics[sample_idx];
            sample_metrics.update_with_match(&match_result);

            // Handle optional metrics tracking that occurs on unmatched barcodes
            if let MatchResult::NoMatch { barcode } = match_result {
                // If sample_barcode_hop checking is in effect, track this barcode
                if let Some(ref mut tracker) =
                    demuxed_group.group_metrics.sample_barcode_hop_tracker
                {
                    tracker.check_barcode(&barcode);
                }
                // If collection of top unmatched barcodes is in effect, track this barcode
                if let Some(ref mut unmatched) = demuxed_group.unmatched {
                    unmatched.push(barcode);
                }
            }

            // Update the header information for the record and set it in its corresponding sample + fastq slot in the output
            for (i, mut read) in output_reads.into_iter().enumerate() {
                let mut header =
                    FastqHeader::try_from(read.head.as_slice()).with_context(|| {
                        format!(
                            "Unable to parse read header: {}",
                            String::from_utf8_lossy(&read.head),
                        )
                    })?;
                if !delim.is_empty() {
                    if let Err(e) = header.set_sample_barcode(&delim) {
                        debug!("{:#}", e);
                    }
                }
                if !umi_barcode.is_empty() {
                    header.set_umi(&umi_barcode.delimited);
                }
                let mut updated_header = Vec::with_capacity(read.head.len());
                header.copy_to_vec(&mut updated_header);
                read.head = updated_header;
                demuxed_group.per_sample_reads[sample_idx].per_fastq_reads[i].push(read);
            }

            demuxed_group.group_metrics.total_templates += 1;
            let sample_metrics = &mut demuxed_group.group_metrics.per_sample_metrics[sample_idx];
            sample_metrics.base_qual_counter.update_with_self(&qual_counter);
        }

        Ok(demuxed_group)
    }

    fn demux_hop_checker(&self) -> Option<&SampleBarcodeHopChecker> {
        self.sample_barcode_hop_checker.as_ref()
    }

    #[inline]
    fn get_undetermined_index(&self) -> usize {
        self.undetermined_sample_index
    }

    #[inline]
    fn samples(&self) -> &[SampleMetadata] {
        self.samples
    }
}

/// Helper struct that has a vector where each element corresponds to an output FASTQ file.
/// Each inner vector is a vector of [`OwnedRecord`]s to write to the specified FASTQ file.
#[derive(Debug)]
pub struct OutputPerSampleReads {
    /// A vec of vecs to hold reads that go to each output FASTQ file.
    pub per_fastq_reads: Vec<Vec<OwnedRecord>>,
}

impl OutputPerSampleReads {
    /// Create a new [`OutputPerSampleReads`] with space preallocated for each of the fastqs.
    pub fn new(number_of_fastqs_for_sample: usize) -> Self {
        Self { per_fastq_reads: vec![vec![]; number_of_fastqs_for_sample] }
    }

    /// Check if all `per_fastq_reads` are empty.
    pub fn is_empty(&self) -> bool {
        self.per_fastq_reads.iter().all(Vec::is_empty)
    }

    /// Get the number of reads that have been added to the `per_fastq_reads` vec.
    pub fn len(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            self.per_fastq_reads[0].len()
        }
    }
}

/// A collection of demultiplexed reads across samples.
///
/// The `sample_records` holds a vector where each element is a collection of reads for a given sample.
#[derive(Debug)]
pub struct DemuxedGroup<'a> {
    /// The [`OutputPerSampleReads`] for each sample.
    pub per_sample_reads: Vec<OutputPerSampleReads>,
    /// The metrics for this [`DemuxedGroup`]. This includes per-sample metrics.
    pub group_metrics: DemuxedGroupMetrics<'a>,
    /// The optional unmatched barcodes.
    pub unmatched: Option<Vec<Vec<u8>>>,
}

impl<'a> DemuxedGroup<'a> {
    /// Create a new [`DemuxedGroup`] with space pre-allocated based on the number of samples and number of output FASTQs per sample.
    pub fn new(
        num_samples: usize,
        num_output_fastqs: usize,
        sample_barcode_hop_checker: Option<&'a SampleBarcodeHopChecker>,
        collect_unmatched: bool,
    ) -> Self {
        let mut per_sample_reads = Vec::with_capacity(num_samples);
        for _ in 0..num_samples {
            per_sample_reads.push(OutputPerSampleReads::new(num_output_fastqs));
        }
        Self {
            per_sample_reads,
            group_metrics: DemuxedGroupMetrics::with_capacity(
                sample_barcode_hop_checker,
                num_samples,
            ),
            unmatched: if collect_unmatched { Some(vec![]) } else { None },
        }
    }
}

#[cfg(test)]
#[allow(clippy::all)]
mod test {
    use bstr::BString;
    use read_structure::{ReadStructure, SegmentType};
    use rstest::rstest;
    use seq_io::fastq::OwnedRecord;
    use std::str::FromStr;

    use crate::{
        demux::{self, Demultiplex, DemuxReadFilterConfig},
        fastq_header::FastqHeader,
        matcher::{
            CachedHammingDistanceMatcher, MatcherKind, PreComputeMatcher, UNDETERMINED_NAME,
        },
        opts::Opts,
        sample_metadata::SampleMetadata,
        utils::test_commons::{
            create_preset_sample_metadata, reads_to_record_set, Fq, SAMPLE_BARCODE_1,
            SAMPLE_BARCODE_4,
        },
    };

    use super::{Demultiplexer, PerFastqRecordSet};

    struct DemuxTestContext {
        read_filter: DemuxReadFilterConfig,
        opts: Opts,
        output_types: Vec<SegmentType>,
        metadata: Vec<SampleMetadata>,
        per_fastq_record_set: PerFastqRecordSet,
        max_mismatches: usize,
        min_delta: usize,
        matcher: MatcherKind,
    }

    impl DemuxTestContext {
        fn new(
            reads: Vec<OwnedRecord>,
            read_structures: Vec<ReadStructure>,
            max_mismatches: usize,
            min_delta: usize,
            max_no_calls: Option<usize>,
            filter_control_reads: bool,
            filter_failing_quality: bool,
            quality_mask_threshold: u8,
            matcher: MatcherKind,
        ) -> Self {
            let opts = Opts {
                read_structures,
                allowed_mismatches: max_mismatches,
                min_delta,
                max_no_calls,
                filter_failing_quality,
                filter_control_reads,
                quality_mask_threshold: vec![quality_mask_threshold],
                ..Opts::default()
            };

            let read_filter = opts.as_read_filter_config();
            let output_types = opts.output_types_to_write().unwrap();
            let metadata = create_preset_sample_metadata();
            let reads = reads.into_iter().map(|r| reads_to_record_set([r].into_iter())).collect();
            let per_fastq_record_set = PerFastqRecordSet::new(reads).unwrap();
            Self {
                read_filter,
                opts,
                output_types,
                metadata,
                per_fastq_record_set,
                max_mismatches,
                min_delta,
                matcher,
            }
        }

        fn make_demuxer<'a>(&'a self) -> Box<dyn Demultiplex + 'a> {
            match self.matcher {
                MatcherKind::CachedHammingDistance => {
                    let matcher = CachedHammingDistanceMatcher::new(
                        &self.metadata[0..self.metadata.len() - 1],
                        self.max_mismatches,
                        self.min_delta,
                        2,
                    );
                    Box::new(
                        Demultiplexer::new(
                            &self.metadata,
                            &self.opts.read_structures,
                            &self.output_types,
                            &self.read_filter,
                            matcher,
                            false,
                            false,
                            self.opts.sample_barcode_in_fastq_header,
                        )
                        .unwrap(),
                    )
                }
                MatcherKind::PreCompute => {
                    let matcher = PreComputeMatcher::new(
                        &self.metadata[0..self.metadata.len() - 1],
                        self.max_mismatches,
                        self.min_delta,
                        2,
                    );
                    Box::new(
                        Demultiplexer::new(
                            &self.metadata,
                            &self.opts.read_structures,
                            &self.output_types,
                            &self.read_filter,
                            matcher,
                            false,
                            false,
                            self.opts.sample_barcode_in_fastq_header,
                        )
                        .unwrap(),
                    )
                }
            }
        }

        fn demux(&self) -> Box<dyn Demultiplex + '_> {
            self.make_demuxer()
        }

        /// Helper function for setting up all the inputs to demux to test different read structure inputs
        fn demux_structures(
            reads: Vec<OwnedRecord>,
            read_structures: Vec<ReadStructure>,
            matcher: MatcherKind,
        ) -> Self {
            Self::new(reads, read_structures, 2, 1, Some(1), false, false, 0, matcher)
        }

        /// Helper function for setting up all the inputs to demux to test different mismatch combos
        fn demux_mismatches(
            reads: Vec<OwnedRecord>,
            read_structures: Vec<ReadStructure>,
            max_mismatches: usize,
            min_delta: usize,
            max_no_calls: Option<usize>,
            matcher: MatcherKind,
        ) -> Self {
            Self::new(
                reads,
                read_structures,
                max_mismatches,
                min_delta,
                max_no_calls,
                false,
                false,
                0,
                matcher,
            )
        }

        /// Helper function for setting up all the inputs to demux to test read filtering and masking
        fn demux_filtering(
            reads: Vec<OwnedRecord>,
            read_structures: Vec<ReadStructure>,
            filter_control_reads: bool,
            filter_failing_quality: bool,
            matcher: MatcherKind,
        ) -> Self {
            Self::new(
                reads,
                read_structures,
                0,
                1,
                Some(1),
                filter_control_reads,
                filter_failing_quality,
                0,
                matcher,
            )
        }

        /// Helper function for setting up all the inputs to demux to test read filtering and masking
        fn demux_masking(
            reads: Vec<OwnedRecord>,
            read_structures: Vec<ReadStructure>,
            quality_mask_threshold: u8,
            matcher: MatcherKind,
        ) -> Self {
            Self::new(
                reads,
                read_structures,
                0,
                1,
                Some(1),
                false,
                false,
                quality_mask_threshold,
                matcher,
            )
        }
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_fragment_reads(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![ReadStructure::from_str("17B100T").unwrap()];
        let fastq_record = Fq { name: "frag", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_structures(vec![fastq_record], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();

        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty),  "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                      1,            "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),      1,            "Sample1 has only a single output fastq expected.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].seq,  &[b'A'; 100], "Sample1 seq has had the barcode removed");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].qual, &[b'!'; 100], "Sample1 quals has had the barcode removed");

        let header = FastqHeader::try_from(result.per_sample_reads[0].per_fastq_reads[0][0].head.as_slice()).unwrap();

        assert!(header.umi().is_none(),                            "Fastq header has no UMI set");
        assert_eq!(header.sample_barcode(), Some(SAMPLE_BARCODE_1),"Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_paired_reads_inline_sample_barcodes(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("8B100T").unwrap(),
            ReadStructure::from_str("9B100T").unwrap(),
        ];

        let fq1 = Fq { name: "frag", bases: &[&SAMPLE_BARCODE_1[0..8], &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();
        let fq2 = Fq { name: "frag", bases: &[&SAMPLE_BARCODE_1[8..],  &[b'T'; 100]].concat(), ..Fq::default() }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1, fq2], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                                1, "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),                2, "Sample1 has 2 output fastq expected.");

        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];
        let r2 = &result.per_sample_reads[0].per_fastq_reads[1][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has had the barcode removed");
        assert_eq!(r2.seq,  &[b'T'; 100], "Sample1 R2 seq has had the barcode removed");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has had the barcode removed");
        assert_eq!(r2.qual, &[b'!'; 100], "Sample1 R2 quals has had the barcode removed");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();
        let r2_header = FastqHeader::try_from(r2.head.as_slice()).unwrap();

        assert!(r1_header.umi().is_none(), "R1 Fastq header has no UMI set");
        assert!(r2_header.umi().is_none(), "R2 Fastq header has no UMI set");

        let expected = [&SAMPLE_BARCODE_1[0..8], &[b'+'], &SAMPLE_BARCODE_1[8..]].concat();
        assert_eq!(r1_header.sample_barcode(), Some(expected.as_slice()), "R1 Set barcode matches Sample1");
        assert_eq!(r2_header.sample_barcode(), Some(expected.as_slice()), "R2 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_dual_indexed_paired_end_reads(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("8B").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
            ReadStructure::from_str("9B").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
        ];

        let fq1 = Fq { name: "pair", bases: &SAMPLE_BARCODE_1[0..8], ..Fq::default() }.to_owned_record();
        let fq2 = Fq { name: "pair", bases: &[b'A'; 100],            ..Fq::default() }.to_owned_record();
        let fq3 = Fq { name: "pair", bases: &SAMPLE_BARCODE_1[8..],  ..Fq::default() }.to_owned_record();
        let fq4 = Fq { name: "pair", bases: &[b'T'; 100],            ..Fq::default() }.to_owned_record();

        let  context = DemuxTestContext::demux_structures(vec![fq1, fq2, fq3, fq4], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                                1, "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),                2, "Sample1 has 2 output fastq expected.");

        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];
        let r2 = &result.per_sample_reads[0].per_fastq_reads[1][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r2.seq,  &[b'T'; 100], "Sample1 R2 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");
        assert_eq!(r2.qual, &[b'!'; 100], "Sample1 R2 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();
        let r2_header = FastqHeader::try_from(r2.head.as_slice()).unwrap();

        assert!(r1_header.umi().is_none(), "R1 Fastq header has no UMI set");
        assert!(r2_header.umi().is_none(), "R2 Fastq header has no UMI set");

        let expected = [&SAMPLE_BARCODE_1[0..8], &[b'+'], &SAMPLE_BARCODE_1[8..]].concat();
        assert_eq!(r1_header.sample_barcode(), Some(expected.as_slice()), "R1 Set barcode matches Sample1");
        assert_eq!(r2_header.sample_barcode(), Some(expected.as_slice()), "R2 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_very_weird_set_of_reads(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("4B4M8S").unwrap(),
            ReadStructure::from_str("4B100T").unwrap(),
            ReadStructure::from_str("100S3B").unwrap(),
            ReadStructure::from_str("6B1S1M1T").unwrap(),
        ];

        let fq1 = Fq { name: "pair", bases: &[&SAMPLE_BARCODE_1[0..4], b"CCCCGGGGTTTT"].concat(), ..Fq::default() }.to_owned_record();
        let fq2 = Fq { name: "pair", bases: &[&SAMPLE_BARCODE_1[4..8], &[b'A'; 100]].concat(),    ..Fq::default() }.to_owned_record();
        let fq3 = Fq { name: "pair", bases: &[&[b'T'; 100], &SAMPLE_BARCODE_1[8..11]].concat(),   ..Fq::default() }.to_owned_record();
        let fq4 = Fq { name: "pair", bases: &[&SAMPLE_BARCODE_1[11..], b"AAT"].concat(),          ..Fq::default() }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1, fq2, fq3, fq4], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                                1, "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),                2, "Sample1 has 2 output fastq expected.");

        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];
        let r2 = &result.per_sample_reads[0].per_fastq_reads[1][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r2.seq,  &[b'T'; 1],   "Sample1 R2 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");
        assert_eq!(r2.qual, &[b'!'; 1],   "Sample1 R2 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();
        let r2_header = FastqHeader::try_from(r2.head.as_slice()).unwrap();

        let expected_umi =  b"CCCC+A";
        assert_eq!(r1_header.umi(), Some(expected_umi.as_slice()), "R1 Fastq header has UMI set");
        assert_eq!(r2_header.umi(), Some(expected_umi.as_slice()), "R2 Fastq header has UMI set");

        let expected = [&SAMPLE_BARCODE_1[0..4], &[b'+'], &SAMPLE_BARCODE_1[4..8], &[b'+'], &SAMPLE_BARCODE_1[8..11], &[b'+'], &SAMPLE_BARCODE_1[11..]].concat();
        assert_eq!(r1_header.sample_barcode(), Some(expected.as_slice()), "R1 Set barcode matches Sample1");
        assert_eq!(r2_header.sample_barcode(), Some(expected.as_slice()), "R2 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_a_read_structure_with_multiple_independent_template_segments_in_one_read(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B20T20S20T20S20T").unwrap(),
        ];

        let fq1 = Fq {
            name: "frag", bases: &[SAMPLE_BARCODE_1, &[b'A'; 20], &[b'C'; 20], &[b'A'; 20], &[b'C'; 20], &[b'A'; 20]].concat(), ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                                1, "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),                1, "Sample1 has 2 output fastq expected.");

        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];

        assert_eq!(r1.seq,  &[b'A'; 60], "Sample1 R1 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 60], "Sample1 R1 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();

        eprintln!("{}", String::from_utf8_lossy(&r1_header.sample_barcode().unwrap()));
        assert!(r1_header.umi().is_none(),                             "R1 Fastq header has no UMI set");
        assert_eq!(r1_header.sample_barcode(), Some(SAMPLE_BARCODE_1), "R1 Set barcode matches Sample1");
    }

    #[rstest]
    #[should_panic]
    #[rustfmt::skip]
    fn test_demux_fails_if_no_read_structures(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![];

        let fq1 = Fq {
            name: "frag", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1], read_structures, matcher);
        let demuxer = context.demux();
        let _result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();

    }

    #[rstest]
    #[should_panic]
    #[rustfmt::skip]
    fn test_demux_fails_if_no_template_sequence_found(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("100S").unwrap(),
            ReadStructure::from_str("100S").unwrap()
        ];

        let fq1 = Fq {
            name: "frag", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1], read_structures, matcher);
        let demuxer = context.demux();
        let _result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
    }

    #[rstest]
    #[should_panic]
    #[rustfmt::skip]
    fn test_demux_fails_if_too_many_fastq_records_passed_in(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {

        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let fq1 = Fq {
            name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();
        let fq2 = Fq {
            name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();
        let fq3 = Fq {
            name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1, fq2, fq3], read_structures, matcher);
        let demuxer = context.demux();
        let _result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
    }

    #[rstest]
    #[should_panic]
    #[rustfmt::skip]
    fn test_demux_fails_if_too_few_fastq_records_passed_in(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
            ReadStructure::from_str("17B100T").unwrap(),
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let fq1 = Fq {
            name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();
        let fq2 = Fq {
            name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(), ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1, fq2], read_structures, matcher);
        let demuxer = context.demux();
        let _result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_molecular_barcode_should_be_set(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B5M100T").unwrap(),
        ];

        let fq1 = Fq { name: "pair", bases: &[SAMPLE_BARCODE_1, &[b'N'; 5], &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();

        let context = DemuxTestContext::demux_structures(vec![fq1], read_structures, matcher);
                let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                                1, "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),                1, "Sample1 has 1 output fastq expected.");

        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();

        let expected_umi =  b"NNNNN";
        assert_eq!(r1_header.umi(), Some(expected_umi.as_slice()), "R1 Fastq header has UMI set");

        assert_eq!(r1_header.sample_barcode(), Some(SAMPLE_BARCODE_1), "R1 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_assigns_to_unmatched_if_read_matches_to_sample_barcodes_within_mismatch_delta() {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let fq1 = Fq { name: "pair", bases: &[SAMPLE_BARCODE_4, &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();

        let context = DemuxTestContext::demux_mismatches(vec![fq1], read_structures, 2, 3, Some(1), MatcherKind::CachedHammingDistance);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[0..4].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-4 have no reads assigned");
        assert_eq!(result.per_sample_reads[4].len(),                                1, "UNDETERMINED contains a single record.");
        assert_eq!(result.per_sample_reads[4].per_fastq_reads.len(),                1, "UNDETERMINED has 1 output fastq expected.");

        let r1 = &result.per_sample_reads[4].per_fastq_reads[0][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();

        assert!(r1_header.umi().is_none(), "R1 Fastq header has no UMI set");
        assert_eq!(r1_header.sample_barcode(), Some(SAMPLE_BARCODE_4), "R1 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_assignes_to_unmatched_if_too_many_mismatches(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let sample1_with_mismatch = b"AAAAAAAAGATTACAGT"; // last two barcode bases have mismatchs
        let fq1 = Fq { name: "pair", bases: &[sample1_with_mismatch.as_slice(), &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();

        let context = DemuxTestContext::demux_mismatches(vec![fq1], read_structures, 0, 1, Some(1), matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[0..4].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-4 have no reads assigned");
        assert_eq!(result.per_sample_reads[4].len(),                                1, "UNDETERMINED contains a single record.");
        assert_eq!(result.per_sample_reads[4].per_fastq_reads.len(),                1, "UNDETERMINED has 1 output fastq expected.");

        let r1 = &result.per_sample_reads[4].per_fastq_reads[0][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();

        assert!(r1_header.umi().is_none(), "R1 Fastq header has no UMI set");
        assert_eq!(r1_header.sample_barcode(), Some(sample1_with_mismatch.as_slice()), "R1 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_assigns_to_unmatched_if_reads_samp_barcode_has_too_many_ns(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let mut barcode = SAMPLE_BARCODE_1.to_vec();
        barcode[0] = b'N';

        let fq1 = Fq { name: "pair", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_mismatches(vec![fq1], read_structures, 1, 1, Some(0), matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads[0..4].iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-4 have no reads assigned");
        assert_eq!(result.per_sample_reads[4].len(),                                1, "UNDETERMINED contains a single record.");
        assert_eq!(result.per_sample_reads[4].per_fastq_reads.len(),                1, "UNDETERMINED has 1 output fastq expected.");

        let r1 = &result.per_sample_reads[4].per_fastq_reads[0][0];

        assert_eq!(r1.seq,  &[b'A'; 100], "Sample1 R1 seq has expected seq");
        assert_eq!(r1.qual, &[b'!'; 100], "Sample1 R1 quals has expected qual");

        let r1_header = FastqHeader::try_from(r1.head.as_slice()).unwrap();

        assert!(r1_header.umi().is_none(), "R1 Fastq header has no UMI set");
        assert_eq!(r1_header.sample_barcode(), Some(barcode.as_slice()), "R1 Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_failed_true_and_quality_y(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:Y:0:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, true, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_failed_true_and_quality_n(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:0:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, true, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_failed_false_and_quality_y(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:Y:0:1"),..Fq::default() }.to_owned_record();
        let context= DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_failed_false_and_quality_n(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:0:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_control_true_control_field_0(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:0:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, true, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_control_true_control_field_1(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:1:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, true, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_control_false_control_field_0(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:0:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_filter_control_false_control_field_1(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:1:1"),..Fq::default() }.to_owned_record();
        let context = DemuxTestContext::demux_filtering(vec![fq1], read_structures, false, false, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty), "Samples 1-UNDETERMINED have no reads assigned");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_combinations_of_filters(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];
        let barcode = SAMPLE_BARCODE_1.to_vec();

        { // pass QC and is not an internal control
            let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:0:1"),..Fq::default() }.to_owned_record();

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));
        }

        { // Does not pass QC and is not an internal control
            let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:Y:0:1"),..Fq::default() }.to_owned_record();

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));
        }

        { // Passes QC but is an internal control
            let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:N:1:1"),..Fq::default() }.to_owned_record();

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));
    }

    { // Does not pass QC and is an internal control
            let fq1 = Fq { name: "frag", bases: &[barcode.as_slice(), &[b'A'; 100]].concat(), comment: Some("1:Y:1:1"),..Fq::default() }.to_owned_record();

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(!result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, false, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), false, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));

            let context = DemuxTestContext::demux_filtering(vec![fq1.clone()], read_structures.clone(), true, true, matcher);
            let demuxer = context.demux();
            let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
            assert!(result.per_sample_reads.iter().all(demux::OutputPerSampleReads::is_empty));
        }
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_should_mask_bases_less_than_the_quality_threshold(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq {
            name: "frag",
            bases: &[barcode.as_slice(), &[b'A'; 100]].concat(),
            quals: Some(&[[b'?'; 100].as_slice(), [b'!'; 17].as_slice()].concat()),
            ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_masking(vec![fq1], read_structures, 20, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];
        assert_eq!(&r1.seq, &[[b'A'; 83].as_slice(), [b'N'; 17].as_slice()].concat());
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_should_not_mask_bases_when_threshold_is_zero(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
        ];

        let barcode = SAMPLE_BARCODE_1.to_vec();

        let fq1 = Fq {
            name: "frag",
            bases: &[barcode.as_slice(), &[b'A'; 100]].concat(),
            quals: Some(&[[b'?'; 100].as_slice(), [b'!'; 17].as_slice()].concat()),
            ..Fq::default()
        }.to_owned_record();

        let context = DemuxTestContext::demux_masking(vec![fq1], read_structures, 0, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();
        let r1 = &result.per_sample_reads[0].per_fastq_reads[0][0];
        assert_eq!(&r1.seq, &[[b'A'; 83].as_slice(), [b'A'; 17].as_slice()].concat());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn test_demux_simple_multi_read() {
        let num_reads = 1000;

        let r1_reads = vec![
            vec![
                Fq {
                    name: "Sample1R1",
                    bases: b"AAAAAAAAAA",
                    set_bc_to: Some(b"ACTGACTG"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
            vec![
                Fq {
                    name: "Sample2R1",
                    bases: b"TTTTTTTTTT",
                    set_bc_to: Some(b"GTCAGTCA"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
        ];

        let r2_reads = vec![
            vec![
                Fq {
                    name: "Sample1R2",
                    bases: b"AAAAAAAAAA",
                    set_bc_to: Some(b"ACTGACTG"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
            vec![
                Fq {
                    name: "Sample2R2",
                    bases: b"TTTTTTTTTT",
                    set_bc_to: Some(b"GTCAGTCA"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
        ];
        let i1_reads = vec![
            vec![
                Fq {
                    name: "Sample1I1",
                    bases: b"ACTGACTG",
                    set_bc_to: Some(b"ACTGACTG"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
            vec![
                Fq {
                    name: "Sample2I1",
                    bases: b"GTCAGTCA",
                    set_bc_to: Some(b"GTCAGTCA"),
                    ..Fq::default()
                }
                .to_owned_record();
                num_reads
            ],
        ];

        let mut metadata = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTGACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("GTCAGTCA"), 1, 3).unwrap(),
        ];

        let output_types = vec![SegmentType::Template];

        let read_structures: Vec<ReadStructure> = vec!["+T", "+T", "8B"]
            .into_iter()
            .map(|s| ReadStructure::from_str(s).unwrap())
            .collect();

        metadata.push(
            SampleMetadata::new_allow_invalid_bases(
                String::from(UNDETERMINED_NAME),
                BString::from(vec![b'N'; metadata[0].barcode.len()]),
                metadata.len(),
            )
            .unwrap(),
        );

        let r1_record_set = reads_to_record_set(r1_reads.iter().cloned().flatten());
        let r2_record_set = reads_to_record_set(r2_reads.iter().cloned().flatten());
        let i1_record_set = reads_to_record_set(i1_reads.iter().cloned().flatten());

        let filter_config = DemuxReadFilterConfig {
            quality_mask_thresholds: vec![0, 0, 0],
            ..DemuxReadFilterConfig::default()
        };

        let matcher = PreComputeMatcher::new(&metadata, 1, 0, 2);
        let demuxer = Demultiplexer::new(
            &metadata,
            &read_structures,
            &output_types,
            &filter_config,
            matcher,
            false,
            false,
            false,
        )
        .unwrap();
        let per_fastq_record_set =
            PerFastqRecordSet::new(vec![r1_record_set, r2_record_set, i1_record_set]).unwrap();

        let demux_result = demuxer.demultiplex(&per_fastq_record_set).unwrap();

        assert!(demux_result.per_sample_reads[2].is_empty(), "UNDETERMINED should be empty");
        assert_eq!(
            demux_result.per_sample_reads[0].len(),
            num_reads,
            "Sample1 should have the expected number of reads"
        );
        assert_eq!(
            demux_result.per_sample_reads[1].len(),
            num_reads,
            "Sample2 should have the expected number of reads"
        );

        assert_eq!(
            demux_result.per_sample_reads[0].per_fastq_reads[0], r1_reads[0],
            "Sample1 R1 reads should match expected"
        );
        assert_eq!(
            demux_result.per_sample_reads[0].per_fastq_reads[1], r2_reads[0],
            "Sample1 R2 reads should match expected"
        );
        assert_eq!(
            demux_result.per_sample_reads[1].per_fastq_reads[0], r1_reads[1],
            "Sample2 R1 reads should match expected"
        );
        assert_eq!(
            demux_result.per_sample_reads[1].per_fastq_reads[1], r2_reads[1],
            "Sample2 R2 reads should match expected"
        );
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_paired_reads_different_read_names(
        #[values(MatcherKind::CachedHammingDistance, MatcherKind::PreCompute)]
        matcher: MatcherKind
    ) {
        let read_structures = vec![
            ReadStructure::from_str("8B100T").unwrap(),
            ReadStructure::from_str("9B100T").unwrap(),
        ];

        let mut fq1: OwnedRecord = Fq { name: "frag1", bases: &[&SAMPLE_BARCODE_1[0..8], &[b'A'; 100]].concat(), ..Fq::default() }.to_owned_record();
        let fq2 = Fq { name: "frag2", bases: &[&SAMPLE_BARCODE_1[8..],  &[b'T'; 100]].concat(), ..Fq::default() }.to_owned_record();

        // mess with the FQ1 read name
        fq1.head = vec![b'f'];
        fq1.head.append(&mut fq2.head.clone());

        let context = DemuxTestContext::demux_structures(vec![fq1, fq2], read_structures, matcher);
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set);
        assert!(result.is_err());
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_fragment_reads_sample_barcode_from_fastq_header_single_index() {
        let read_structures = vec![ReadStructure::from_str("100T").unwrap()];
        let fastq_record = Fq { name: "frag",  bases: &[b'A'; 100], set_bc_to: Some(SAMPLE_BARCODE_1), ..Fq::default() }.to_owned_record();
        let mut context = DemuxTestContext::demux_structures(vec![fastq_record], read_structures, MatcherKind::CachedHammingDistance);
        context.opts.sample_barcode_in_fastq_header = true;
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();

        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty),  "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                      1,            "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),      1,            "Sample1 has only a single output fastq expected.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].seq,  &[b'A'; 100], "Sample1 seq has had the barcode removed");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].qual, &[b'!'; 100], "Sample1 quals has had the barcode removed");

        let header = FastqHeader::try_from(result.per_sample_reads[0].per_fastq_reads[0][0].head.as_slice()).unwrap();

        assert!(header.umi().is_none(),                            "Fastq header has no UMI set");
        assert_eq!(header.sample_barcode(), Some(SAMPLE_BARCODE_1),"Set barcode matches Sample1");
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_demux_fragment_reads_sample_barcode_from_fastq_header_dual_index() {
        let read_structures = vec![ReadStructure::from_str("100T").unwrap()];
        let (index1, index2) = SAMPLE_BARCODE_1.split_at(10);
        let delimited_sample_barcode = [index1, &[b'+'], index2].concat();
        let fastq_record = Fq { name: "frag",  bases: &[b'A'; 100], set_bc_to: Some(&delimited_sample_barcode), ..Fq::default() }.to_owned_record();
        let mut context = DemuxTestContext::demux_structures(vec![fastq_record], read_structures, MatcherKind::CachedHammingDistance);
        context.opts.sample_barcode_in_fastq_header = true;
        let demuxer = context.demux();
        let result = demuxer.demultiplex(&context.per_fastq_record_set).unwrap();

        assert!(result.per_sample_reads[1..].iter().all(demux::OutputPerSampleReads::is_empty),  "Samples 2-4 and UNDETERMINED have no reads assigned");
        assert_eq!(result.per_sample_reads[0].len(),                      1,            "Sample1 contains a single record.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads.len(),      1,            "Sample1 has only a single output fastq expected.");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].seq,  &[b'A'; 100], "Sample1 seq has had the barcode removed");
        assert_eq!(result.per_sample_reads[0].per_fastq_reads[0][0].qual, &[b'!'; 100], "Sample1 quals has had the barcode removed");

        let header = FastqHeader::try_from(result.per_sample_reads[0].per_fastq_reads[0][0].head.as_slice()).unwrap();
        assert!(header.umi().is_none(),                            "Fastq header has no UMI set");
        assert_eq!(header.sample_barcode(), Some(delimited_sample_barcode.as_slice()),"Set barcode matches Sample1");
    }
}
