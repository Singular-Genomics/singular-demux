//! Functionality pertaining to the collection of metrics during demultiplexing.
//!
//! For each [`crate::demux::DemuxedGroup`] that is created in demultiplexing, a corresponding [`DemuxedGroupMetrics`]
//! is also created. The [`DemuxedGroupMetrics`] object collects metrics for that group of demuxed FASTQ records
//! and has an [`DemuxedGroupMetrics::update_with`] method that allows group metrics to be aggregated.
//!
//! Each group metrics holds a count of the number of control reads seen, number of failed filter reads seen, total
//! templates seen, as well as various per sample metrics and an optional [`SampleBarcodeHopTracker`].
//!
//! The per-sample metrics [`DemuxedGroupSampleMetrics`] keep track of the number of perfect, off by one, and total matches
//! for the given sample. Additionally they also track the quality metrics via a [`BaseQualCounter`].
//!
//! The optional [`SampleBarcodeHopTracker`] is used only in the case where the input read structure indicates dual-index reads.
//! Any barcodes that don't match to any sample are instead tested for index hopping and then tracked by the tracker.
//!
//! Lastly, an `UnmatchedCounterThread` can be created to keep track of the most frequently seen unmatched barcodes.
//! Unmatched barcodes from a [`crate::demux::DemuxedGroup`] are collected and sent over a channel to the unmatched
//! counter thread which adds the barcodes to hash. Each time the number of keys in that hash exceeds a set limit
//! the unmatched barcodes are sorted by most to least frequent, taking only the top N barcodes and and dropping
//! the rest.
//!
//! All metrics are writable to files.

use std::thread::JoinHandle;
use std::{collections::HashMap, path::Path};

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use bstr::BString;
use fgoxide::io::DelimFile;
use flume::{Receiver, Sender};
use itertools::Itertools;
use read_structure::{ReadSegment, ReadStructure, SegmentType};
use serde::{Deserialize, Serialize};

use crate::{matcher::MatchResult, sample_metadata::SampleMetadata};

/// The max number of barcode sets that can be in the channel at one time
pub const DEFAULT_UNMATCHED_CHANNEL_SIZE: usize = 1_000;
/// The max number of keys the hash can retain before re-sizing
pub const DEFAULT_UNMATCHED_MAX_COUNTER_SIZE: usize = 5_000_000;
/// The number of keys to downsize to
pub const DEFAULT_UNMATCHED_DOWNSIZE_TO: usize = 5_000;

/// Helper type to reduce type complexity in channels.
type RawBarcode = Vec<u8>;

/// The [`UnmatchedCounterThread`] tracks the most frequently seen unmatched barcodes.
///
/// This is done with a heuristic. Every time the inner hashmap grows over a set number of keys
/// the hashmap is sorted from most to least frequent and only the top N keys are kept, dropping the
/// remainder. This prevents the hashmap from growing infinitely large.
pub struct UnmatchedCounterThread {
    /// The join handle for the thread the [`UnmatchedCounterThread`] is running on.
    pub handle: Option<JoinHandle<UnmatchedCounter>>,
    /// The `Sender` end of the channel to send barcodes to the counter.
    pub tx: Option<Sender<Vec<RawBarcode>>>,
    /// A oneshot channel to telling the counter thread to shutdown.
    shutdown_tx: Option<Sender<()>>,
}

impl UnmatchedCounterThread {
    /// Create a new [`UnmatchedCounterThread`].
    ///
    /// This spins up the thread in the background.
    ///
    /// # Arguments
    /// - `channel_size` - The max number of barcode sets that can be in the channel at one time, defaults to [`DEFAULT_UNMATCHED_CHANNEL_SIZE`].
    /// - `max_counter_size` - The max number of keys the hash can retain before re-sizing, defaults to [`DEFAULT_UNMATCHED_MAX_COUNTER_SIZE`].
    /// - `downsize_to` - The number of keys to downsize to, defaults to [`DEFAULT_UNMATCHED_DOWNSIZE_TO`].
    pub fn new(
        channel_size: Option<usize>,
        max_counter_size: Option<usize>,
        downsize_to: Option<usize>,
    ) -> Result<Self> {
        let (tx, rx): (Sender<Vec<RawBarcode>>, Receiver<Vec<RawBarcode>>) =
            flume::bounded(channel_size.unwrap_or(DEFAULT_UNMATCHED_CHANNEL_SIZE));
        let (shutdown_tx, shutdown_rx) = flume::unbounded();
        let handle = std::thread::spawn(move || {
            let mut collector = UnmatchedCounter::new(
                max_counter_size.unwrap_or(DEFAULT_UNMATCHED_MAX_COUNTER_SIZE),
                downsize_to.unwrap_or(DEFAULT_UNMATCHED_DOWNSIZE_TO),
            );
            while let Ok(barcodes) = rx.recv() {
                Self::process_barcodes(barcodes, &mut collector);

                // Shutdown has been requested. Drain anything currently in the queue and drop the
                // receiver. Any further attempts to send will fail.
                if shutdown_rx.is_disconnected() {
                    // Drain anything left in the queue and shutdown
                    let remainder = rx.drain().collect::<Vec<_>>();
                    drop(rx);
                    for barcodes in remainder {
                        Self::process_barcodes(barcodes, &mut collector);
                    }
                    break;
                }
            }
            collector
        });
        Ok(Self { handle: Some(handle), tx: Some(tx), shutdown_tx: Some(shutdown_tx) })
    }

    /// Process each barcode from the collection of barcodes that has been sent over the channel.
    fn process_barcodes(barcodes: Vec<RawBarcode>, collector: &mut UnmatchedCounter) {
        for barcode in barcodes {
            collector.insert(barcode);
        }
    }

    /// Stop the collector thread and return the [`UnmatchedCounter`].
    ///
    /// Any further attempts to send over the [`UnmatchedCounterThread`] channel will fail.
    pub fn collect(mut self) -> UnmatchedCounter {
        let tx = self.tx.take().unwrap();
        drop(tx);

        let shutdown = self.shutdown_tx.take().unwrap();
        drop(shutdown);

        // If the thread was panicking, resume the panic
        let counter = match self.handle.take().unwrap().join() {
            Ok(result) => result,
            Err(e) => std::panic::resume_unwind(e),
        };
        counter
    }
}

/// Container for tracking the number of times each unmatched barcode has been seen.
pub struct UnmatchedCounter {
    /// The unmatched barcode counter.
    unmatched_counter: AHashMap<RawBarcode, i64>,
    /// The max number of keys that can be held before downsizing.
    max_counter_size: usize,
    /// The number of keys to retain when downsizing.
    downsize_to: usize,
}

impl UnmatchedCounter {
    /// Create a new [`UnmatchedCounter`].
    ///
    /// # Argument
    /// - `max_counter_size` - the max number of keys the internal hashmap can contain before downsizing.
    /// - `downsize_to` - the number of keys to retain when downsizing.
    pub fn new(max_counter_size: usize, downsize_to: usize) -> Self {
        Self { unmatched_counter: AHashMap::new(), max_counter_size, downsize_to }
    }

    /// Insert a barcode in the collector and check if the criteria has been met for downsizing.
    pub fn insert(&mut self, barcode: RawBarcode) {
        // If we have hit the max_counter_size, and we are about to add one more, downsize to as to not hit
        // max_counter_size + 1
        if self.unmatched_counter.len() == self.max_counter_size {
            self.downsize();
        }
        let counter = self.unmatched_counter.entry(barcode).or_insert(0);
        *counter += 1;
    }

    /// Downsize the counter to `downsize_to` keys.
    ///
    /// Note, this is requires sorting they map by values and creation of a new hashmap.
    pub fn downsize(&mut self) {
        // TODO: would it be faster to remove the individual keys that fall below the downsize_to?
        let new = AHashMap::with_capacity(self.downsize_to);
        let previous = std::mem::replace(&mut self.unmatched_counter, new);
        self.unmatched_counter.extend(
            previous
                .into_iter()
                .sorted_unstable_by_key(|(_k, count)| -count)
                .take(self.downsize_to),
        );
    }

    /// Write the top `n` unmatched barcodes to a `most_frequent_unmatched.tsv` file in the specified directory.
    pub fn to_file<P: AsRef<Path>>(self, output_dir: P, n: usize, prefix: &str) -> Result<()> {
        let filename = [prefix.to_string(), "most_frequent_unmatched.tsv".to_string()].concat();
        let output_path = output_dir.as_ref().join(filename);
        let delim = DelimFile::default();
        delim.write_tsv(
            &output_path,
            self.unmatched_counter
                .into_iter()
                .map(|(barcode, count)| BarcodeCount::new(barcode.into(), count as isize))
                .sorted_unstable_by_key(|b| -b.count)
                .take(n),
        )?;
        Ok(())
    }
}

/// Container for aggregating metrics across [`crate::demux::DemuxedGroup`]s.
#[derive(Debug, Default)]
pub struct DemuxedGroupMetrics<'a> {
    /// The number of reads that were discarded due to having the failed quality filter set in the FASTQ header.
    pub num_reads_failing_filters: usize,
    /// The number of reads that were discarded due to being control reads.
    pub num_reads_filtered_as_control: usize,
    /// The total number of templates seen in this group.
    pub total_templates: usize,
    /// The per-sample level metrics.
    pub per_sample_metrics: Vec<DemuxedGroupSampleMetrics>,
    /// The optional [`SampleBarcodeHopTracker`].
    pub sample_barcode_hop_tracker: Option<SampleBarcodeHopTracker<'a>>,
}

/// Write the per sample metrics.
fn write_per_group_metrics<P: AsRef<Path>>(
    samples: &[SampleMetadata],
    per_group_metrics: &[DemuxedGroupSampleMetrics],
    output_dir: &P,
    filename: &str,
) -> Result<()> {
    let total_templates = per_group_metrics.iter().map(|m| m.total_matches).sum();
    // Get the group with the highest template count. Don't include unmatched when determining
    // the "best" barcode.
    let best_barcode_count = per_group_metrics[0..per_group_metrics.len() - 1]
        .iter()
        .max_by_key(|s| s.total_matches)
        .map_or(0, |s| s.total_matches);

    // Build per group metrics to output.  We need to derive (e.g. ratio, percentage) values here.
    let per_group_metrics: Vec<SampleMetricsProcessed> = samples
        .iter()
        .zip(per_group_metrics.iter())
        .map(|(sample, metrics)| metrics.as_processed(total_templates, best_barcode_count, sample))
        .collect();

    // Finally output it to a file
    let output_path = output_dir.as_ref().join(filename);
    let delim = DelimFile::default();
    delim.write_tsv(&output_path, per_group_metrics)?;

    Ok(())
}

/// Builds a new set of "samples" and demux metrics by grouping input samples by project.  Each
/// new sample corresponds to a unique project, and the associated metric are the sum of values
/// for samples with the same project.  The undetermined sample is not grouped with any other
/// sample.
fn build_per_project_metrics(
    samples: &[SampleMetadata],
    per_sample_metrics: &[DemuxedGroupSampleMetrics],
) -> Result<(Vec<SampleMetadata>, Vec<DemuxedGroupSampleMetrics>)> {
    let mut project_to_metric: AHashMap<Option<&String>, DemuxedGroupSampleMetrics> =
        AHashMap::new();
    let mut project_to_sample: HashMap<Option<&String>, SampleMetadata> = HashMap::new();
    let num_per_group_metrics = per_sample_metrics.len();

    // Iterate over all samples, except the undetermined sample
    for (sample, metric) in
        samples.iter().zip(per_sample_metrics.iter()).take(num_per_group_metrics - 1)
    {
        let key = sample.project.as_ref();
        match project_to_metric.get_mut(&key) {
            None => {
                // Build a new dummy sample for this project
                let project = sample.project.clone().unwrap_or_else(|| String::from("None"));
                let sample = SampleMetadata::new_allow_invalid_bases(
                    project,
                    BString::from(vec![b'N'; samples[0].barcode.len()]),
                    project_to_metric.len(),
                    None,
                )?;
                // Add the sample and current metric to the maps
                project_to_sample.insert(key, sample);
                project_to_metric.insert(key, metric.clone());
            }
            Some(project_metric) => project_metric.update_with(metric),
        };
    }

    let mut project_metrics: Vec<DemuxedGroupSampleMetrics> =
        Vec::with_capacity(project_to_metric.len() + 1);
    let mut project_samples: Vec<SampleMetadata> = Vec::with_capacity(project_to_metric.len() + 1);

    // Sort by project name
    for key in project_to_metric.keys().sorted() {
        let metric = project_to_metric.get(key).unwrap();
        let sample = project_to_sample.remove(key).unwrap();
        project_metrics.push(metric.clone());
        project_samples.push(sample);
    }

    // add the undetermined sample to keep it seperate
    {
        let undetermined_sample = &samples[per_sample_metrics.len() - 1];
        let undetermined_metric = &per_sample_metrics[per_sample_metrics.len() - 1];
        project_samples.push(undetermined_sample.clone());
        project_metrics.push(undetermined_metric.clone());
    }

    Ok((project_samples, project_metrics))
}

impl<'a> DemuxedGroupMetrics<'a> {
    /// Create a [`DemuxedGroupMetrics`] with a reserved capacity equal to the number of samples that it will track.
    pub fn with_capacity(
        sample_barcode_hop_checker: Option<&'a SampleBarcodeHopChecker>,
        cap: usize,
    ) -> Self {
        Self {
            per_sample_metrics: vec![DemuxedGroupSampleMetrics::default(); cap],
            sample_barcode_hop_tracker: sample_barcode_hop_checker
                .map(SampleBarcodeHopTracker::new),
            ..DemuxedGroupMetrics::default()
        }
    }

    /// Update this group with the results from another group.
    ///
    /// **Note** that no checking is done to ensure that both [`DemuxedGroupMetrics`] originate from the same set of
    /// samples. It is assumed that they line up correctly.
    pub fn update_with(&mut self, other: Self) {
        self.num_reads_failing_filters += other.num_reads_failing_filters;
        self.num_reads_filtered_as_control += other.num_reads_filtered_as_control;
        self.total_templates += other.total_templates;
        for (s, o) in self.per_sample_metrics.iter_mut().zip(other.per_sample_metrics.into_iter()) {
            s.update_with(&o);
        }

        if let Some(ref mut sample_barcode_hop_tracker) = self.sample_barcode_hop_tracker {
            if let Some(other_sample_barcode_hop_tracker) = other.sample_barcode_hop_tracker {
                sample_barcode_hop_tracker.update_with_self(other_sample_barcode_hop_tracker);
            }
        }
    }

    /// Write the metrics files associated with the [`DemuxedGroupMetrics`].
    ///
    /// This will create a `per_sample_metrics.tsv` file, a `per_project_metrics.tsv`,
    /// a `metrics.tsv` file, and optionally an `sample_barcode_hop_metrics.tsv`
    /// file (if dual-indexed) in the provided `output_dir`.
    pub fn write_metrics_files<P: AsRef<Path>>(
        self,
        samples: &[SampleMetadata],
        output_dir: P,
        prefix: &str,
    ) -> Result<()> {
        // Write the top level metrics file
        let run_metrics = RunMetrics {
            control_reads_omitted: self.num_reads_filtered_as_control,
            failing_reads_omitted: self.num_reads_failing_filters,
            total_templates: self.total_templates,
        };

        let filename = [prefix.to_string(), "metrics.tsv".to_string()].concat();
        let output_path = output_dir.as_ref().join(filename);
        let delim = DelimFile::default();
        delim.write_tsv(&output_path, std::iter::once(run_metrics))?;

        // Write the per sample metrics
        let filename = [prefix.to_string(), "per_sample_metrics.tsv".to_string()].concat();
        write_per_group_metrics(samples, &self.per_sample_metrics, &output_dir, &filename)?;

        // Build then write the per-project metrics
        let (per_project_samples, per_project_metrics) =
            build_per_project_metrics(samples, &self.per_sample_metrics)?;
        let filename = [prefix.to_string(), "per_project_metrics.tsv".to_string()].concat();
        write_per_group_metrics(
            &per_project_samples,
            &per_project_metrics,
            &output_dir,
            &filename,
        )?;

        // Optionally write the index hop file
        if let Some(sample_barcode_hop_tracker) = self.sample_barcode_hop_tracker {
            let filename =
                [prefix.to_string(), "sample_barcode_hop_metrics.tsv".to_string()].concat();
            let output_path = output_dir.as_ref().join(filename);

            let b1_len = sample_barcode_hop_tracker.checker.index1_length;
            delim.write_tsv(
                &output_path,
                sample_barcode_hop_tracker
                    .count_tracker
                    .into_iter()
                    .map(|(barcode, count)| {
                        let b1 = BString::from(&barcode[..b1_len]);
                        let b2 = BString::from(&barcode[b1_len..]);
                        let barcode = BString::from(format!("{}+{}", b1, b2));
                        BarcodeCount::new(barcode, count as isize)
                    })
                    .sorted_unstable_by_key(|b| -b.count),
            )?;
        }

        Ok(())
    }
}

/// The high level metrics for all reads that were demuxed.
#[derive(Debug, Serialize, Deserialize)]
pub struct RunMetrics {
    /// The number of reads that were omitted for being control reads.
    pub(crate) control_reads_omitted: usize,
    /// The number of reads that were omitted for having failed QC.
    pub(crate) failing_reads_omitted: usize,
    /// The total number of template reads that were output.
    pub(crate) total_templates: usize,
}

/// The per-sample metrics
#[derive(Debug, Default, Clone)]
pub struct DemuxedGroupSampleMetrics {
    /// Metrics about the base qualities of both the index and template bases for reads attributed to this sample.
    pub base_qual_counter: BaseQualCounter,
    /// The number of reads with perfect barcode matches.
    pub perfect_matches: usize,
    /// The number of reads with only a single mismatch.
    pub one_mismatch_matches: usize,
    /// The total number of reads attributed to this sample.
    pub total_matches: usize,
}

impl DemuxedGroupSampleMetrics {
    /// Update this [`DemuxedGroupSampleMetrics`] with another [`DemuxedGroupSampleMetrics`].
    pub fn update_with(&mut self, other: &Self) {
        self.base_qual_counter.update_with_self(&other.base_qual_counter);
        self.perfect_matches += other.perfect_matches;
        self.one_mismatch_matches += other.one_mismatch_matches;
        self.total_matches += other.total_matches;
    }

    /// Update the match counters with a [`MatchResult`].
    pub fn update_with_match(&mut self, match_result: &MatchResult) {
        if let MatchResult::Match { hamming_dist, .. } = *match_result {
            if hamming_dist == 0 {
                self.perfect_matches += 1;
            } else if hamming_dist == 1 {
                self.one_mismatch_matches += 1;
            }
        }
        self.total_matches += 1;
    }

    /// Convert [`DemuxedGroupSampleMetrics`] into [`SampleMetricsProcessed`].
    fn as_processed(
        &self,
        total_templates: usize,
        best_barcode_template_count: usize,
        sample_metadata: &SampleMetadata,
    ) -> SampleMetricsProcessed {
        SampleMetricsProcessed {
            barcode_name: sample_metadata.sample_id.clone(),
            library_name: sample_metadata.sample_id.clone(),
            barcode: sample_metadata.get_semantic_barcode().to_string(),
            templates: self.total_matches,
            perfect_matches: self.perfect_matches,
            one_mismatch_matches: self.one_mismatch_matches,
            q20_bases: self.base_qual_counter.q20_bases,
            q30_bases: self.base_qual_counter.q30_bases,
            total_number_of_bases: self.base_qual_counter.bases,
            fraction_matches: self.total_matches as f64 / total_templates as f64,
            ratio_this_barcode_to_best_barcode: self.total_matches as f64
                / best_barcode_template_count as f64,
            frac_q20_bases: self.base_qual_counter.q20_bases as f64
                / self.base_qual_counter.bases as f64,
            frac_q30_bases: self.base_qual_counter.q30_bases as f64
                / self.base_qual_counter.bases as f64,
            mean_index_base_quality: self.base_qual_counter.index_bases_qual_sum as f64
                / self.base_qual_counter.index_bases_total_seen as f64,
        }
    }
}

/// A helper object for tracking the quality scores for index and template bases.
#[derive(Debug, Default, Copy, Clone)]
pub struct BaseQualCounter {
    /// The count of bases with quality >= 20
    pub q20_bases: u64,
    /// The count of bases with quality >= 30
    pub q30_bases: u64,
    // The sum of the quality scores of the index bases for this sample
    pub index_bases_qual_sum: u64,
    // The total number of index bases seen per sample
    pub index_bases_total_seen: u64,
    // total template bases
    pub bases: u64,
}

impl BaseQualCounter {
    /// Create a new [`BaseQualCounter`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Update self based on another instance of [`BaseQualCounter`].
    pub fn update_with_self(&mut self, other: &Self) {
        self.q20_bases += other.q20_bases;
        self.q30_bases += other.q30_bases;
        self.index_bases_qual_sum += other.index_bases_qual_sum;
        self.index_bases_total_seen += other.index_bases_total_seen;
        self.bases += other.bases;
    }

    /// Update the counts based on a quality string.
    pub fn update(&mut self, quals: &[u8], kind: &SegmentType) {
        if *kind == SegmentType::SampleBarcode {
            self.index_bases_total_seen += quals.len() as u64;
            for q in quals {
                let q = q - 33;
                self.index_bases_qual_sum += q as u64;
            }
        } else if *kind == SegmentType::Template {
            self.bases += quals.len() as u64;
            for q in quals {
                let q = q - 33;
                if q >= 30 {
                    self.q30_bases += 1;
                    self.q20_bases += 1;
                } else if q >= 20 {
                    self.q20_bases += 1;
                }
            }
        }
    }

    /// Update the counter by attaching to an iterator, this elides the need to iterate over
    /// the quals twice in the event that they are being used.
    pub fn update_with_iter<'a, 's>(
        &'s mut self,
        quals: impl Iterator<Item = &'a u8> + 'a,
        kind: &SegmentType,
    ) -> Box<dyn Iterator<Item = &'a u8> + 'a>
    where
        's: 'a,
    {
        if *kind == SegmentType::SampleBarcode {
            Box::new(quals.map(|q| {
                let q_adjusted = q - 33;
                self.index_bases_qual_sum += q_adjusted as u64;
                self.index_bases_total_seen += 1;
                q
            }))
        } else if *kind == SegmentType::Template {
            Box::new(quals.map(|q| {
                let q_adjusted = q - 33;
                self.bases += 1;
                if q_adjusted >= 30 {
                    self.q30_bases += 1;
                    self.q20_bases += 1;
                } else if q_adjusted >= 20 {
                    self.q20_bases += 1;
                }
                q
            }))
        } else {
            Box::new(quals)
        }
    }
}

/// [`SampleBarcodeHopChecker`] checks to see if an unmatched barcode for a dual index read pair
/// is actually a concatenation of an exact match for an index1 from one sample and an index2
/// from a separate channel.
///
/// This is done by creating a lookup of index1 barcodes to sample and index2 barcodes to sample.
/// An unmatched read is then split and each half is checked against each lookup.
#[derive(Debug)]
pub struct SampleBarcodeHopChecker {
    index1_length: usize,
    index2_length: usize,
    // A lookup of index1 to the set of samples that use that index1.
    index1_map: AHashMap<Vec<u8>, AHashSet<usize>>,
    // A lookup of index2 to the set of samples that use that index2.
    index2_map: AHashMap<Vec<u8>, AHashSet<usize>>,
}

impl SampleBarcodeHopChecker {
    /// Attempt to create an [`SampleBarcodeHopChecker`].
    ///
    /// This will fail if there are any repeat index1 sequences or repeat index2 sequences across samples.
    pub fn try_new(
        samples: &[SampleMetadata],
        read_structures: &[ReadStructure],
    ) -> Result<Option<Self>> {
        // Collect all sample instances of a sample barcode
        let sample_barcodes: Vec<Vec<ReadSegment>> = read_structures
            .iter()
            .map(|r| r.sample_barcodes().copied().collect())
            .filter(|read: &Vec<ReadSegment>| !read.is_empty())
            .collect();

        // Check if this is a dual index setup
        if sample_barcodes.len() == 2 {
            let index1_length = sample_barcodes[0]
                .iter()
                .map(|s| s.length().expect("Sample barcode segment must have a length"))
                .sum();
            let index2_length = sample_barcodes[1]
                .iter()
                .map(|s| s.length().expect("Sample barcode segment must have a length"))
                .sum();

            // For each sample barcode in `SampleMetadata` split the barcode based on the above lengths and put each half in the hashmap
            let mut index1_map = AHashMap::with_capacity(samples.len());
            let mut index2_map = AHashMap::with_capacity(samples.len());
            for (i, sample) in samples.iter().enumerate() {
                let index1 = sample.barcode[0..index1_length].to_vec();
                let index2 = sample.barcode[index1_length..].to_vec();

                let index1_entry = index1_map.entry(index1).or_insert_with(AHashSet::new);
                index1_entry.insert(i);
                let index2_entry = index2_map.entry(index2).or_insert_with(AHashSet::new);
                index2_entry.insert(i);
            }

            Ok(Some(Self { index1_length, index2_length, index1_map, index2_map }))
        } else {
            Ok(None)
        }
    }
}

/// [`SampleBarcodeHopTracker`] keeps track of the count of each barcode that is identified as an indexhop.
#[derive(Debug)]
pub struct SampleBarcodeHopTracker<'a> {
    /// The [`SampleBarcodeHopChecker`] instance for this run.
    checker: &'a SampleBarcodeHopChecker,
    /// The hashmap that tracks the number of times each indexhop occurs
    count_tracker: AHashMap<Vec<u8>, usize>,
}

impl<'a> SampleBarcodeHopTracker<'a> {
    /// Create a new [`SampleBarcodeHopTracker`].
    pub fn new(checker: &'a SampleBarcodeHopChecker) -> Self {
        Self { checker, count_tracker: AHashMap::new() }
    }

    /// Update this instance of [`SampleBarcodeHopTracker`] with the contents of another, merging counters.
    pub fn update_with_self(&mut self, other: Self) {
        for (key, value) in other.count_tracker {
            let counter = self.count_tracker.entry(key).or_insert(0);
            *counter += value;
        }
    }

    /// Get the number keys (barcodes) in the [`SampleBarcodeHopTracker`].
    pub fn len(&self) -> usize {
        self.count_tracker.len()
    }

    /// Check if the [`SampleBarcodeHopTracker`] is empty.
    pub fn is_empty(&self) -> bool {
        self.count_tracker.is_empty()
    }

    /// Check for index hopping.
    ///
    /// Check for occurrences where index1 is an exact match to a known barcode and index2 is an exact match to a different
    /// known barcode.
    ///
    /// **Note**: This assumes that the barcode being passed in has already been checked for exact matches against the known barcodes.
    /// This method does not check if the barcode "halves" both match the same sample (which would not be a hop, that would just be
    /// an exact match).
    pub fn check_barcode(&mut self, barcode: &[u8]) {
        // Verify that the barcode length matches the expected length
        if barcode.len() == self.checker.index1_length + self.checker.index2_length {
            let index1a = &barcode[0..self.checker.index1_length];
            let index2a = &barcode[self.checker.index1_length..];

            let index1b = &barcode[self.checker.index2_length..];
            let index2b = &barcode[0..self.checker.index2_length];

            if (self.checker.index1_map.contains_key(index1a)
                && self.checker.index2_map.contains_key(index2a))
                || (self.checker.index1_map.contains_key(index1b)
                    && self.checker.index2_map.contains_key(index2b))
            {
                let counter = self.count_tracker.entry(barcode.to_vec()).or_insert(0);
                *counter += 1;
            }
        }
    }
}

/// A helper struct fot serializing and deserializing barcode counts.
#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
pub(crate) struct BarcodeCount {
    /// The barcode.
    pub(crate) barcode: String,
    /// the count of the barcode.
    pub(crate) count: isize,
}

impl BarcodeCount {
    /// Create a new [`BarcodeCount`] instance.
    fn new(barcode: BString, count: isize) -> Self {
        Self { barcode: barcode.to_string(), count }
    }
}

/// The final result of aggregating the [`DemuxedGroupSampleMetrics`] and computing metrics from raw counts.
#[derive(Debug, Serialize, Deserialize)]
pub struct SampleMetricsProcessed {
    /// The name for the sample barcode, typically the same name from the SampleSheet.
    pub(crate) barcode_name: String,
    /// The name of the library, typically the library identifier from the SampleSheet.
    pub(crate) library_name: String,
    /// The sample barcode bases. Dual index barcodes will have two sample barcode sequences delimited by a `+`.
    pub(crate) barcode: String,
    /// The total number of templates matching the given barcode.
    pub(crate) templates: usize,
    /// The number of templates that match perfectly the given barcode.
    pub(crate) perfect_matches: usize,
    /// The number of pass-filter templates that match the given barcode with exactly one mismatch.
    pub(crate) one_mismatch_matches: usize,
    /// The number of bases in a template with a quality score 20 or above.
    pub(crate) q20_bases: u64,
    /// The number of bases in a template with a quality score 30 or above.
    pub(crate) q30_bases: u64,
    /// The total number of bases in the templates combined.
    pub(crate) total_number_of_bases: u64,
    /// The fraction of all templates that match the given barcode.
    pub(crate) fraction_matches: f64,
    /// The rate of all templates matching this barcode to all template reads matching the most prevalent barcode.
    /// For the most prevalent barcode this will be 1, for all others it will be less than 1 (except for the possible
    /// exception of when there are more unmatched templates than for any other barcode, in which case the value may
    /// be arbitrarily large). One over the lowest number in this column gives you the fold-difference in representation
    /// between barcodes.
    pub(crate) ratio_this_barcode_to_best_barcode: f64,
    /// The fraction of bases in a template with a quality score 20 or above.
    pub(crate) frac_q20_bases: f64,
    /// The fraction of bases in a template with a quality score 30 or above.
    pub(crate) frac_q30_bases: f64,
    /// The mean quality of index bases.
    pub(crate) mean_index_base_quality: f64,
}
