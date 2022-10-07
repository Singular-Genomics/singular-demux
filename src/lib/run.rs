use std::{fs::File, io::BufWriter, num::NonZeroUsize, path::PathBuf, sync::Arc, vec::Vec};

use anyhow::{anyhow, bail, ensure, Context, Result};
use clap::Parser;
use env_logger::Env;
use gzp::BUFSIZE;
use itertools::Itertools;
use log::{debug, info};
use parking_lot::Mutex;
use pooled_writer::{bgzf::BgzfCompressor, Pool};
use rayon::prelude::*;
use read_structure::{ReadStructure, SegmentType};

use crate::{
    demux::{Demultiplex, Demultiplexer, DemuxReadFilterConfig, DemuxedGroup, PerFastqRecordSet},
    matcher::{CachedHammingDistanceMatcher, MatcherKind, PreComputeMatcher, UNDETERMINED_NAME},
    metrics::{DemuxedGroupMetrics, UnmatchedCounterThread},
    pooled_sample_writer::PooledSampleWriter,
    sample_metadata::{self},
    thread_reader::ThreadReader,
    utils::{built_info, check_bgzf, filenames, MultiZip},
};

static LOGO: &str = "
███████╗██╗███╗   ██╗ ██████╗ ██╗   ██╗██╗      █████╗ ██████╗
██╔════╝██║████╗  ██║██╔════╝ ██║   ██║██║     ██╔══██╗██╔══██╗
███████╗██║██╔██╗ ██║██║  ███╗██║   ██║██║     ███████║██████╔╝
╚════██║██║██║╚██╗██║██║   ██║██║   ██║██║     ██╔══██║██╔══██╗
███████║██║██║ ╚████║╚██████╔╝╚██████╔╝███████╗██║  ██║██║  ██║
╚══════╝╚═╝╚═╝  ╚═══╝ ╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝

 ██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗██╗ ██████╗███████╗
██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██║██╔════╝██╔════╝
██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║██║██║     ███████╗
██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██║██║     ╚════██║
╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║██║╚██████╗███████║
 ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝
";

static SHORT_USAGE: &str = "Performs sample demultiplexing on block-compressed (BGZF) FASTQs.";

static LONG_USAGE: &str = "
███████╗██╗███╗   ██╗ ██████╗ ██╗   ██╗██╗      █████╗ ██████╗
██╔════╝██║████╗  ██║██╔════╝ ██║   ██║██║     ██╔══██╗██╔══██╗
███████╗██║██╔██╗ ██║██║  ███╗██║   ██║██║     ███████║██████╔╝
╚════██║██║██║╚██╗██║██║   ██║██║   ██║██║     ██╔══██║██╔══██╗
███████║██║██║ ╚████║╚██████╔╝╚██████╔╝███████╗██║  ██║██║  ██║
╚══════╝╚═╝╚═╝  ╚═══╝ ╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝

 ██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗██╗ ██████╗███████╗
██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██║██╔════╝██╔════╝
██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║██║██║     ███████╗
██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██║██║     ╚════██║
╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║██║╚██████╗███████║
 ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝

Performs sample demultiplexing on block-compressed (BGZF) FASTQs.

Input FASTQs must be block compressed (e.g. with `bgzip`).  A single bgzipped FASTQ file
should be provided per instrument read.  One read structure should be provided per input FASTQ.

The output directory specified with --output must exist.  Per-sample files with suffixes like
_R1.fastq.gz will be written to the output directory.

The sample metadata file must be a two-column CSV file with headers.  The `Sample_ID` column
must contain a unique, non-empty identifier for each sample.  The `Sample_Barcode` column must
contain the unique set of sample barcode bases for the sample(s).

Example invocation:

sgdemux \\
  --fastqs R1.fq.gz R2.fq.gz I1.fq.gz \\
  --read-structures +T +T 8B \\
  --sample-metadata samples.csv \\
  --output demuxed-fastqs/

For complete documentation see: https://github.com/Singular-Genomics/singular-demux
For support please contact: care@singulargenomics.com
";

#[derive(Parser, Debug)]
#[clap(name = "sgdemux", version = built_info::VERSION.as_str(), about=SHORT_USAGE, long_about=LONG_USAGE, term_width=0)]
pub struct Opts {
    /// Path to the input FASTQs.
    #[clap(long, short = 'f', display_order = 1, required = true, multiple_values = true)]
    pub fastqs: Vec<PathBuf>,

    /// Path to the sample metadata.
    #[structopt(long, short = 's', display_order = 2)]
    pub sample_metadata: PathBuf,

    /// Read structures, one per input FASTQ.
    #[clap(long, short = 'r', display_order = 3, required = true, multiple_values = true)]
    pub read_structures: Vec<ReadStructure>,

    /// The directory to write outputs, the directory must exist.
    ///
    /// This tool will overwrite existing files.
    #[clap(long, short, display_order = 4)]
    pub output_dir: PathBuf,

    /// Number of allowed mismatches between the observed barcode and the expected barcode.
    #[clap(long, short = 'm', default_value = "1", display_order = 11)]
    pub allowed_mismatches: usize,

    /// The minimum allowed difference between an observed barcode and the second closest expected
    /// barcode.
    #[clap(long, short = 'd', default_value = "2", display_order = 11)]
    pub min_delta: usize,

    /// Number of N's to allow in a barcode without counting against the allowed_mismatches
    #[clap(long, short = 'F', default_value = "1", display_order = 11)]
    pub free_ns: usize,

    /// Max no-calls (N's) in a barcode before it is considered unmatchable.
    ///
    /// A barcode with total N's greater than `max_no_call` will be considered unmatchable.
    ///
    /// [default: None]
    #[clap(long, short = 'N', display_order = 11)]
    pub max_no_calls: Option<usize>,

    /// Mask bases with quality score less than or equal to this value.
    ///
    /// If set to zero, no masking will be performed.
    #[clap(long, short = 'M', default_value = "0", display_order = 11)]
    pub quality_mask_threshold: u8,

    /// Filter out control reads.
    #[clap(long, short = 'C', display_order = 11)]
    pub filter_control_reads: bool,

    /// Filter reads failing quality filter.
    #[clap(long, short = 'Q', display_order = 11)]
    pub filter_failing_quality: bool,

    /// The types of output FASTQs to write.
    ///
    /// These may be any of the following:
    /// - `T` - Template bases
    /// - `B` - Sample barcode bases
    /// - `M` - Molecular barcode bases
    /// - `S` - Skip bases
    ///
    /// For each read structure, all segment types listed by `--output-types` will be output to a
    /// FASTQ file.
    // TODO: add FromStr implementation for SegmentType so this can be a Vec<SegmentType>.  See: https://github.com/fulcrumgenomics/read-structure/issues/3
    #[clap(long, short = 'T', default_value = "T", verbatim_doc_comment, display_order = 21)]
    pub output_types: String,

    /// The sample name for undetermined reads (reads that do not match an expected barcode).
    #[clap(long, short = 'u', default_value = UNDETERMINED_NAME, display_order = 21)]
    pub undetermined_sample_name: String,

    /// Output the most frequent "unmatched" barcodes up to this number.
    ///
    /// If set to 0 unmatched barcodes will not be collected, improving overall performance.
    #[clap(long, short = 'U', default_value = "1000", display_order = 31)]
    pub most_unmatched_to_output: usize,

    /// Size of the channel for the unmatched barcode counter.
    ///
    /// Each item in the channel is a `Vec` of the barcodes of the unmatched reads. If all reads were unmatched
    /// the max number of barcodes in the Vec would be equal to `chunksize`.
    #[clap(long, default_value = "1000", display_order = 31, hide = true)]
    pub most_unmatched_channel_size: usize,

    /// Max number of keys the most unmatched hash map is allowed to contain.
    #[clap(long, default_value = "5000000", display_order = 31, hide = true)]
    pub most_unmatched_max_map_size: usize,

    /// Number of keys to shrink the most unmatched hash map down to when it hits the `most_unmatched_max_map_size`.
    #[clap(long, default_value = "5000", display_order = 31, hide = true)]
    pub most_unmatched_downsize_to: usize,

    /// The number of reads to extract from a FASTQ at one time.
    ///
    /// A "chunk" is the unit of parallelization for all of demultiplexing and defines how many reads are operated on at one time.
    #[clap(long, short = 'c', default_value = "1000", display_order = 31, hide = true)]
    pub chunksize: NonZeroUsize,

    /// Number of threads for demultiplexing.
    ///
    /// The number of threads to use for the process of determining which input reads should be assigned to which sample.
    #[clap(long, short = 't', default_value = "4", display_order = 31)]
    pub demux_threads: usize,

    /// Number of threads for compression the output reads.
    ///
    /// The number of threads to use for compressing reads that are queued for writing.
    #[clap(long, default_value = "12", display_order = 31)]
    pub compressor_threads: usize,

    /// Number of threads for writing compressed reads to output.
    ///
    /// The number of threads to have writing reads to their individual output files.
    #[clap(long, default_value = "5", display_order = 31)]
    pub writer_threads: usize,

    /// The number of threads to use for decompression for each reader.
    #[clap(long, default_value = "4", display_order = 31, hide = true)]
    pub decompression_threads_per_reader: usize,

    /// Override the matcher heuristic.
    ///
    /// If the sample barcodes are > 12 bp long, a cached hamming distance matcher is used.
    /// If the barcodes are less than or equal to 12 bp long, all possible matches are precomputed.
    ///
    /// This option allows for overriding that heuristic.
    ///
    /// [default: None]
    #[clap(long, possible_values=MatcherKind::possible_values(), display_order = 31)]
    pub override_matcher: Option<MatcherKind>,
}

impl Opts {
    /// Extract a [`DemuxReadFilterConfig`] from the CLI opts.
    pub fn as_read_filter_config(&self) -> DemuxReadFilterConfig {
        DemuxReadFilterConfig::new(
            self.filter_control_reads,
            self.filter_failing_quality,
            self.quality_mask_threshold,
            self.max_no_calls,
        )
    }

    /// Get the [`SegmentType`]s to write to file from the CLI opts.
    pub fn output_types_to_write(&self) -> Result<Vec<SegmentType>> {
        let mut output_types_to_write = vec![];
        for b in self.output_types.as_bytes() {
            output_types_to_write.push(SegmentType::try_from(*b).map_err(|msg| anyhow!(msg))?);
        }
        Ok(output_types_to_write)
    }
}

/// Implement defaults that match the CLI options to allow for easier testing.
///
/// Note that these defaults exist only within test code.
#[cfg(test)]
impl Default for Opts {
    fn default() -> Self {
        Self {
            fastqs: vec![],
            read_structures: vec![],
            allowed_mismatches: 2,
            min_delta: 1,
            free_ns: 1,
            filter_control_reads: false,
            filter_failing_quality: false,
            most_unmatched_to_output: 1_000,
            most_unmatched_channel_size: 1_000,
            most_unmatched_max_map_size: 5_000_000,
            most_unmatched_downsize_to: 5_000,
            max_no_calls: Some(1),
            quality_mask_threshold: 0,
            output_types: String::from("T"),
            undetermined_sample_name: UNDETERMINED_NAME.to_string(),
            sample_metadata: PathBuf::default(),
            chunksize: NonZeroUsize::new(500).unwrap(),
            demux_threads: 16,
            compressor_threads: 4,
            writer_threads: 4,
            decompression_threads_per_reader: 4,
            override_matcher: None,
            output_dir: PathBuf::default(),
        }
    }
}

/// Run demultiplexing.
#[allow(clippy::too_many_lines)]
pub fn run(opts: Opts) -> Result<(), anyhow::Error> {
    eprint!("{}", LOGO);

    let read_filter_config = opts.as_read_filter_config();

    let output_types_to_write = opts.output_types_to_write()?;
    let samples = sample_metadata::from_path(
        opts.sample_metadata,
        Some(opts.allowed_mismatches),
        true,
        &opts.undetermined_sample_name,
    )?;

    // Preflight checks
    ensure!(
        opts.output_dir.exists(),
        "Output directory does not exist: {}",
        &opts.output_dir.to_string_lossy()
    );
    ensure!(!opts.read_structures.is_empty(), "At least one read structure must be specified");
    ensure!(
        opts.read_structures.iter().any(|s| s.sample_barcodes().count() > 0),
        "No sample barcodes found in read structures"
    );
    ensure!(
        opts.read_structures.iter().any(|s| s.templates().count() > 0),
        "No templates found in read structures"
    );
    ensure!(
        opts.fastqs.len() == opts.read_structures.len(),
        "Same number of read structures should be given as FASTQs"
    );
    ensure!(
        samples.iter().all(|s| s.barcode.len() == samples[0].barcode.len()),
        "Sample metadata barcodes are unequal lengths."
    );
    ensure!(
        opts.read_structures
            .iter()
            .map(|s| s.sample_barcodes().map(|b| b.length().unwrap_or(0)).sum::<usize>())
            .sum::<usize>()
            == samples[0].barcode.len(),
        "The number of sample barcode bases in read structures does not match sample metadata"
    );

    // Check that the input FASTQs are in BGZF format
    for fastq in &opts.fastqs {
        check_bgzf(fastq)?;
    }

    info!("Creating writer threads");
    let writers: Result<Vec<_>> = samples
        .iter()
        .sorted_by(|a, b| a.ordinal.cmp(&b.ordinal))
        .map(|s| {
            let writers: Result<Vec<_>> =
                filenames(s, &opts.output_dir, &opts.read_structures, &output_types_to_write)
                    .into_iter()
                    .map(|name| {
                        File::create(&name)
                            .with_context(|| {
                                format!("Unable to create file: {}", name.to_string_lossy())
                            })
                            .map(|f| BufWriter::with_capacity(BUFSIZE, f))
                    })
                    .collect();
            writers
        })
        .collect();
    let writers = writers?;
    let group_size = writers[0].len();
    let (mut pool, pooled_writers) = Pool::new::<_, BgzfCompressor>(
        opts.writer_threads,
        opts.compressor_threads,
        2,
        writers.into_iter().flatten().collect(),
    )?;
    let mut writers = vec![];
    for grouped_writers in &pooled_writers.into_iter().chunks(group_size) {
        writers.push(Arc::new(Mutex::new(PooledSampleWriter::new(
            grouped_writers.collect::<Vec<_>>(),
        )?)));
    }

    info!("Creating reader threads");
    let mut readers = opts
        .fastqs
        .iter()
        .map(|f| {
            ThreadReader::new(f.clone(), opts.chunksize, opts.decompression_threads_per_reader)
        })
        .collect::<Vec<_>>();
    let rpool = rayon::ThreadPoolBuilder::new().num_threads(opts.demux_threads).build().unwrap();

    let unmatched_counter = if opts.most_unmatched_to_output > 0 {
        Some(UnmatchedCounterThread::new(
            Some(opts.most_unmatched_channel_size),
            Some(opts.most_unmatched_max_map_size),
            Some(opts.most_unmatched_downsize_to),
        )?)
    } else {
        None
    };

    info!("Processing data");
    let demuxer: Box<dyn Demultiplex> = if samples[0].barcode.len() <= 12
        || opts.override_matcher.map_or(false, |m| m == MatcherKind::PreCompute)
    {
        debug!("Using PreComputeMatcher");
        let matcher = PreComputeMatcher::new(
            &samples[0..samples.len() - 1],
            opts.allowed_mismatches,
            opts.min_delta,
            opts.free_ns,
        );
        let demuxer = Demultiplexer::new(
            &samples,
            &opts.read_structures,
            &output_types_to_write,
            &read_filter_config,
            matcher,
            true,
        )?;
        Box::new(demuxer)
    } else {
        debug!("Using CachedHammingDistanceMatcher");
        let matcher = CachedHammingDistanceMatcher::new(
            &samples[0..samples.len() - 1],
            opts.allowed_mismatches,
            opts.min_delta,
            opts.free_ns,
        );
        let demuxer = Demultiplexer::new(
            &samples,
            &opts.read_structures,
            &output_types_to_write,
            &read_filter_config,
            matcher,
            true,
        )?;
        Box::new(demuxer)
    };

    let metrics: Result<DemuxedGroupMetrics> = rpool.install(|| {
        // This unwrap can't fail since the tx will exist until we call `collect`.
        let unmatched_counter_tx = unmatched_counter
            .as_ref()
            .map(|unmatched_counter| unmatched_counter.tx.as_ref().unwrap().clone());

        let iterators = readers.iter_mut().map(|r| r.rx.iter()).collect();
        MultiZip::new(iterators)
            .par_bridge()
            .map(|reads_super_chunks| {
                let grouped_set = PerFastqRecordSet::new(reads_super_chunks)?;
                demuxer.demultiplex(&grouped_set)
            })
            .map(|records: Result<DemuxedGroup>| {
                match records {
                    Ok(records) => {
                        for (sample_index, sample_reads) in
                            records.per_sample_reads.into_iter().enumerate()
                        {
                            let writer = &mut writers[sample_index].lock();
                            writer.write_records(sample_reads).with_context(|| {
                                format!("Failed to write reads to sample number: {}", sample_index)
                            })?;
                        }
                        if let Some(unmatched) = records.unmatched {
                            // The records.unmatched will only be Some if we are collecting `unmatched_counter`, so the unwrap here is safe
                            unmatched_counter_tx
                                .as_ref()
                                .unwrap()
                                .send(unmatched)
                                .context("Failed to send to unmatched counter thread.")?;
                        }
                        Ok(records.group_metrics)
                    }
                    Err(err) => Err(err),
                }
            })
            .fold(
                || {
                    Ok(DemuxedGroupMetrics::with_capacity(
                        demuxer.demux_hop_checker(),
                        samples.len(),
                    ))
                },
                |all: Result<DemuxedGroupMetrics>, d: Result<DemuxedGroupMetrics>| match (all, d) {
                    (Ok(mut all), Ok(d)) => {
                        all.update_with(d);
                        Ok(all)
                    }
                    (Ok(_), Err(d)) => Err(d),
                    (Err(a), Err(_) | Ok(_)) => Err(a),
                },
            )
            .reduce(
                || {
                    Ok(DemuxedGroupMetrics::with_capacity(
                        demuxer.demux_hop_checker(),
                        samples.len(),
                    ))
                },
                |all: Result<DemuxedGroupMetrics>, d: Result<DemuxedGroupMetrics>| match (all, d) {
                    (Ok(mut all), Ok(d)) => {
                        all.update_with(d);
                        Ok(all)
                    }
                    (Ok(_), Err(d)) => Err(d),
                    (Err(a), Err(_) | Ok(_)) => Err(a),
                },
            )
    });
    let metrics = metrics?;

    info!("Joining reader threads");
    for reader in readers {
        match reader.handle.join() {
            Ok(result) => result?,
            Err(e) => std::panic::resume_unwind(e),
        }
    }

    info!("Joining writer threads and collecting stats");
    for w in writers {
        let w = match Arc::try_unwrap(w) {
            Ok(w) => w,
            Err(_e) => bail!("Writer lock is till locked"), // This should never happen
        };
        let w = w.into_inner();
        w.finish().context("Failed to flush and finish writing.")?;
    }

    info!("Stopping pool");
    pool.stop_pool()?;

    info!("Writing stats");
    metrics.write_metrics_files(&samples, &opts.output_dir)?;
    if let Some(unmatched_counter) = unmatched_counter {
        unmatched_counter.collect().to_file(&opts.output_dir, opts.most_unmatched_to_output)?;
    }
    Ok(())
}

/// Parse args and set up logging / tracing
pub fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}

#[cfg(test)]
mod test {
    use std::{
        fs::create_dir,
        path::{Path, PathBuf},
        str::FromStr,
    };

    use fgoxide::io::{DelimFile, Io};
    use read_structure::ReadStructure;
    use rstest::rstest;
    use seq_io::{fastq::OwnedRecord, BaseRecord};

    use crate::{
        fastq_header::FastqHeader,
        matcher::{MatcherKind, UNDETERMINED_NAME},
        metrics::{BarcodeCount, RunMetrics, SampleMetricsProcessed},
        sample_metadata::{self, SampleMetadata},
        utils::{
            segment_kind_to_fastq_id,
            test_commons::{
                create_preset_sample_metadata_file, slurp_fastq, write_reads_to_file, Fq,
                SAMPLE_BARCODE_1, SAMPLE_BARCODE_4,
            },
        },
    };

    use super::{run, Opts};

    const SINGLE_FASTQ_QUAL_RANGE: &[u8] = &[
        35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
        58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    ];

    fn fastq_path_single(dir: impl AsRef<Path>) -> PathBuf {
        let path = dir.as_ref().join("test.fastq.gz");
        let record = Fq {
            name: "frag1",
            bases: &[SAMPLE_BARCODE_1, &[b'A'; 39]].concat(),
            quals: Some(&[&[b'?'; 17], SINGLE_FASTQ_QUAL_RANGE].concat()),
            ..Fq::default()
        }
        .to_owned_record();
        write_reads_to_file(std::iter::once(record), &path);
        path
    }

    fn test_end_to_end_with_quality_threshold(
        min_base_qual_for_masking: u8,
        threads: usize,
        matcher: Option<MatcherKind>,
    ) -> Vec<OwnedRecord> {
        let dir = tempfile::tempdir().unwrap();
        let read_structure = ReadStructure::from_str("17B39T").unwrap();
        let input = fastq_path_single(&dir.path());
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![input],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures: vec![read_structure],
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            quality_mask_threshold: min_base_qual_for_masking,
            override_matcher: matcher,
            ..Opts::default()
        };

        run(opts).unwrap();

        let fastq = output.join(format!("{}_R1.fastq.gz", "Sample1"));
        slurp_fastq(&fastq)
    }

    #[rstest]
    #[case(0, 1, &[b'A'; 39])]
    #[case(10, 1, &[[b'N'; 9].as_slice(), &[b'A'; 30]].concat())]
    #[case(0, 2, &[b'A'; 39])]
    #[case(10, 2, &[[b'N'; 9].as_slice(), &[b'A'; 30]].concat())]
    #[case(39, 2, &[[b'N'; 38].as_slice(), &[b'A'; 1]].concat())]
    #[case(100, 2, &[b'N'; 39])]
    fn test_run_end_to_end_with_quality_threshold(
        #[case] quality: u8,
        #[case] threads: usize,
        #[case] expected: &[u8],
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        // https://github.com/fulcrumgenomics/fgbio/blob/ed906cd29a7f38c9da0583c8307dc7fbcd94880b/src/test/scala/com/fulcrumgenomics/fastq/DemuxFastqsTest.scala#L470
        // Note that this tool masks anything less than or equal to the quality threshold, but fgbio masks only less than
        // so this will mask one more base than fgbio
        let records = test_end_to_end_with_quality_threshold(quality, threads, matcher);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, expected);
    }

    #[rstest]
    #[should_panic(expected = "Same number of read structures should be given as FASTQs")]
    #[case(vec![ReadStructure::from_str("8B100T").unwrap(), ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(expected = "No sample barcodes found in read structures")]
    #[case(vec![ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(expected = "No templates found in read structures")]
    #[case(vec![ReadStructure::from_str("8B").unwrap()])]
    #[should_panic(expected = "Same number of read structures should be given as FASTQs")]
    #[case(vec![ReadStructure::from_str("8B92T").unwrap(), ReadStructure::from_str("100T").unwrap(), ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(
        expected = "The number of sample barcode bases in read structures does not match sample metadata"
    )]
    #[case(vec![ReadStructure::from_str("18B100T").unwrap()])]
    fn test_read_structure_failures(
        #[case] read_structures: Vec<ReadStructure>,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let input = fastq_path_single(&dir.path());
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![input],
            output_dir: output,
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            override_matcher: matcher,
            ..Opts::default()
        };
        run(opts).unwrap();
    }

    #[rstest]
    #[should_panic(
        expected = " Invalid barcode sequence for Sample5 `` - Barcode is an empty string. Line 6"
    )]
    fn test_missing_barcode_in_sample_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let input = fastq_path_single(&dir.path());
        let output = dir.path().join("output");
        let read_structures = vec![ReadStructure::from_str("17B100T").unwrap()];
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());
        let io = Io::default();
        let mut metadata_lines = io.read_lines(&metadata).unwrap();
        metadata_lines.push(String::from("Sample5,"));
        let buggy_metadata = dir.path().join("buggy_metadata.csv");
        io.write_lines(&buggy_metadata, metadata_lines).unwrap();

        let opts = Opts {
            fastqs: vec![input],
            output_dir: output,
            sample_metadata: buggy_metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            ..Opts::default()
        };
        run(opts).unwrap();
    }

    fn fastq_path(dir: impl AsRef<Path>) -> PathBuf {
        let path = dir.as_ref().join("test.fastq.gz");
        let records = vec![
            Fq {
                name: "frag1", // matches the first sample -> first sample
                comment: Some("1:Y:0:SampleNumber"),
                bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag2", // matches the first sample, one mismatch -> first sample
                comment: Some("2:Y:0:SampleNumber"),
                bases: &[b"AAAAAAAAGATTACAGT".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag3", // matches the first sample, three mismatches -> unmatched
                comment: Some("3:Y:0:SampleNumber"),
                bases: &[b"AAAAAAAAGATTACTTT".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag4", // matches the fourth sample perfectly and 3rd barcode with two mismathces, delta too small -> unmatched
                comment: Some("4:Y:0:SampleNumber"),
                bases: &[SAMPLE_BARCODE_4, &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag5", // matches the first sample, too many Ns -> unmatched
                comment: Some("5:Y:0:SampleNumber"),
                bases: &[b"AAAAAAAAGANNNNNNN".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
        ];

        write_reads_to_file(records.into_iter(), &path);
        path
    }

    fn names_and_barcodes(records: &[OwnedRecord]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        records
            .iter()
            .map(|r| {
                let header = FastqHeader::try_from(r.head()).unwrap();
                let name =
                    header.comment.as_ref().unwrap().other.as_ref().unwrap().clone().into_owned();
                let barcode = header.comment.unwrap().info.sample_barcode.into_owned();
                (name, barcode)
            })
            .unzip()
    }

    #[rstest]
    fn test_end_to_end_simple(
        #[values(1, 2)] threads: usize,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
        #[values("T", "B", "TB")] output_types: String,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let read_structure = ReadStructure::from_str("17B100T").unwrap();
        let input = fastq_path(&dir.path());
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![input],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures: vec![read_structure],
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            output_types: output_types,
            ..Opts::default()
        };

        let samples = sample_metadata::from_path(
            &opts.sample_metadata,
            Some(opts.allowed_mismatches),
            true,
            &opts.undetermined_sample_name,
        )
        .unwrap();

        let output_types_to_write = opts.output_types_to_write().unwrap();

        run(opts).unwrap();

        // Check outputs
        for name in samples.iter().map(|s| &s.sample_id) {
            for output_type_to_write in &output_types_to_write {
                let fastq_id = segment_kind_to_fastq_id(output_type_to_write);
                let fastq = output.join(format!("{}_{}1.fastq.gz", name, fastq_id));
                let records = slurp_fastq(&fastq);
                let (names, barcodes) = names_and_barcodes(&records);

                if *name == "Sample1" {
                    assert_eq!(records.len(), 2);
                    assert!(names.contains(&b"frag1".to_vec()));
                    assert!(names.contains(&b"frag2".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGATTACAGA".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGATTACAGT".to_vec()));
                } else if *name == UNDETERMINED_NAME {
                    assert_eq!(records.len(), 3);
                    assert!(names.contains(&b"frag3".to_vec()));
                    assert!(names.contains(&b"frag4".to_vec()));
                    assert!(names.contains(&b"frag5".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGATTACTTT".to_vec()));
                    assert!(barcodes.contains(&b"GGGGGGTTGATTACAGA".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGANNNNNNN".to_vec()));
                } else {
                    assert!(records.is_empty());
                }
            }
        }

        // Check metrics
        let per_sample_metrics = output.join("per_sample_metrics.tsv");
        let run_metrics = output.join("metrics.tsv");
        let most_unmatched = output.join("most_frequent_unmatched.tsv");

        let delim = DelimFile::default();
        let per_sample_metrics: Vec<SampleMetricsProcessed> =
            delim.read_tsv(&per_sample_metrics).unwrap();
        assert_eq!(per_sample_metrics.len(), 5);

        for (metric, sample) in per_sample_metrics.into_iter().zip(samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == samples[0].sample_id {
                assert_eq!(metric.templates, 2);
                assert_eq!(metric.perfect_matches, 1);
                assert_eq!(metric.one_mismatch_matches, 1);
            } else if metric.barcode_name == samples[4].sample_id {
                assert_eq!(metric.templates, 3);
                assert_eq!(metric.perfect_matches, 0);
                assert_eq!(metric.one_mismatch_matches, 0);
            } else {
                assert_eq!(metric.templates, 0);
                assert_eq!(metric.perfect_matches, 0);
                assert_eq!(metric.one_mismatch_matches, 0);
            }
        }

        let run_metrics: Vec<RunMetrics> = delim.read_tsv(&run_metrics).unwrap();
        assert_eq!(run_metrics.len(), 1);
        assert_eq!(run_metrics[0].control_reads_omitted, 0);
        assert_eq!(run_metrics[0].failing_reads_omitted, 0);
        assert_eq!(run_metrics[0].total_templates, 5);

        let most_frequent_unmatched: Vec<BarcodeCount> = delim.read_tsv(&most_unmatched).unwrap();
        assert_eq!(most_frequent_unmatched.len(), 3);
        assert!(most_frequent_unmatched.iter().all(|b| b.count == 1));
    }

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_end_to_end(
        #[values(1, 2)] threads: usize,
        #[values(true, false)] omit_failing_reads: bool,
        #[values(true, false)] omit_control_reads: bool,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
            ReadStructure::from_str("117T").unwrap(),
        ];
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let filter_flag = if omit_failing_reads { "Y" } else { "N" };
        let control_flag = if omit_control_reads { "1" } else { "0" };

        let path = dir.as_ref().join("test.fastq.gz");
        let records = vec![
            Fq {
                name: "1", // matches the first sample -> first sample
                comment: Some(&format!("1:{}:{}:SampleNumber", filter_flag, control_flag)),
                bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "2", // matches the first sample, one mismatch -> first sample
                comment: Some("2:N:0:SampleNumber"),
                bases: &[b"AAAAAAAAGATTACAGT".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "3", // matches the first sample, three mismatches -> unmatched
                comment: Some(&format!("3:N:{}:SampleNumber", control_flag)),
                bases: &[b"AAAAAAAAGATTACTTT".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "4", // matches the fourth sample perfectly and 3rd barcode with two mismatches, delta too small -> unmatched
                comment: Some("4:N:0:SampleNumber"),
                bases: &[SAMPLE_BARCODE_4, &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "5", // matches the first sample, too many Ns -> unmatched
                comment: Some("5:N:0:SampleNumber"),
                bases: &[b"AAAAAAAAGANNNNNNN".as_slice(), &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record(),
        ];

        let barcodes_per_sample = vec![
            if omit_failing_reads || omit_control_reads {
                vec![b"AAAAAAAAGATTACAGT".as_slice()]
            } else {
                vec![SAMPLE_BARCODE_1, b"AAAAAAAAGATTACAGT".as_slice()]
            }, // sample1
            vec![], // sample2
            vec![], // sample3
            vec![], // sample4
            if omit_control_reads {
                vec![SAMPLE_BARCODE_4, b"AAAAAAAAGANNNNNNN".as_slice()]
            } else {
                vec![
                    b"AAAAAAAAGATTACTTT".as_slice(),
                    SAMPLE_BARCODE_4,
                    b"AAAAAAAAGANNNNNNN".as_slice(),
                ]
            }, // undetermined
        ];

        let assignments_per_sample = vec![
            if omit_failing_reads || omit_control_reads { vec!["2"] } else { vec!["1", "2"] }, // sample1
            vec![],                                                                // sample2
            vec![],                                                                // sample3
            vec![],                                                                // sample4
            if omit_control_reads { vec!["4", "5"] } else { vec!["3", "4", "5"] }, // undetermined
        ];

        write_reads_to_file(records.into_iter(), &path);
        let input = path;
        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![input.clone(), input],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            filter_control_reads: omit_control_reads,
            filter_failing_quality: omit_failing_reads,
            ..Opts::default()
        };

        let samples = sample_metadata::from_path(
            &opts.sample_metadata,
            Some(opts.allowed_mismatches),
            true,
            &opts.undetermined_sample_name,
        )
        .unwrap();

        run(opts).unwrap();

        for (i, sample) in samples.iter().enumerate() {
            let barcodes: Vec<_> = barcodes_per_sample[i].iter().cloned().collect();
            let assignment: Vec<_> =
                assignments_per_sample[i].iter().map(|s| (*s).as_bytes()).collect();
            let sample_id = sample.sample_id.clone();
            let r1_fastq = output.join(format!("{}_R1.fastq.gz", &sample_id));
            let r2_fastq = output.join(format!("{}_R2.fastq.gz", &sample_id));
            let r1_records = slurp_fastq(&r1_fastq);
            let r2_records = slurp_fastq(&r2_fastq);
            let (r1_names, r1_barcodes) = names_and_barcodes(&r1_records);
            let (r2_names, r2_barcodes) = names_and_barcodes(&r2_records);

            // Ensure R1 and R2 line up
            assert_eq!(r1_names, r2_names);
            assert_eq!(r1_names, assignment);
            assert_eq!(r2_names, assignment);
            assert_eq!(r1_barcodes, r2_barcodes);
            assert_eq!(r1_barcodes, barcodes);
            assert_eq!(r2_barcodes, barcodes);
        }
    }

    #[rstest]
    fn test_demux_fragment_reads_with_standard_qualities(
        #[values(1, 2)] threads: usize,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let read_structures = vec![ReadStructure::from_str("4B+T").unwrap()];
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let fq = Fq {
            name: "frag",
            bases: b"TCGTGACAGGAATCAAATGAAAACACTTGGT",
            quals: Some(b"1>1>11AFF?FFGGGGGGCBDGBGCFGH11B"),
            ..Fq::default()
        }
        .to_owned_record();
        let input = dir.path().join("test.fastq.gz");
        write_reads_to_file(std::iter::once(fq.clone()), &input);
        let delim = DelimFile::default();
        let metadata = dir.path().join("sample_metadata.csv");
        delim
            .write_csv(
                &metadata,
                SampleMetadata::new(String::from("foo"), b"TCGT".as_slice().into(), 0),
            )
            .unwrap();

        let opts = Opts {
            fastqs: vec![input],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            ..Opts::default()
        };

        run(opts).unwrap();

        let record = &slurp_fastq(output.join("foo_R1.fastq.gz"))[0];
        let header = FastqHeader::try_from(record.head()).unwrap();
        assert_eq!(header.comment.unwrap().other.unwrap().into_owned(), b"frag");
        assert_eq!(record.seq(), fq.seq().iter().copied().skip(4).collect::<Vec<u8>>());
        assert_eq!(record.qual, fq.qual.iter().copied().skip(4).collect::<Vec<u8>>());
    }

    #[rstest]
    #[should_panic(
        expected = "Unequal number of reads in each record set (likely uneven input FASTQs)"
    )]
    fn test_demux_should_fail_if_one_fastq_has_fewer_records_than_the_other(
        #[values(1, 2)] threads: usize,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        fn to_fq(i: usize) -> OwnedRecord {
            Fq {
                name: &format!("frag{}", i),
                bases: &[SAMPLE_BARCODE_1, &[b'A'; 100]].concat(),
                ..Fq::default()
            }
            .to_owned_record()
        }

        let dir = tempfile::tempdir().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("17B100T").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
        ];
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let r1_input = dir.path().join("test_r1.fastq.gz");
        let r2_input = dir.path().join("test_r2.fastq.gz");
        write_reads_to_file([to_fq(1), to_fq(2)].into_iter(), &r1_input);
        write_reads_to_file([to_fq(1)].into_iter(), &r2_input);

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![r1_input, r2_input],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            ..Opts::default()
        };

        run(opts).unwrap();
    }

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_demux_dual_index_paired_end_reads(
        #[values(1, 2)] threads: usize,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let fq1_path = dir.path().join("fq1.fastq.gz");
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: b"AAAAAAAA", ..Fq::default() }.to_owned_record(),
            ),
            &fq1_path,
        );
        let fq2_path = dir.path().join("fq2.fastq.gz");
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: &[b'A'; 100], ..Fq::default() }.to_owned_record(),
            ),
            &fq2_path,
        );
        let fq3_path = dir.path().join("fq3.fastq.gz");
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: &[b'T'; 100], ..Fq::default() }.to_owned_record(),
            ),
            &fq3_path,
        );
        let fq4_path = dir.path().join("fq4.fastq.gz");
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: b"GATTACAGA", ..Fq::default() }.to_owned_record(),
            ),
            &fq4_path,
        );

        let read_structures = vec![
            ReadStructure::from_str("8B").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
            ReadStructure::from_str("9B").unwrap(),
        ];

        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let opts = Opts {
            fastqs: vec![fq1_path, fq2_path, fq3_path, fq4_path],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            filter_control_reads: false,
            ..Opts::default()
        };

        let samples = sample_metadata::from_path(
            &opts.sample_metadata,
            Some(opts.allowed_mismatches),
            true,
            &opts.undetermined_sample_name,
        )
        .unwrap();
        run(opts).unwrap();

        for (i, sample) in samples.iter().enumerate() {
            let sample_id = sample.sample_id.clone();
            let r1_fastq = output.join(format!("{}_R1.fastq.gz", &sample_id));
            let r2_fastq = output.join(format!("{}_R2.fastq.gz", &sample_id));
            let r1_records = slurp_fastq(&r1_fastq);
            let r2_records = slurp_fastq(&r2_fastq);

            let (r1_names, r1_barcodes) = names_and_barcodes(&r1_records);
            let (r2_names, r2_barcodes) = names_and_barcodes(&r2_records);

            // Ensure R1 and R2 line up
            assert_eq!(r1_names, r2_names);
            assert_eq!(r1_barcodes, r2_barcodes);

            if i == 0 {
                assert_eq!(r1_records.len(), 1);
                assert_eq!(r2_records.len(), 1);
                assert!(r1_barcodes.iter().all(|b| b
                    .iter()
                    .filter(|&&c| c != b'+')
                    .copied()
                    .collect::<Vec<_>>()
                    == SAMPLE_BARCODE_1));
                assert!(r2_barcodes.iter().all(|b| b
                    .iter()
                    .filter(|&&c| c != b'+')
                    .copied()
                    .collect::<Vec<_>>()
                    == SAMPLE_BARCODE_1));
            } else {
                assert!(r1_records.is_empty());
                assert!(r2_records.is_empty());
            }
        }

        let per_sample_metrics = output.join("per_sample_metrics.tsv");
        let delim = DelimFile::default();
        let per_sample_metrics: Vec<SampleMetricsProcessed> =
            delim.read_tsv(&per_sample_metrics).unwrap();
        assert_eq!(per_sample_metrics.len(), 5);
        for (metric, sample) in per_sample_metrics.into_iter().zip(samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == samples[0].sample_id {
                assert_eq!(metric.templates, 1);
                assert_eq!(metric.perfect_matches, 1);
                assert_eq!(metric.one_mismatch_matches, 0);
            } else if metric.barcode_name == samples[4].sample_id {
                assert_eq!(metric.templates, 0);
                assert_eq!(metric.perfect_matches, 0);
                assert_eq!(metric.one_mismatch_matches, 0);
            } else {
                assert_eq!(metric.templates, 0);
                assert_eq!(metric.perfect_matches, 0);
                assert_eq!(metric.one_mismatch_matches, 0);
            }
        }
    }

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_end_to_end_metrics(
        #[values(1, 2)] threads: usize,
        #[values(true, false)] omit_failing_reads: bool,
        #[values(true, false)] omit_control_reads: bool,
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().join("output");
        create_dir(&output).unwrap();
        let fq_path = dir.path().join("fq.fastq.gz");
        let fastqs = vec![
            Fq {
                name: "frag1",
                bases: &[SAMPLE_BARCODE_1, &[b'A'; 39]].concat(),
                quals: Some(&[&[b'?'; 17], SINGLE_FASTQ_QUAL_RANGE].concat()),
                comment: Some("1:Y:0:SampleNumber"),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag2",
                bases: &[SAMPLE_BARCODE_1, &[b'G'; 39]].concat(),
                quals: Some(&[&[b'?'; 17], SINGLE_FASTQ_QUAL_RANGE].concat()),
                comment: Some("2:N:0:SampleNumber"),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag3",
                bases: &[SAMPLE_BARCODE_1, &[b'G'; 39]].concat(),
                quals: Some(&[&[b'?'; 17], SINGLE_FASTQ_QUAL_RANGE].concat()),
                comment: Some("3:Y:1:SampleNumber"),
                ..Fq::default()
            }
            .to_owned_record(),
            Fq {
                name: "frag4",
                bases: &[SAMPLE_BARCODE_1, &[b'G'; 39]].concat(),
                quals: Some(&[&[b'?'; 17], SINGLE_FASTQ_QUAL_RANGE].concat()),
                comment: Some("4:N:1:SampleNumber"),
                ..Fq::default()
            }
            .to_owned_record(),
        ];
        write_reads_to_file(fastqs.into_iter(), &fq_path);
        let metadata = create_preset_sample_metadata_file(&dir.path());
        let opts = Opts {
            fastqs: vec![fq_path],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures: vec![ReadStructure::from_str("17B39T").unwrap()],
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            override_matcher: matcher,
            filter_control_reads: omit_control_reads,
            filter_failing_quality: omit_failing_reads,
            ..Opts::default()
        };

        let samples = sample_metadata::from_path(
            &opts.sample_metadata,
            Some(opts.allowed_mismatches),
            true,
            &opts.undetermined_sample_name,
        )
        .unwrap();

        run(opts).unwrap();

        // Setup expected
        // Note, if a read is both control and failed quality, it is counted as failed quality
        let (
            templates,
            failing_reads_ommited,
            control_reads_omitted,
            total_num_bases,
            q30_above,
            q20_above,
        ) = if omit_failing_reads && omit_control_reads {
            (1, 2, 1, 39, 11, 21)
        } else if omit_failing_reads && !omit_control_reads {
            (2, 2, 0, 78, 22, 42)
        } else if !omit_failing_reads && omit_control_reads {
            (2, 0, 2, 78, 22, 42)
        } else if !omit_failing_reads && !omit_control_reads {
            (4, 0, 0, 156, 44, 84)
        } else {
            unreachable!()
        };

        // Check metrics
        let per_sample_metrics = output.join("per_sample_metrics.tsv");
        let run_metrics = output.join("metrics.tsv");
        let most_unmatched = output.join("most_frequent_unmatched.tsv");

        let delim = DelimFile::default();
        let per_sample_metrics: Vec<SampleMetricsProcessed> =
            delim.read_tsv(&per_sample_metrics).unwrap();
        assert_eq!(per_sample_metrics.len(), 5);

        for (metric, sample) in per_sample_metrics.into_iter().zip(samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == samples[0].sample_id {
                assert_eq!(metric.templates, templates);
                assert_eq!(metric.perfect_matches, templates);
                assert_eq!(metric.one_mismatch_matches, 0);
                assert_eq!(metric.total_number_of_bases, total_num_bases);
                assert_eq!(metric.q30_bases, q30_above);
                assert_eq!(metric.q20_bases, q20_above);
                assert_eq!(metric.frac_q20_bases, q20_above as f64 / total_num_bases as f64);
                assert_eq!(metric.frac_q30_bases, q30_above as f64 / total_num_bases as f64);
            } else {
                assert_eq!(metric.templates, 0);
                assert_eq!(metric.perfect_matches, 0);
                assert_eq!(metric.one_mismatch_matches, 0);
                assert_eq!(metric.total_number_of_bases, 0);
                assert_eq!(metric.q30_bases, 0);
                assert_eq!(metric.q20_bases, 0);
                assert!(metric.frac_q20_bases.is_nan());
                assert!(metric.frac_q30_bases.is_nan());
            }
        }

        let run_metrics: Vec<RunMetrics> = delim.read_tsv(&run_metrics).unwrap();
        eprintln!("omit-failing {:?} omit-control {:?}", omit_failing_reads, omit_control_reads);
        eprintln!("run_metrics {:?} ", run_metrics);
        assert_eq!(run_metrics.len(), 1);
        assert_eq!(run_metrics[0].control_reads_omitted, control_reads_omitted);
        assert_eq!(run_metrics[0].failing_reads_omitted, failing_reads_ommited);
        assert_eq!(run_metrics[0].total_templates, templates);

        let most_frequent_unmatched: Vec<BarcodeCount> = delim.read_tsv(&most_unmatched).unwrap();
        assert_eq!(most_frequent_unmatched.len(), 0);
    }
}
