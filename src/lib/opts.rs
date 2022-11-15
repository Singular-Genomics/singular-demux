#![forbid(unsafe_code)]

use std::{num::NonZeroUsize, path::PathBuf, str::FromStr, vec::Vec};

use anyhow::{anyhow, Result};
use clap::Parser;
use env_logger::Env;
use read_structure::{ReadSegment, ReadStructure, SegmentType};

use crate::{
    demux::DemuxReadFilterConfig,
    matcher::{MatcherKind, UNDETERMINED_NAME},
    utils::{built_info, infer_fastq_sequence_length},
};

pub static LOGO: &str = "
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

pub static TOOL_NAME: &str = "sgdemux";

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

Per-sample files with suffixes like _R1.fastq.gz will be written to the output directory specified with --output.

The sample metadata file may be a Sample Sheet or a simple two-column CSV file with headers.
The Sample Sheet may haave a `[Demux]` section for command line options, and must have a `[Data]`
section for sample information.  The `Sample_ID` column must contain a unique, non-empty identifier
for each sample.  Both `Index1_Sequence` and `Index2_Sequence` must be present with values for
indexed runs.  For non-indexed runs, a single sample must be given with an empty value for the 
`Index1_Sequence` and `Index2_Sequence` columns.  For the simple two-column CSV, the 
`Sample_Barcode` column must contain the unique set of sample barcode bases for the sample(s).

Example invocation:

sgdemux \\
  --fastqs R1.fq.gz R2.fq.gz I1.fq.gz \\
  --read-structures +T +T 8B \\
  --sample-metadata samples.csv \\
  --output demuxed-fastqs/

For complete documentation see: https://github.com/Singular-Genomics/singular-demux
For support please contact: care@singulargenomics.com
";

#[derive(Parser, Debug, Clone)]
#[clap(name = TOOL_NAME, version = built_info::VERSION.as_str(), about=SHORT_USAGE, long_about=LONG_USAGE, term_width=0)]
pub struct Opts {
    /// Path to the input FASTQs or path prefix.
    #[clap(long, short = 'f', display_order = 1, required = true, multiple_values = true)]
    pub fastqs: Vec<PathBuf>,

    /// Path to the sample metadata.
    #[structopt(long, short = 's', display_order = 2)]
    pub sample_metadata: PathBuf,

    /// Read structures, one per input FASTQ.
    #[clap(long, short = 'r', display_order = 3, required = true, multiple_values = true)]
    pub read_structures: Vec<ReadStructure>,

    /// The directory to write outputs.
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

    /// Mask template bases with quality scores less than specified value(s).
    ///
    /// Sample /// barcode/index and UMI bases are never masked. If provided either a single value,
    /// or one value per FASTQ must be provided.
    #[clap(long, short = 'M', required = false, multiple = true, display_order = 11)]
    pub quality_mask_threshold: Vec<u8>,

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
        // Make a vec of quality mask thresholds that is always the same length as the number
        // of files being demultiplexed.
        let num_reads = self.read_structures.len();
        let num_values = self.quality_mask_threshold.len();

        // This should be validated on construction so panic here if incorrect number of
        // masking thresholds.
        if num_values != 0 && num_values != 1 && num_values != num_reads {
            panic!("Opts.quality_mask_thresholds must have 0, 1 or num_fastqs values.");
        }

        let thresholds: Vec<u8> = if num_values == num_reads {
            self.quality_mask_threshold.clone()
        } else if num_values == 1 {
            vec![self.quality_mask_threshold[0]; num_reads]
        } else {
            vec![0; num_reads]
        };

        DemuxReadFilterConfig::new(
            self.filter_control_reads,
            self.filter_failing_quality,
            thresholds,
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

    /// Builds a new copy of this `Opt` where read structures with variable length sample barcodes
    /// are replaced with fixed length sample barcodes.  The fixed length of each sample barcode
    /// is inferred from the FASTQ corresponding to the given read structure by examining the
    /// length of the first read, then removing any other fixed length segments, with the remaining
    /// length used as the fixed sample barcode length.
    ///
    /// If there is no read structure with a variable length sample barcode, the `Opt` itself is
    /// returned.
    pub fn with_fixed_sample_barcodes(self) -> Result<Self> {
        if self
            .read_structures
            .iter()
            .all(|s| s.sample_barcodes().all(read_structure::ReadSegment::has_length))
        {
            Ok(self)
        } else {
            // Go through each read structure and try to infer the length of any variable length
            // sample barcode
            let mut read_structures: Vec<ReadStructure> = vec![];
            for (read_structure, fastq) in self.read_structures.iter().zip(self.fastqs.iter()) {
                if read_structure
                    .iter()
                    .all(|s| s.kind != SegmentType::SampleBarcode || s.has_length())
                {
                    // No variable length, so just clone the current read structure
                    read_structures.push(read_structure.clone());
                } else {
                    // Get the read length from the FASTQ, subtract all non variable length
                    // segments (must be a sample barcode if we've gone this far).
                    let read_length: usize = infer_fastq_sequence_length(fastq.clone())?;
                    let fixed_length: usize =
                        read_structure.iter().map(|s| s.length().unwrap_or(0)).sum();

                    let read_segments: Vec<ReadSegment> = read_structure
                        .iter()
                        .filter_map(|s| {
                            if s.has_length() {
                                Some(*s)
                            } else if fixed_length == read_length {
                                None
                            } else {
                                Some(
                                    ReadSegment::from_str(&format!(
                                        "{}{}",
                                        read_length - fixed_length,
                                        s.kind.value()
                                    ))
                                    .unwrap(),
                                )
                            }
                        })
                        .collect();
                    read_structures.push(ReadStructure::new(read_segments)?);
                }
            }
            Ok(Opts { read_structures, ..self })
        }
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
            quality_mask_threshold: vec![0],
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
        fs::File,
        io::{BufWriter, Write},
        path::PathBuf,
        str::FromStr,
    };

    use gzp::{BgzfSyncWriter, Compression};
    use read_structure::ReadStructure;
    use tempfile::tempdir;

    use super::Opts;

    struct FastqDef {
        read_structure: ReadStructure,
        fastq: PathBuf,
        read_length: usize,
        expected: ReadStructure,
    }

    #[test]
    fn test_with_fixed_sample_barcodes() {
        let dir = tempdir().unwrap();

        let fastq_defs = vec![
            // no change: no sample barcode, template is variable length
            FastqDef {
                read_structure: ReadStructure::from_str("+T").unwrap(),
                fastq: dir.path().join("1.fastq.gz"),
                read_length: 150,
                expected: ReadStructure::from_str("+T").unwrap(),
            },
            // no change: no sample barcode, fixed molecluar barcode, template is variable length
            FastqDef {
                read_structure: ReadStructure::from_str("24M+T").unwrap(),
                fastq: dir.path().join("3.fastq.gz"),
                read_length: 150,
                expected: ReadStructure::from_str("24M+T").unwrap(),
            },
            // no change: fixed sample barcode length
            FastqDef {
                read_structure: ReadStructure::from_str("+B").unwrap(),
                fastq: dir.path().join("2.fastq.gz"),
                read_length: 20,
                expected: ReadStructure::from_str("20B").unwrap(),
            },
            // updated: one fixed length and one variable length sample barcode
            FastqDef {
                read_structure: ReadStructure::from_str("8B10S+B").unwrap(),
                fastq: dir.path().join("4.fastq.gz"),
                read_length: 30,
                expected: ReadStructure::from_str("8B10S12B").unwrap(),
            },
            // updated: read length is short enough that variable length sample barcode is removed
            FastqDef {
                read_structure: ReadStructure::from_str("8B10S+B").unwrap(),
                fastq: dir.path().join("5.fastq.gz"),
                read_length: 18,
                expected: ReadStructure::from_str("8B10S").unwrap(),
            },
        ];

        // Create output FASTQ files
        for fastq_def in &fastq_defs {
            let writer = BufWriter::new(File::create(&fastq_def.fastq).unwrap());
            let mut gz_writer = BgzfSyncWriter::new(writer, Compression::new(3));
            let bases = String::from_utf8(vec![b'A'; fastq_def.read_length]).unwrap();
            let quals = String::from_utf8(vec![b'I'; fastq_def.read_length]).unwrap();
            let data = format!("@NAME\n{}\n+\n{}\n", bases, quals);
            gz_writer.write_all(data.as_bytes()).unwrap();
            gz_writer.flush().unwrap();
            drop(gz_writer);
        }

        // Create the opt, and update it
        let opt = Opts {
            read_structures: fastq_defs.iter().map(|f| f.read_structure.clone()).collect(),
            fastqs: fastq_defs.iter().map(|f| f.fastq.clone()).collect(),
            ..Opts::default()
        }
        .with_fixed_sample_barcodes()
        .unwrap();

        // Check the new read structures
        let expected_read_structures: Vec<ReadStructure> =
            fastq_defs.iter().map(|f| f.expected.clone()).collect();
        assert_eq!(opt.read_structures.len(), expected_read_structures.len());
        for (actual, expected) in opt.read_structures.iter().zip(expected_read_structures.iter()) {
            assert_eq!(actual.to_string(), expected.to_string());
        }
    }
}
