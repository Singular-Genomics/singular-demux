#![forbid(unsafe_code)]

use std::{
    cmp::Ordering, collections::HashMap, fs::File, io::BufReader, num::NonZeroUsize, path::PathBuf,
    str::FromStr, vec::Vec,
};

use anyhow::{anyhow, bail, ensure, Context, Result};
use bstr::ByteSlice;
use clap::Parser;
use env_logger::Env;
use gzp::{BgzfSyncReader, BUFSIZE};
use itertools::Itertools;
use read_structure::{ReadSegment, ReadStructure, SegmentType, ANY_LENGTH_STR};
use seq_io::{fastq, BaseRecord};

use crate::{
    demux::DemuxReadFilterConfig,
    matcher::{MatcherKind, UNDETERMINED_NAME},
    utils::{built_info, infer_fastq_sequence_length, InputFastq},
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

pub static DOC_LINK_AND_SUPPORT_EMAIL: &str = "
For complete documentation see: https://github.com/Singular-Genomics/singular-demux
For support please contact: care@singulargenomics.com
";

#[derive(Parser, Debug, Clone)]
#[clap(name = TOOL_NAME, version = built_info::VERSION.as_str(), about=SHORT_USAGE, long_about=LONG_USAGE, term_width=0)]
pub struct Opts {
    /// Path to the input FASTQs, or path prefix if not a file.
    #[clap(long, short = 'f', display_order = 1, required = true, multiple_values = true)]
    pub fastqs: Vec<PathBuf>,

    /// Path to the sample metadata.
    #[structopt(long, short = 's', display_order = 2)]
    pub sample_metadata: PathBuf,

    /// Read structures, one per input FASTQ. Do not provide when using a path prefix for FASTQs.
    #[clap(long, short = 'r', display_order = 3, multiple_values = true)]
    pub read_structures: Vec<ReadStructure>,

    /// The directory to write outputs to.
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
    /// Sample barcode/index and UMI bases are never masked. If provided either a single value,
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

    /// If this is true, then all the read names across FASTQs will not be enforced to be the same.
    /// This may be useful when the read names are known to be the same and performance matters.
    /// Regardless, the first read name in each FASTQ will always be checked.
    #[clap(long, display_order = 31)]
    pub skip_read_name_check: bool,

    /// If this is true, then the sample barcode is expected to be in the FASTQ read header.  For
    /// dual indexed data, the barcodes must be `+` (plus) delimited.  Additionally, if true, then
    /// neither index FASTQ files nor sample barcode segments in the read structure may be
    /// specified.
    #[clap(long, display_order = 32)]
    pub sample_barcode_in_fastq_header: bool,
}

impl Opts {
    /// Extracts one or more [`FastqsAndReadStructure`]s given the opts.
    ///
    /// If the input FASTQs are all files then there must be one read structure per FASTQ.  The
    /// groups returned will have a single FASTQ and read structure per.
    ///
    /// If the input FASTQs are path prefix(es) then read structure must be empty (not given).  The
    /// FASTQ files are gathered using the path prefix(es), with their kind (e.g. sample barcode or
    /// template), kind number (e.g. read one or index two), and lane inferred from the file name.
    /// Reads across differnt lanes with the same kind and kind number are grouped, and a
    /// corresponding read structure is built assuming all the bases in that read are of the given
    /// kind.  The groups are then sorted by kind and kind number to proper ordering of sample
    /// barcodes, UMIs, and template reads (i.e. read pairs).
    ///
    /// If the input FASTQS are a mix of files and path prefixes, this function will return an
    /// error.
    pub fn from(
        fastqs: Vec<PathBuf>,
        read_structures: Vec<ReadStructure>,
        sample_barcode_in_fastq_header: bool,
    ) -> Result<Vec<FastqsAndReadStructure>> {
        // must have FASTQs
        ensure!(!fastqs.is_empty(), "No FASTQs or path prefixes found with --fastq");

        // If the input FASTQs are in fact path prefixes, then glob in the FASTQs and create the
        // read structure based on the kind/kind-number inferred from the file name.
        let input_fastq_group: Vec<FastqsAndReadStructure> = if fastqs.iter().all(|f| f.is_file()) {
            // must have read structures, and match the # of FASTQs
            ensure!(
                fastqs.len() == read_structures.len(),
                "Same number of read structures should be given as FASTQs"
            );

            // Ensure that if the sample barcodes are to be extracted from the FASTQ header then
            // no sample barcode segments should be given
            if sample_barcode_in_fastq_header {
                ensure!(read_structures.iter().all(|g| g.sample_barcodes().count() == 0),
                    "Read structures may not contain sample barcode segment(s) when extracting sample barcodes from the FASTQ header");
            }

            FastqsAndReadStructure::zip(&fastqs, &read_structures)
        } else if fastqs.iter().all(|f| !f.is_file()) {
            ensure!(
                read_structures.is_empty(),
                "Read Structure must not be given when the input FASTQs are a path prefix."
            );

            let input_fastq_group = FastqsAndReadStructure::from_prefixes(&fastqs);

            // Ensure that if the sample barcodes are to be extracted from the FASTQ header then
            // no index reads (sample barcode reads) should be found.
            if sample_barcode_in_fastq_header {
                for fastqs_and_read_structure in &input_fastq_group {
                    if fastqs_and_read_structure.read_structure.sample_barcodes().count() > 0 {
                        bail!("Index reads found when extracting sample barcodes from the FASTQ header: {}",
                        fastqs_and_read_structure
                            .fastqs
                            .iter()
                            .map(|fastq| format!("{:?}", fastq.to_string_lossy())).join(", ")
                            );
                    }
                }
            }

            // must find some FASTQs since we ensured it above
            ensure!(
                !input_fastq_group.is_empty(),
                "Bug: No FASTQS found for prefix(es): {:?}.",
                fastqs
            );

            input_fastq_group
        } else {
            bail!("Input FASTQS (--fastqs) must either all be files or all path prefixes, not a mixture of both")
        };

        // Ensure that every group has the same number of FASTQs!
        for group in &input_fastq_group {
            ensure!(
                group.fastqs.len() == input_fastq_group[0].fastqs.len(),
                "Different # of FASTQs per group\nGroup 1: {:?}\nGroup 2: {:?}",
                input_fastq_group[0].fastqs,
                group.fastqs
            );
        }

        // Ensure that all FASTQs have the same read name in the header.  This checks only the
        // first read in each FASTQ
        let mut first_head: Vec<u8> = vec![];
        let mut end_index: usize = 0;
        let mut first_file: String = String::from("");
        for group in &input_fastq_group {
            for file in &group.fastqs {
                // Open the reader
                let reader = BufReader::with_capacity(
                    BUFSIZE,
                    File::open(&file)
                        .with_context(|| format!("Failed to open {}", file.to_string_lossy()))?,
                );
                let mut reader = fastq::Reader::with_capacity(BgzfSyncReader::new(reader), BUFSIZE);

                // Check the first record
                match reader.next() {
                    Some(Ok(record)) => {
                        if first_head.is_empty() {
                            first_head = record.head().to_vec();
                            end_index = first_head.find_byte(b' ').unwrap_or(first_head.len());
                            first_file = file.to_string_lossy().to_string();
                        } else {
                            let head = record.head();
                            let ok = head.len() == end_index
                                || (head.len() > end_index && head[end_index] == b' ');
                            let ok = ok && first_head[0..end_index] == head[0..end_index];
                            ensure!(
                                ok,
                                "Mismatching read names in the FASTQS:\n{} in {}\n{} in {}",
                                std::str::from_utf8(&first_head)?,
                                first_file,
                                std::str::from_utf8(record.head())?,
                                file.to_string_lossy()
                            );
                        }
                    }
                    _ => bail!("Empty input FASTQ: {}", file.to_string_lossy()),
                }
            }
        }

        // If there is a read structure that's all sample barcode, we need to replace it with the
        // expected length to enable index hopping metrics.  Do so by inspecting the first read in the
        // corresponding FASTQ
        input_fastq_group.iter().map(|g| g.clone().with_fixed_sample_barcodes()).collect()
    }

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
                    let read_length: usize = infer_fastq_sequence_length(&fastq)?;
                    let fixed_length: usize =
                        read_structure.iter().map(|s| s.length().unwrap_or(0)).sum();
                    match fixed_length.cmp(&read_length) {
                        Ordering::Greater => bail!(
                            "Read length ({}) is too short ({}) for the read structure {} in {:?}",
                            read_length,
                            fixed_length,
                            read_structure,
                            fastq
                        ),
                        Ordering::Equal => bail!("Variable length sample barcode would be zero (read length: {}, sum of fixed segment lengths: {}) in FASTQ: {:?}",
                            read_length, fixed_length, fastq
                        ),
                        Ordering::Less => {
                            let read_segments: Vec<ReadSegment> = read_structure
                                .iter()
                                .map(|s| {
                                    if s.has_length() {
                                        *s
                                    } else {
                                        ReadSegment::from_str(&format!(
                                            "{}{}",
                                            read_length - fixed_length,
                                            s.kind.value()
                                        ))
                                        .unwrap()
                                    }
                                })
                                .collect();
                            read_structures.push(ReadStructure::new(read_segments)?);
                        }
                    }
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
            skip_read_name_check: false,
            sample_barcode_in_fastq_header: false,
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

/// A set of FASTQs with the same kind (e.g. sample barcode, template) and kind number (e.g. read
/// one or index two).  This FASTQs differ based on path prefix and lane number.
#[derive(Debug, Clone)]
pub struct FastqsAndReadStructure {
    pub fastqs: Vec<PathBuf>,
    pub read_structure: ReadStructure,
}

impl FastqsAndReadStructure {
    /// Zips the path to FASTQs and associated read structures into a `FastqsAndReadStructure`.  
    /// The number of FASTQs and read structures must be the same, but is not validated here.
    pub fn zip(
        fastqs: &[PathBuf],
        read_structures: &[ReadStructure],
    ) -> Vec<FastqsAndReadStructure> {
        fastqs
            .iter()
            .zip(read_structures.iter())
            .map(|(f, r)| FastqsAndReadStructure {
                fastqs: vec![f.clone()],
                read_structure: r.clone(),
            })
            .collect()
    }

    /// Finds all FASTQs with a given file pattern (see `INPUT_FASTQ_REGEX` in `utils.rs`) and
    /// groups them by kind (`ReadSegment`) and kind number (e.g. read 1 or 2).
    ///
    /// After the FASTQs are discovered, the FASTQs are grouped and sorted by kind and kind number.
    /// Then each group is sorted by prefix and lane.  These steps ensure that (1) the group of
    /// FASTQS with the same kind and kind number are read in serially, (2) the FASTQS are
    /// synchronized across groups (e.g. the same lane/prefix across kind/kind-number will be read
    /// at the same time), and (3) we emit the read pairs in the correct order (e.g. R1/R2).
    pub fn from_prefixes(prefixes: &[PathBuf]) -> Vec<FastqsAndReadStructure> {
        // glob all the FASTQs
        let input_fastqs: Vec<InputFastq> =
            prefixes.iter().flat_map(|prefix| InputFastq::glob(&prefix).unwrap()).collect();

        // Collect all FASTQs by kind and kind number
        let mut unsorted_groups: HashMap<(SegmentType, u32), Vec<InputFastq>> = HashMap::new();
        for f in input_fastqs {
            let key = (f.kind, f.kind_number);
            unsorted_groups.entry(key).or_insert_with(|| Vec::with_capacity(1)).push(f);
        }

        // Sort across each group of FASTQS by kind and kind number.  For each group, sort by
        // prefix and lane.
        unsorted_groups
            .into_iter()
            .sorted_by(
                |((left_kind, left_kind_number), _), ((right_kind, right_kind_number), _)| {
                    left_kind.cmp(right_kind).then(left_kind_number.cmp(right_kind_number))
                },
            )
            .map(|((kind, _), fastqs)| {
                // this should always succeed
                let read_structure =
                    ReadStructure::from_str(&format!("{}{}", ANY_LENGTH_STR, kind.value()))
                        .unwrap();

                // sort FASTQs by prefix and lane for consistent ordering across groups
                let fastqs = fastqs
                    .iter()
                    .sorted_by(|left, right| {
                        left.prefix.cmp(&right.prefix).then(left.lane.cmp(&right.lane))
                    })
                    .map(|f| f.path.clone())
                    .collect();

                FastqsAndReadStructure { fastqs, read_structure }
            })
            .collect()
    }

    /// Builds a new copy of this `FastqsAndReadStructure` where a read structure with variable
    /// length sample barcodes is replaced with fixed length sample barcode.  The fixed length of
    /// the sample barcode is inferred from the first FASTQ by examining the length of the first
    /// read, then removing any other fixed length segments, with the remaining length used as
    /// the fixed sample barcode length.
    ///
    /// If there is no read structure with a variable length sample barcode, itself is
    /// returned.
    pub fn with_fixed_sample_barcodes(self) -> Result<Self> {
        if self.read_structure.sample_barcodes().all(read_structure::ReadSegment::has_length) {
            Ok(self)
        } else {
            // Get the read length from the FASTQ, subtract all non variable length
            // segments (must be a sample barcode if we've gone this far).
            let read_length: usize = infer_fastq_sequence_length(&self.fastqs[0])?;
            let fixed_length: usize =
                self.read_structure.iter().map(|s| s.length().unwrap_or(0)).sum();
            let read_structure = match fixed_length.cmp(&read_length) {
                Ordering::Greater => bail!(
                    "Read length ({}) is too short ({}) for the read structure {} in {:?}",
                    read_length,
                    fixed_length,
                    self.read_structure,
                    self.fastqs[0]
                ),
                Ordering::Equal => bail!("Variable length sample barcode would be zero (read length: {}, sum of fixed segment lengths: {}) in FASTQ: {:?}",
                    read_length, fixed_length,  self.fastqs[0]
                ),
                Ordering::Less => {
                    let read_segments: Vec<ReadSegment> = self.read_structure
                        .iter()
                        .map(|s| {
                            if s.has_length() {
                                *s
                            } else {
                                ReadSegment::from_str(&format!(
                                    "{}{}",
                                    read_length - fixed_length,
                                    s.kind.value()
                                ))
                                .unwrap()
                            }
                        })
                        .collect();
                    ReadStructure::new(read_segments)?
                }
            };
            Ok(FastqsAndReadStructure { fastqs: self.fastqs, read_structure })
        }
    }
}

#[cfg(test)]
mod test {
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufWriter, Write},
        path::{Path, PathBuf},
        str::FromStr,
    };

    use gzp::{BgzfSyncWriter, Compression};
    use read_structure::{ReadStructure, SegmentType};
    use rstest::rstest;
    use tempfile::tempdir;

    use crate::utils::segment_type_to_fastq_kind;

    use super::Opts;

    fn write(fastq: &Path, read_length: usize) {
        let writer = BufWriter::new(File::create(fastq).unwrap());
        let mut gz_writer = BgzfSyncWriter::new(writer, Compression::new(3));
        let bases = String::from_utf8(vec![b'A'; read_length]).unwrap();
        let quals = String::from_utf8(vec![b'I'; read_length]).unwrap();
        let data = format!("@NAME\n{}\n+\n{}\n", bases, quals);
        gz_writer.write_all(data.as_bytes()).unwrap();
        gz_writer.flush().unwrap();
    }

    fn write_with_read_name(fastq: &Path, read_name: &str) {
        let writer = BufWriter::new(File::create(fastq).unwrap());
        let mut gz_writer = BgzfSyncWriter::new(writer, Compression::new(3));
        let bases = String::from_utf8(vec![b'A'; 10]).unwrap();
        let quals = String::from_utf8(vec![b'I'; 10]).unwrap();
        let data = format!("@{}\n{}\n+\n{}\n", read_name, bases, quals);
        gz_writer.write_all(data.as_bytes()).unwrap();
        gz_writer.flush().unwrap();
    }

    struct FastqDef {
        read_structure: ReadStructure,
        fastq: PathBuf,
        read_length: usize,
        expected: ReadStructure,
    }

    impl FastqDef {
        fn write(&self) {
            write(&self.fastq, self.read_length);
        }
    }

    #[test]
    fn test_with_fixed_sample_barcodes_ok() {
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
        ];

        // Create output FASTQ files
        for fastq_def in &fastq_defs {
            fastq_def.write();
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

    #[test]
    fn test_with_fixed_sample_barcodes_error() {
        let dir = tempdir().unwrap();

        let fastq_defs = vec![
            // error: read length is short enough that variable length sample barcode is zero
            FastqDef {
                read_structure: ReadStructure::from_str("8B10S+B").unwrap(),
                fastq: dir.path().join("5.fastq.gz"),
                read_length: 18,
                expected: ReadStructure::from_str("+B").unwrap(), // ignore
            },
            // error: read length is short enough that leading fixed lengths are too short
            FastqDef {
                read_structure: ReadStructure::from_str("8B10S+B").unwrap(),
                fastq: dir.path().join("5.fastq.gz"),
                read_length: 17,
                expected: ReadStructure::from_str("+B").unwrap(), // ignore
            },
        ];

        // Create output FASTQ files
        for fastq_def in &fastq_defs {
            fastq_def.write();

            // Create the opt, and update it
            let result = Opts {
                read_structures: vec![fastq_def.read_structure.clone()],
                fastqs: vec![fastq_def.fastq.clone()],
                ..Opts::default()
            }
            .with_fixed_sample_barcodes();
            assert!(result.is_err());
        }
    }

    #[test]
    fn test_opts_from_error_mix_of_files_and_paths() {
        let dir = tempdir().unwrap();

        let file = dir.path().join("foo.fastq.gz");
        let prefix = dir.path().join("bar");

        // actually touch the path
        std::fs::File::create(file.clone()).unwrap();

        let opts = Opts {
            read_structures: vec![ReadStructure::from_str("+B").unwrap()],
            fastqs: vec![file, prefix],
            ..Opts::default()
        };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error.to_string().contains("must either all be files or all path prefixes"));
        }
    }

    #[test]
    fn test_opts_from_error_no_fastqs() {
        let opts = Opts {
            read_structures: vec![ReadStructure::from_str("+B").unwrap()],
            fastqs: vec![],
            ..Opts::default()
        };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error.to_string().contains("No FASTQs or path prefixes found"));
        }
    }

    #[test]
    fn test_opts_from_error_same_number_of_fastqs_as_read_structures() {
        let dir = tempdir().unwrap();

        let file1 = dir.path().join("foo.1.fastq.gz");
        let file2 = dir.path().join("foo.2.fastq.gz");

        // actually touch the path
        std::fs::File::create(file1.clone()).unwrap();
        std::fs::File::create(file2.clone()).unwrap();

        let opts = Opts {
            read_structures: vec![ReadStructure::from_str("+B").unwrap()],
            fastqs: vec![file1, file2],
            ..Opts::default()
        };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error
                .to_string()
                .contains("Same number of read structures should be given as FASTQs"));
        }
    }

    #[test]
    fn test_opts_from_error_read_structures_with_prefix() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("prefix");
        let opts = Opts {
            read_structures: vec![ReadStructure::from_str("+B").unwrap()],
            fastqs: vec![prefix],
            ..Opts::default()
        };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error.to_string().contains("Read Structure must not be given"));
        }
    }

    #[test]
    fn test_opts_from_error_no_fastqs_found_with_prefix() {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("prefix");
        let opts = Opts { read_structures: vec![], fastqs: vec![prefix], ..Opts::default() };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error.to_string().contains("No FASTQS found for prefix"));
        }
    }

    #[rstest]
    #[case("pre", & ["L001_R1_001.fastq.gz", "L002_R1_001.fastq.gz", "L001_R2_001.fastq.gz"])] // # lanes differ: R1 has L1, an R2 had both L1 & L2
    #[case("pre", & ["sub1_L001_R1_001.fastq.gz","sub1_L002_R1_001.fastq.gz","sub2_L001_R2_001.fastq.gz"])] // # of lanes differ: R1 has L1 & L2 (sub1), R2 has L1 (sub2)
    #[case("pre", & ["sub1_L001_R1_001.fastq.gz","sub2_L001_R1_001.fastq.gz","sub2_L001_I2_001.fastq.gz"])] // # of prefixes differ: R1 has L1 (sub1 & sub2), I1 has L1 (sub2)
    fn test_oppts_from_error_differing_number_of_fastqs_per_group_with_prefix(
        #[case] prefix: String,
        #[case] suffixes: &[&str],
    ) {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join(prefix);

        let mut fastqs: Vec<PathBuf> = vec![];
        for suffix in suffixes.iter() {
            let fastq = PathBuf::from(prefix.to_string_lossy().to_string() + "_" + suffix);
            std::fs::File::create(fastq.clone()).unwrap();
            fastqs.push(fastq);
        }

        let opts =
            Opts { read_structures: vec![], fastqs: vec![prefix.clone()], ..Opts::default() };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);
        assert!(result.is_err());
        if let Err(error) = result {
            assert!(error.to_string().contains("Different # of FASTQs per group"));
        }
    }

    fn compare_fastq_file_names(actual: &[PathBuf], expected: &[String]) {
        assert_eq!(actual.len(), expected.len());
        for (act, exp) in actual.iter().zip(expected.iter()) {
            assert_eq!(act.file_name().unwrap().to_string_lossy().to_string(), exp.to_string());
        }
    }

    #[test]
    fn test_optsfrom_ok_complex() {
        let dir = tempdir().unwrap();

        let mut prefixes: Vec<PathBuf> = vec![];
        let mut fastqs_map: HashMap<(SegmentType, u32), Vec<String>> = HashMap::new();
        for prefix_i in 1..=4 {
            let prefix = format!("prefix{}", prefix_i);
            prefixes.push(dir.path().join(prefix));
            for lane_i in 1..=prefix_i {
                for kind in [SegmentType::SampleBarcode, SegmentType::Template] {
                    for kind_number in 1..=2 {
                        let file_name = format!(
                            "prefix{}_L00{}_{}{}_001.fastq.gz",
                            prefix_i,
                            lane_i,
                            segment_type_to_fastq_kind(&kind),
                            kind_number
                        );
                        fastqs_map
                            .entry((kind, kind_number))
                            .or_insert_with(|| Vec::with_capacity(1))
                            .push(file_name.clone());
                        let fastq = dir.path().join(file_name);
                        write(&fastq, 50);
                    }
                }
            }
        }

        let opts = Opts { read_structures: vec![], fastqs: prefixes, ..Opts::default() };

        let results = Opts::from(opts.fastqs, opts.read_structures, false).unwrap();
        assert_eq!(results.len(), 4);

        // I1
        let result = &results[0];
        assert_eq!(result.read_structure.to_string(), "50B");
        compare_fastq_file_names(
            &result.fastqs,
            fastqs_map.get(&(SegmentType::SampleBarcode, 1)).unwrap(),
        );
        // I2
        let result = &results[1];
        assert_eq!(result.read_structure.to_string(), "50B");
        compare_fastq_file_names(
            &result.fastqs,
            fastqs_map.get(&(SegmentType::SampleBarcode, 2)).unwrap(),
        );

        // R1
        let result = &results[2];
        assert_eq!(result.read_structure.to_string(), "+T");
        compare_fastq_file_names(
            &result.fastqs,
            fastqs_map.get(&(SegmentType::Template, 1)).unwrap(),
        );
        // R2
        let result = &results[3];
        assert_eq!(result.read_structure.to_string(), "+T");
        compare_fastq_file_names(
            &result.fastqs,
            fastqs_map.get(&(SegmentType::Template, 2)).unwrap(),
        );
    }

    #[rstest]
    #[case(&"foo", &"foo", "L001_R1_001.fastq.gz", "L002_R1_001.fastq.gz", true)] // same FASTQ group, read names match
    #[case(&"foo", &"foo", "L001_R1_001.fastq.gz", "L001_R2_001.fastq.gz", true)] // different FASTQ group, read names match
    #[case(&"fo", &"foo", "L001_R1_001.fastq.gz", "L001_R2_001.fastq.gz", false)] // read names different length
    #[case(&"foo", &"fao", "L001_R1_001.fastq.gz", "L002_R1_001.fastq.gz", false)] // read names same length, mismatch
    fn test_opts_from_fastqs_check_first_read_name(
        #[case] fq1_read_name: &str,
        #[case] fq2_read_name: &str,
        #[case] fq1_suffix: &str,
        #[case] fq2_suffix: &str,
        #[case] ok: bool,
    ) {
        let dir = tempdir().unwrap();

        let fastq_r1 = dir.path().join(format!("prefix_{}", fq1_suffix));
        write_with_read_name(&fastq_r1, fq1_read_name);

        let fastq_r2 = dir.path().join(format!("prefix_{}", fq2_suffix));
        write_with_read_name(&fastq_r2, fq2_read_name);

        let opts = Opts {
            read_structures: vec![],
            fastqs: vec![dir.path().to_path_buf()],
            ..Opts::default()
        };

        let result = Opts::from(opts.fastqs, opts.read_structures, false);

        assert_eq!(result.is_ok(), ok);
        if let Err(error) = result {
            assert!(error.to_string().contains("Mismatching read names in the FASTQS"));
        }
    }

    // 2. FASTQ prefix, index FASTQs, sample_barcode_in_fastq_header: true => error

    #[test]
    fn test_opts_from_sample_barcode_in_fastq_header_is_error() {
        let dir = tempdir().unwrap();
        let fastq = dir.path().join("foo_L001_I1_001.fastq.gz");
        let prefix = dir.path().join("foo");
        let read_structure = ReadStructure::from_str("8B10S12B").unwrap();
        write(&fastq, 50);

        // FASTQ file
        // Read structure has a sample barcode, but extracting sample barcode from FASTQ header,
        // so should fail
        let result = Opts::from(vec![fastq.clone()], vec![read_structure.clone()], true);
        assert!(result.is_err());

        // FASTQ file
        // Read structure has a sample barcode and not extracting sample barcode from FASTQ header,
        // so should be ok
        let result = Opts::from(vec![fastq], vec![read_structure], false);
        assert!(result.is_ok());

        // Path prefix
        // Index FASTQ found, but extracting sample barcode from FASTQ header,
        // so should fail
        let result = Opts::from(vec![prefix.clone()], vec![], true);
        assert!(result.is_err());

        // Path prefix
        //Index FASTQ found and not extracting sample barcode from FASTQ header,
        // so should be ok
        let result = Opts::from(vec![prefix], vec![], false);
        assert!(result.is_ok());
    }
}
