//! Utility functions.
use std::{
    fs::File,
    io::{BufReader, Read},
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::sample_metadata::SampleMetadata;
use ahash::AHashMap;
use anyhow::{anyhow, Context};
use core::fmt::Display;
use gzp::{
    deflate::Bgzf, par::decompress::ParDecompressBuilder, BlockFormatSpec, GzpError, BUFSIZE,
};
use itertools::Itertools;
use lazy_static::lazy_static;
use read_structure::{ReadStructure, SegmentType, ANY_LENGTH_STR};
use regex::Regex;
use seq_io::{fastq, BaseRecord};
use std::cmp::Ordering;

lazy_static! {
    /// Return the number of cpus as a String
    pub static ref NUM_CPU: String = num_cpus::get().to_string();
}

pub mod built_info {
    use lazy_static::lazy_static;
    include!(concat!(env!("OUT_DIR"), "/built.rs"));

    /// Get a software version string including
    ///   - Git commit hash
    ///   - Git dirty info (whether the repo had uncommitted changes)
    ///   - Cargo package version if no git info found
    fn get_software_version() -> String {
        let prefix = if let Some(s) = GIT_COMMIT_HASH {
            format!("{}-{}", PKG_VERSION, s[0..8].to_owned())
        } else {
            // This shouldn't happen
            PKG_VERSION.to_string()
        };
        let suffix = match GIT_DIRTY {
            Some(true) => "-dirty",
            _ => "",
        };
        format!("{}{}", prefix, suffix)
    }

    lazy_static! {
        /// Version of the software with git hash
        pub static ref VERSION: String = get_software_version();
    }
}

pub fn s(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).to_string()
}

/// Returns the character representation of this segment to be used in a FASTQ file name.
pub fn segment_kind_to_fastq_id(kind: &SegmentType) -> char {
    match kind {
        SegmentType::Template => 'R',
        SegmentType::SampleBarcode => 'I',
        SegmentType::MolecularBarcode => 'U',
        SegmentType::Skip => 'S',
        _ => unimplemented!(),
    }
}

/// Determine the output file name for the given sample, segment type, and type number.
pub fn filename(sample: &SampleMetadata, kind: &SegmentType, type_number: u32) -> String {
    format!(
        "{}_S{}_L00{}_{}{}_001.fastq.gz",
        sample.sample_id,
        sample.ordinal + 1,
        1,
        segment_kind_to_fastq_id(kind),
        type_number
    )
}

/// Determine the output filenames for sample given the [`SampleMetadata`] and the [`ReadStructure`].
/// If the sample has an associated project, the FASTQs will be written in project sub-directories.
pub fn filenames<P: AsRef<Path>>(
    sample: &SampleMetadata,
    output_dir: P,
    read_structures: &[ReadStructure],
    output_types_to_write: &[SegmentType],
) -> Vec<PathBuf> {
    let mut output_paths = vec![];
    let mut counter = AHashMap::new();
    let output_dir = match &sample.project {
        Some(project) => output_dir.as_ref().join(project.to_string()),
        None => output_dir.as_ref().to_path_buf(),
    };

    for structure in read_structures.iter() {
        for kind in structure.iter().map(|segment| segment.kind).sorted().dedup() {
            if output_types_to_write.contains(&kind) {
                let type_number = counter.entry(kind).or_insert(1);
                output_paths.push(output_dir.join(filename(sample, &kind, *type_number)));
                *type_number += 1;
            }
        }
    }
    output_paths
}

/// A `MultiZip` object allows for zipping over N items.
///
/// For example, if you have a `Vec` of length 10 of `Vec`s this will pull one item from
/// each of the 10 inner vecs and return a `Vec` of length 10 with those items.
///
/// This will stop iteration as soon as the first of the inner vecs runs out of items.
pub struct MultiZip<T>(Vec<T>);

impl<T> MultiZip<T> {
    /// Create a new [`MultiZip`] iterator over a `Vec` of items.
    #[must_use]
    pub fn new(items: Vec<T>) -> Self {
        Self(items)
    }
}

impl<T> Iterator for MultiZip<T>
where
    T: Iterator,
{
    type Item = Vec<T::Item>;
    fn next(&mut self) -> Option<Self::Item> {
        self.0.iter_mut().map(Iterator::next).collect()
    }
}

/// Checks if the file is a BGZF file
pub fn check_bgzf(file: &Path) -> Result<(), anyhow::Error> {
    let mut reader = match File::open(&file) {
        Ok(f) => BufReader::with_capacity(BUFSIZE, f),
        Err(error) => {
            return Err(error).with_context(|| format!("Failed to open {}", file.to_string_lossy()))
        }
    };
    let mut bytes = vec![0; Bgzf::HEADER_SIZE];
    match reader.read_exact(&mut bytes) {
        Err(error) => {
            // not enough bytes read
            let message =
                format!("Error reading from: {}\nIs it truncated?", file.to_string_lossy());
            Err(anyhow!(message).context(error))
        }
        Ok(()) => {
            if bytes[0] != 31 || bytes[1] != 139 || bytes[2] != 8 {
                // not a valid GZIP file
                report_bgzf_error(file, GzpError::InvalidHeader("Header not in GZIP format"), false)
            } else if bytes[3] & 4 != 4 || bytes[12] != b'B' || bytes[13] != b'C' {
                // non-BGZF GZIP file
                report_bgzf_error(
                    file,
                    GzpError::InvalidHeader("Header in GZIP but not BGZF format"),
                    true,
                )
            } else {
                Ok(())
            }
        }
    }
}

/// Creates an error message when the file does not look like a BGZF file.
fn report_bgzf_error<C>(file: &Path, context: C, is_gzip: bool) -> Result<(), anyhow::Error>
where
    C: Display + Send + Sync + 'static,
{
    let filename = file.to_string_lossy();
    let format = if is_gzip { "gzip" } else { "unknown" };
    let message = format!(
        "
Error reading from: {}

The input must be in BGZF (bgzip) format!

The input was found in a {} format.

To re-compress a GZIP file with bgzip:
  1. install with `conda install -c bioconda htslib`
     or from http://www.htslib.org/download/
  2. `gunzip -c {} > tmp.fastq`
  3. `bgzip --stdout --threads tmp.fastq > {}`

To compress an uncompressed FASTQ file with bgzip:
  1. install with `conda install -c bioconda htslib`
     or from http://www.htslib.org/download/
  2. `bgzip --threads {}`
",
        filename, format, filename, filename, filename,
    );
    Err(anyhow!(message).context(context))
}

pub static INPUT_FASTQ_SUFFIX: &str = "_001.fastq.gz";

lazy_static! {
    /// <*>_L00#_<R# or I#>_001.fastq.gz
    static ref INPUT_FASTQ_REGEX: Regex = Regex::new(r"^(.*)_L00(\d{1})_([RIUS])(\d{1})_001.fastq.gz$").unwrap();
}

/// Contains information about a FASTQ that has been inferred from the file name.  This includes
/// lane, segment type (e.g. template/read, sample barcode/index, molecular barcode/index, and
/// skip bases), and segment number (e.g. R1 or R2).
#[derive(Debug, Clone)]
pub struct InputFastq {
    pub path: PathBuf,
    pub prefix: String,
    pub lane: usize,
    pub kind: SegmentType,
    pub kind_number: i32,
}

impl InputFastq {
    /// Create a new `InputFastq` inferring information from the file name.  This must match the
    /// `INPUT_FASTQ_REGEX` pattern.
    pub fn new<P: AsRef<Path>>(path: P) -> Option<InputFastq> {
        let file_name = path.as_ref().file_name().unwrap().to_string_lossy();
        match INPUT_FASTQ_REGEX.captures(&file_name) {
            None => None,
            Some(captures) => {
                let prefix = captures.get(1).unwrap().as_str().to_string();
                let lane = captures.get(2).unwrap().as_str().parse::<usize>().unwrap();
                let kind: SegmentType = match captures.get(3).unwrap().as_str() {
                    "R" => SegmentType::Template,
                    "I" => SegmentType::SampleBarcode,
                    "U" => SegmentType::MolecularBarcode,
                    "S" => SegmentType::Skip,
                    knd => panic!("Could not determine kind from {}", knd),
                };
                let kind_number = captures.get(4).unwrap().as_str().parse::<i32>().unwrap();

                Some(InputFastq {
                    path: path.as_ref().to_path_buf(),
                    prefix,
                    lane,
                    kind,
                    kind_number,
                })
            }
        }
    }

    /// Identifies all FASTQs that share the common path prefix and match the
    /// `INPUT_FASTQ_REGEX` pattern.  The path prefix may also be a directory.  The FASTQS are
    /// returned in sorted order.
    pub fn slurp<P: AsRef<Path>>(path_prefix: P) -> Vec<InputFastq> {
        let (parent, prefix) = if path_prefix.as_ref().is_file() {
            (
                path_prefix.as_ref().parent().unwrap(),
                path_prefix.as_ref().file_name().unwrap().to_str().unwrap(),
            )
        } else {
            (path_prefix.as_ref(), "")
        };

        let mut fastqs: Vec<InputFastq> = std::fs::read_dir(parent)
            .unwrap()
            .map(|res| res.map(|e| e.path()).unwrap())
            .filter(|p| {
                let file_name = p.file_name().unwrap().to_str().unwrap();
                p.is_file()
                    && file_name.starts_with(prefix)
                    && file_name.ends_with(INPUT_FASTQ_SUFFIX)
            })
            .filter_map(InputFastq::new)
            .collect();
        fastqs.sort();
        fastqs
    }

    pub fn read_structure(&self) -> ReadStructure {
        ReadStructure::from_str(&format!("{}{}", ANY_LENGTH_STR, self.kind.value())).unwrap()
    }
}

/// Returns a numeric value in priority order for segment types, to aid in ordering input FASTQs.
///
/// The main motivation is so when the FASTQs are ordered, that for a dual-index run the FASTQs are
/// ordered assuming a read structure of `+B +T +T +B`.
fn kind_to_order_num(kind: SegmentType) -> usize {
    match kind {
        SegmentType::SampleBarcode => 0,
        SegmentType::Skip => 1,
        SegmentType::MolecularBarcode => 2,
        SegmentType::Template => 3,
        kind => panic!("Could not determine kind from {:?}", kind),
    }
}

impl Ord for InputFastq {
    /// Defines a partial ordering based on:
    /// 1. The FASTQ prefix
    /// 2. The lane
    /// 3. The read number (kind number)
    /// 4. The kind.  If the kind number is odd, orders by sample barcode, skip, molecular barcode,
    ///    then template.  The opposite order if the kind is even.  This is so if the four FASTQS
    ///    for I1, R1, R2, I2, they are remain in that order after sorting.
    ///
    /// The main motivation is so when the FASTQs are ordered, that for a dual-index run the FASTQs are
    /// ordered assuming a read structure of `+B +T +T +B`.
    ///
    /// Read numbers beyond two are supported, but are non-sensical at this time.
    fn cmp(&self, other: &Self) -> Ordering {
        // Prefix
        let mut res = self.prefix.cmp(&other.prefix);
        if res.is_ne() {
            return res;
        };

        // Lane
        res = self.lane.cmp(&other.lane);
        if res.is_ne() {
            return res;
        };

        // Kind number
        res = self.kind_number.cmp(&other.kind_number);
        if res.is_ne() {
            return res;
        };

        // Kind, conditional on kind number
        res = kind_to_order_num(self.kind).cmp(&kind_to_order_num(other.kind));
        if self.kind_number % 2 == 0 {
            res.reverse()
        } else {
            res
        }
    }
}

impl PartialOrd for InputFastq {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for InputFastq {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for InputFastq {}

/// Infers the read length contained in the given FASTQ by examining the length of the first read.
pub fn infer_fastq_sequence_length(file: PathBuf) -> Result<usize, anyhow::Error> {
    let reader = BufReader::with_capacity(
        BUFSIZE,
        File::open(&file).with_context(|| format!("Failed to open {}", file.to_string_lossy()))?,
    );

    let mut reader = fastq::Reader::with_capacity(
        ParDecompressBuilder::<Bgzf>::new()
            .num_threads(1)
            .with_context(|| {
                format!(
                    "Error in setting threads when creating decompressor for {}",
                    file.to_string_lossy()
                )
            })?
            .from_reader(reader),
        BUFSIZE,
    );

    match reader.next() {
        Some(Ok(record)) => Ok(record.seq().len()),
        _ => Err(anyhow!(
            "Could determine sample barcode length from input FASTQ: {}",
            file.to_string_lossy()
        )),
    }
}

#[cfg(not(tarpaulin_include))]
#[cfg(test)]
pub mod test_commons {
    //! Common utility methods for testing demultiplexing.

    use std::{
        fs::File,
        io::{BufReader, BufWriter, Write},
        path::{Path, PathBuf},
    };

    use crate::{opts::Opts, sample_metadata, sample_sheet::SampleSheet};
    use bgzf::CompressionLevel;
    use sample_metadata::SampleMetadata;
    use seq_io::{
        fastq::{OwnedRecord, Reader, RecordSet},
        BaseRecord,
    };
    use tempfile::{tempdir, NamedTempFile};

    pub const SAMPLE_BARCODE_1: &[u8] = b"AAAAAAAAGATTACAGA";
    pub const SAMPLE_BARCODE_2: &[u8] = b"CCCCCCCCGATTACAGA";
    pub const SAMPLE_BARCODE_3: &[u8] = b"GGGGGGGGGATTACAGA";
    pub const SAMPLE_BARCODE_4: &[u8] = b"GGGGGGTTGATTACAGA";

    /// Helper method to write a standard 4 samples to a metadata file.
    pub fn create_preset_sample_metadata_file(dir: impl AsRef<Path>) -> PathBuf {
        let file_contents = format!(
            "Sample_ID,Sample_Barcode\n\
        Sample1,{}\n\
        Sample2,{}\n\
        Sample3,{}\n\
        Sample4,{}",
            String::from_utf8_lossy(SAMPLE_BARCODE_1),
            String::from_utf8_lossy(SAMPLE_BARCODE_2),
            String::from_utf8_lossy(SAMPLE_BARCODE_3),
            String::from_utf8_lossy(SAMPLE_BARCODE_4)
        );

        let output = dir.as_ref().join("sample_metadata.csv");
        std::fs::write(&output, file_contents).expect("Failed to write sample metadata to file.");
        output
    }

    /// Create a collection of the preset sample metadata WITH an undetermined sample added on.
    pub fn create_preset_sample_metadata() -> Vec<SampleMetadata> {
        let dir = tempdir().unwrap();
        let file = create_preset_sample_metadata_file(dir.path());
        let opts = Opts { sample_metadata: file, ..Opts::default() };
        SampleSheet::from_path(opts).unwrap().samples
    }

    /// Configuration struct for creating a FASTQ read
    #[derive(Debug, Default, Clone, Copy)]
    pub struct Fq<'a> {
        pub name: &'a str,
        pub bases: &'a [u8],
        pub set_bc_to: Option<&'a [u8]>,
        pub quals: Option<&'a [u8]>,
        pub comment: Option<&'a str>,
        pub read_number: Option<usize>,
        pub is_filtered: Option<&'a str>,
        pub control_number: Option<&'a str>,
    }

    impl<'a> Fq<'a> {
        /// Convert the configuration into an [`OwnedRecord`].
        pub fn to_owned_record(&self) -> OwnedRecord {
            if self.comment.is_some()
                && (self.read_number.is_some()
                    || self.is_filtered.is_some()
                    || self.control_number.is_some())
            {
                panic!("Set comment fields in comment directly in comment");
            }

            let default_read_name = "H00233:4:AAAFGW3HV:1:1101:59586:1000";
            let readname = if let Some(comment) = self.comment {
                format!("{} {} {}", default_read_name, comment, self.name)
            } else {
                format!(
                    "{} {}:{}:{}:{} {}",
                    default_read_name,
                    self.read_number.unwrap_or(0),
                    self.is_filtered.unwrap_or("N"),
                    self.control_number.unwrap_or("0"),
                    String::from_utf8_lossy(self.set_bc_to.unwrap_or(b"NNNNNNNN")),
                    self.name
                )
            };

            let qual = if let Some(qual) = self.quals {
                assert_eq!(qual.len(), self.bases.len());
                qual.to_vec()
            } else {
                vec![b'!'; self.bases.len()]
            };

            OwnedRecord { head: readname.as_bytes().to_vec(), seq: self.bases.to_vec(), qual }
        }
    }

    /// Write a set of fastq reads to a file, returning the number of reads written.
    ///
    /// If the file extension is `gz` the reads will be compressed
    pub fn write_reads_to_file(
        reads: impl Iterator<Item = OwnedRecord>,
        file: impl AsRef<Path>,
    ) -> usize {
        let mut num_written = 0;
        let mut writer: Box<dyn Write> =
            if file.as_ref().extension().map_or(false, |ext| ext == "gz") {
                Box::new(bgzf::Writer::new(
                    BufWriter::new(File::create(file).unwrap()),
                    CompressionLevel::new(3).unwrap(),
                ))
            } else {
                Box::new(BufWriter::new(File::create(file).unwrap()))
            };
        for read in reads {
            read.write(&mut writer).unwrap();
            num_written += 1;
        }
        writer.flush().unwrap();
        num_written
    }

    /// Convert a collection of [`OwnedRecord`]s into a [`RecordSet`].
    pub fn reads_to_record_set(reads: impl Iterator<Item = OwnedRecord>) -> RecordSet {
        let file = NamedTempFile::new().unwrap();
        let num_written = write_reads_to_file(reads, file.path());
        let mut reader = Reader::from_path(file.path()).unwrap();
        let mut record_set = RecordSet::default();
        reader.read_record_set_exact(&mut record_set, num_written).unwrap();
        record_set
    }

    /// Slurp all records out of a BGZF FASTQ file
    pub fn slurp_fastq(file: impl AsRef<Path>) -> Vec<OwnedRecord> {
        let reader = bgzf::Reader::new(BufReader::new(
            File::open(&file).unwrap_or_else(|_| panic!("Unable to open {:?}", &file.as_ref())),
        ));
        let mut reader = Reader::new(reader);
        let mut records = vec![];
        for r in reader.records() {
            let r = r.unwrap();
            records.push(r);
        }
        records
    }

    /// A more generic way to create reads that allows for unique name for each read
    pub fn generate_reads(
        sample_num: usize,
        read_name: &str,
        num_reads: usize,
    ) -> Vec<OwnedRecord> {
        let mut reads = vec![];
        for i in 1..=num_reads {
            let name = format!("s{}_read_{}_of_{}_{}", sample_num, i, num_reads, read_name);
            reads.push(OwnedRecord {
                head: name.as_bytes().to_vec(),
                seq: vec![b'A'; 150],
                qual: vec![b'!'; 150],
            });
        }
        reads
    }
}

#[cfg(test)]
mod test {
    use std::{
        fs::File,
        io::{BufWriter, Write},
        iter::IntoIterator,
        path::PathBuf,
        str::FromStr,
        vec::IntoIter,
    };

    use gzp::{BgzfSyncWriter, Compression, MgzipSyncWriter};
    use read_structure::{ReadStructure, SegmentType};

    use crate::{sample_metadata::SampleMetadata, utils::filename};
    use strum::IntoEnumIterator;

    use crate::utils::{kind_to_order_num, segment_kind_to_fastq_id, INPUT_FASTQ_SUFFIX};

    use super::{check_bgzf, filenames, infer_fastq_sequence_length, InputFastq, MultiZip};
    use tempfile::tempdir;

    #[test]
    fn test_filenames_simple() {
        let sample = SampleMetadata::new(
            String::from("Sample1"),
            "ACTGACTG".as_bytes().to_vec().into(),
            1,
            2,
        )
        .unwrap();
        let output_dir = "/tmp";
        let read_structures: Vec<ReadStructure> = vec!["+T", "+T", "10M", "10B", "10S"]
            .into_iter()
            .map(|s| ReadStructure::from_str(s).unwrap())
            .collect();
        let output_types_to_write = vec![SegmentType::Template, SegmentType::SampleBarcode];

        let filenames = filenames(&sample, output_dir, &read_structures, &output_types_to_write);
        let expected = vec![
            PathBuf::from("/tmp").join(filename(&sample, &SegmentType::Template, 1)),
            PathBuf::from("/tmp").join(filename(&sample, &SegmentType::Template, 2)),
            PathBuf::from("/tmp").join(filename(&sample, &SegmentType::SampleBarcode, 1)),
        ];
        assert_eq!(filenames, expected);
    }

    #[test]
    fn test_filenames_with_project() {
        let output_dir = "/tmp";
        let read_structures: Vec<ReadStructure> = vec!["+T", "+T", "10M", "10B", "10S"]
            .into_iter()
            .map(|s| ReadStructure::from_str(s).unwrap())
            .collect();
        let output_types_to_write = vec![SegmentType::Template, SegmentType::SampleBarcode];
        let mut sample = SampleMetadata::new(
            String::from("Sample1"),
            "ACTGACTG".as_bytes().to_vec().into(),
            1,
            2,
        )
        .unwrap();

        // Sample with project
        sample.project = Some(String::from("Project1"));
        let output_fastqs =
            filenames(&sample, output_dir, &read_structures, &output_types_to_write);
        let expected = vec![
            PathBuf::from("/tmp/Project1/Sample1_S2_L001_R1_001.fastq.gz"),
            PathBuf::from("/tmp/Project1/Sample1_S2_L001_R2_001.fastq.gz"),
            PathBuf::from("/tmp/Project1/Sample1_S2_L001_I1_001.fastq.gz"),
        ];
        assert_eq!(output_fastqs, expected);

        // Sample without project
        sample.project = None;
        let output_fastqs =
            filenames(&sample, output_dir, &read_structures, &output_types_to_write);
        let expected = vec![
            PathBuf::from("/tmp/Sample1_S2_L001_R1_001.fastq.gz"),
            PathBuf::from("/tmp/Sample1_S2_L001_R2_001.fastq.gz"),
            PathBuf::from("/tmp/Sample1_S2_L001_I1_001.fastq.gz"),
        ];
        assert_eq!(output_fastqs, expected);
    }

    #[test]
    fn test_multizip_equal() {
        let input = vec![vec![1; 10], vec![1; 10], vec![1; 10], vec![1; 10]];
        let multi: MultiZip<IntoIter<i32>> =
            MultiZip::new(input.into_iter().map(IntoIterator::into_iter).collect());
        let sum = multi.map(|zipped| zipped.into_iter().sum::<i32>()).sum::<i32>();
        assert_eq!(40, sum);
    }

    #[test]
    fn test_multizip_unequal() {
        let input = vec![vec![1; 10], vec![1; 8], vec![1; 10], vec![1; 10]];
        let multi: MultiZip<IntoIter<i32>> =
            MultiZip::new(input.into_iter().map(IntoIterator::into_iter).collect());
        let sum = multi.map(|zipped| zipped.into_iter().sum::<i32>()).sum::<i32>();
        assert_eq!(32, sum);
    }

    #[test]
    fn test_check_bgzf_path_does_not_exist_fail() {
        let dir = tempdir().unwrap();
        let file = dir.path().join("does_not_exist.txt");
        assert!(check_bgzf(&file).is_err());
    }

    #[test]
    fn test_check_bgzf_on_plaintext_fail() {
        let dir = tempdir().unwrap();

        // Create output file
        let file = dir.path().join("plaintext.txt");
        let mut writer = BufWriter::new(File::create(&file).unwrap());
        writer.write_all(b"@NAME\nGATTACA\n+\nIIIIIII\n").unwrap();
        drop(writer);

        assert!(check_bgzf(&file).is_err());
    }

    #[test]
    fn test_check_bgzf_on_gzip_fail() {
        let dir = tempdir().unwrap();

        // Create output file
        let file = dir.path().join("fastq.gz");
        let writer = BufWriter::new(File::create(&file).unwrap());
        let mut gz_writer = MgzipSyncWriter::new(writer, Compression::new(3));
        gz_writer.write_all(b"@NAME\nGATTACA\n+\nIIIIIII\n").unwrap();
        gz_writer.flush().unwrap();
        drop(gz_writer);

        assert!(check_bgzf(&file).is_err());
    }

    #[test]
    fn test_check_bgzf_ok() {
        let dir = tempdir().unwrap();

        // Create output file
        let file = dir.path().join("fastq.gz");
        let writer = BufWriter::new(File::create(&file).unwrap());
        let mut gz_writer = BgzfSyncWriter::new(writer, Compression::new(3));
        gz_writer.write_all(b"@NAME\nGATTACA\n+\nIIIIIII\n").unwrap();
        gz_writer.flush().unwrap();
        drop(gz_writer);

        assert!(check_bgzf(&file).is_ok());
    }

    #[test]
    fn test_input_fastqs_does_not_match() {
        let dir = tempdir().unwrap();

        // generic
        assert!(InputFastq::new(dir.path().join("fastq.gz")).is_none());

        // lane issues
        assert!(InputFastq::new(dir.path().join("foo_R1_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_R1_L000_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_R1_L01_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_R1_L0001_001.fastq.gz")).is_none());

        // segment type issues
        assert!(InputFastq::new(dir.path().join("foo_L001_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_L001_T1_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_L001_1_001.fastq.gz")).is_none());

        // segment num
        assert!(InputFastq::new(dir.path().join("foo_L001_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_L001_R01_001.fastq.gz")).is_none());
        assert!(InputFastq::new(dir.path().join("foo_L001_RB_001.fastq.gz")).is_none());
    }

    #[test]
    fn test_input_fastqs_matches() {
        let dir = tempdir().unwrap();

        // create the FASZTQ paths
        let mut expected_fastqs = vec![];
        for prefix in ["foo", "bar"] {
            for lane in 1..4 {
                for kind_number in 1..10 {
                    for kind in SegmentType::iter() {
                        let name = format!(
                            "{}_L00{}_{}{}{}",
                            prefix.to_string(),
                            lane,
                            segment_kind_to_fastq_id(&kind),
                            kind_number,
                            INPUT_FASTQ_SUFFIX
                        );
                        let path = dir.path().join(name);
                        let fq = InputFastq::new(path.clone()).unwrap();
                        assert_eq!(
                            fq,
                            InputFastq {
                                path: path.clone(),
                                prefix: prefix.to_string(),
                                lane,
                                kind,
                                kind_number
                            }
                        );
                        // actually touch the path
                        std::fs::File::create(path.clone()).unwrap();
                        expected_fastqs.push(path);
                    }
                }
            }
        }

        // list the FASTQs from the directory
        let actual_paths: Vec<PathBuf> =
            InputFastq::slurp(dir.path()).iter().map(|p| p.path.clone()).collect();
        // check that all the expected FASTQs were found
        for expected_fastq in &expected_fastqs {
            assert!(actual_paths.contains(expected_fastq));
        }
        // check that all the found FASTQS were expected
        for actual_path in actual_paths {
            assert!(expected_fastqs.contains(&actual_path));
        }
    }

    #[test]
    fn test_input_fastqs_ordering() {
        let dir = tempdir().unwrap();

        // different prefix
        let lower = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            lane: 1,
            kind: SegmentType::SampleBarcode,
            kind_number: 1,
        };
        let upper =
            InputFastq { path: dir.path().to_path_buf(), prefix: "aab".to_string(), ..lower };
        assert!(lower.lt(&upper));
        assert!(upper.gt(&lower));
        assert!(lower.eq(&lower));

        // same prefix, different lane
        let upper = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            lane: 2,
            ..lower
        };
        assert!(lower.lt(&upper));
        assert!(upper.gt(&lower));
        assert!(lower.eq(&lower));

        // same prefix and lane, different kind number
        let upper = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            kind_number: 2,
            ..lower
        };
        assert!(lower.lt(&upper));
        assert!(upper.gt(&lower));
        assert!(lower.eq(&lower));

        // same prefix, lane, and kind number (odd), but different kind
        // For read 1, we have SampelBarcode < Skip < MolecularBarcode < Template
        let upper = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            kind: SegmentType::Template,
            ..lower
        };
        assert!(lower.lt(&upper));
        assert!(upper.gt(&lower));
        assert!(lower.eq(&lower));

        // same prefix, lane, and kind number (even), but different kind
        // For read 2, we have SampelBarcode > Skip > MolecularBarcode > Template
        let upper = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            kind_number: 2,
            kind: SegmentType::SampleBarcode,
            ..lower
        };
        let lower = InputFastq {
            path: dir.path().to_path_buf(),
            prefix: "aaa".to_string(),
            kind: SegmentType::Template,
            ..upper
        };
        assert!(lower.lt(&upper));
        assert!(upper.gt(&lower));
        assert!(lower.eq(&lower));
    }

    #[test]
    fn test_kind_to_order_num() {
        // check regression
        assert_eq!(kind_to_order_num(SegmentType::SampleBarcode), 0);
        assert_eq!(kind_to_order_num(SegmentType::Skip), 1);
        assert_eq!(kind_to_order_num(SegmentType::MolecularBarcode), 2);
        assert_eq!(kind_to_order_num(SegmentType::Template), 3);
    }

    #[test]
    fn test_infer_fastq_sequence_length() {
        let dir = tempdir().unwrap();

        // Create output file
        let file = dir.path().join("fastq.gz");
        let writer = BufWriter::new(File::create(&file).unwrap());
        let mut gz_writer = BgzfSyncWriter::new(writer, Compression::new(3));
        gz_writer.write_all(b"@NAME\nGATTACA\n+\nIIIIIII\n").unwrap();
        gz_writer.write_all(b"@NAME\nGATTACAA\n+\nIIIIIIII\n").unwrap();
        gz_writer.flush().unwrap();
        drop(gz_writer);

        // infers from the first read if only one read
        assert_eq!(infer_fastq_sequence_length(file).unwrap(), 7);
    }
}
