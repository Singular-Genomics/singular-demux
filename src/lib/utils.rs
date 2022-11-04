//! Utility functions.
use std::{
    fs::File,
    io::{BufReader, Read},
    path::{Path, PathBuf},
};

use crate::sample_metadata::SampleMetadata;
use ahash::AHashMap;
use anyhow::{anyhow, Context};
use core::fmt::Display;
use gzp::{deflate::Bgzf, BlockFormatSpec, GzpError, BUFSIZE};
use itertools::Itertools;
use lazy_static::lazy_static;
use read_structure::{ReadStructure, SegmentType};

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

/// Determine the output filenames for sample given the [`SampleMetadata`] and the [`ReadStructure`].
pub fn filenames<P: AsRef<Path>>(
    sample: &SampleMetadata,
    output_dir: P,
    read_structures: &[ReadStructure],
    output_types_to_write: &[SegmentType],
) -> Vec<PathBuf> {
    let mut output_paths = vec![];
    let mut counter = AHashMap::new();

    for structure in read_structures.iter() {
        for kind in structure.iter().map(|segment| segment.kind).sorted().dedup() {
            if output_types_to_write.contains(&kind) {
                let type_number = counter.entry(kind).or_insert(1);
                output_paths.push(output_dir.as_ref().join(format!(
                    "{}_{}{}.fastq.gz",
                    sample.sample_id,
                    segment_kind_to_fastq_id(&kind),
                    type_number
                )));
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
        let opts = Opts::default();
        SampleSheet::from_path(&file, &opts).unwrap().samples
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

    use crate::sample_metadata::SampleMetadata;

    use super::{check_bgzf, filenames, MultiZip};
    use tempfile::tempdir;

    #[test]
    fn test_filenames_simple() {
        let sample =
            SampleMetadata::new(String::from("Sample1"), "ACTGACTG".as_bytes().to_vec().into(), 1)
                .unwrap();
        let output_dir = "/tmp";
        let read_structures: Vec<ReadStructure> = vec!["+T", "+T", "10M", "10B", "10S"]
            .into_iter()
            .map(|s| ReadStructure::from_str(s).unwrap())
            .collect();
        let output_types_to_write = vec![SegmentType::Template, SegmentType::SampleBarcode];

        let filenames = filenames(&sample, output_dir, &read_structures, &output_types_to_write);
        let expected = vec![
            PathBuf::from("/tmp/Sample1_R1.fastq.gz"),
            PathBuf::from("/tmp/Sample1_R2.fastq.gz"),
            PathBuf::from("/tmp/Sample1_I1.fastq.gz"),
        ];
        assert_eq!(filenames, expected);
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
}
