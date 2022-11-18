//! A [`ThreadReader`] allows for pushing the reading of a compressed FASTQ file onto a separate thread.
//!
//! This will spun up pooled decompressor with a default of 4 threads for decompression and will also
//! do some initial parsing of the FASTQ file into N sized chunks of FASTQ records.

use std::{fs::File, io::BufReader, num::NonZeroUsize, path::PathBuf, thread::JoinHandle};

use anyhow::{Context, Result};
use flume::{bounded, Receiver};
use gzp::{
    deflate::Bgzf, par::decompress::ParDecompress, par::decompress::ParDecompressBuilder, BUFSIZE,
};
use seq_io::fastq::{self, Reader, RecordSet};

/// The number of chunks to allow in the reader channel at one time.
const READER_CHANNEL_SIZE: usize = 100;

/// A struct to hold onto the handle for a thread that is reading chunks of FASTQ from a file.
pub struct ThreadReader {
    /// The [`JoinHandle`] for the thread that is reading.
    pub handle: JoinHandle<Result<()>>,
    /// The channel that will be receiving [`RecordSet`]s.
    pub rx: Receiver<RecordSet>,
}

impl ThreadReader {
    /// Create a new [`ThreadReader`] that will read `chunksize` records at a time from each file
    /// serially.
    ///
    /// # Panics
    /// - The reader thread will panic if the input file can't be opened for reading.
    /// - The reader thread will panic if it is unable to send over the channel.
    // Chunksize can't be zero or `seq_io` would try to read the whole file into one record set.
    pub fn new(
        files: Vec<PathBuf>,
        chunksize: NonZeroUsize,
        decompression_threads_per_reader: usize,
    ) -> Self {
        let (tx, rx) = bounded(READER_CHANNEL_SIZE);
        let handle = std::thread::spawn(move || {
            // Important: check that all files exist right away so that a panic will occur if any
            // of the files are not found
            for file in &files {
                assert!(file.exists(), "Input FASTQ does not exist: {:?}", file);
                assert!(file.is_file(), "Input FASTQ is not a file: {:?}", file);
            }

            for file in &files {
                let reader = BufReader::with_capacity(
                    BUFSIZE,
                    File::open(&file)
                        .with_context(|| format!("Failed to open {}", file.to_string_lossy()))?,
                );

                let mut reader: Reader<ParDecompress<Bgzf>> = fastq::Reader::with_capacity(
                    ParDecompressBuilder::<Bgzf>::new()
                        .num_threads(decompression_threads_per_reader)
                        .with_context(|| {
                            format!(
                                "Error in setting threads when creating decompressor for {}",
                                file.to_string_lossy()
                            )
                        })?
                        .from_reader(reader),
                    BUFSIZE,
                );

                loop {
                    let mut record_set = RecordSet::default();

                    let filled_set = reader
                        .read_record_set_exact(&mut record_set, usize::from(chunksize))
                        .with_context(|| {
                            format!("Failed reading record set from {}", file.to_string_lossy())
                        })?;

                    if !filled_set {
                        break;
                    }
                    tx.send(record_set).context("Failed to send record set from reader")?;
                }
            }

            Ok(())
        });

        Self { handle, rx }
    }
}

#[cfg(test)]
mod test {
    use std::num;

    use rstest::rstest;
    use seq_io::fastq::OwnedRecord;
    use tempfile::tempdir;

    use crate::utils::test_commons::{generate_reads, write_reads_to_file};

    #[rstest]
    #[case(1, 10)]
    #[case(10, 1)]
    #[case(0, 100)]
    #[should_panic]
    #[case(100, 0)] // 0 is an illegal chunksize
    #[case(100, 1)]
    #[case(100_000, 500)] // 500 is the default chunksize
    fn test_thread_reader(#[case] reads_in_file: usize, #[case] chunksize: usize) {
        let dir = tempdir().unwrap();
        let file = dir.path().join("reads.fastq.gz");
        let reads = generate_reads(1, "frag", reads_in_file);
        write_reads_to_file(reads.clone().into_iter(), &file);

        let reader =
            super::ThreadReader::new(vec![file], num::NonZeroUsize::new(chunksize).unwrap(), 4);
        let seen_reads: Vec<OwnedRecord> = reader
            .rx
            .iter()
            .flat_map(|chunk| {
                chunk.into_iter().map(|record| record.to_owned_record()).collect::<Vec<_>>()
            })
            .collect();
        assert_eq!(seen_reads, reads);
    }

    #[test]
    fn test_thread_reader_multiple_fastqs() {
        let dir = tempdir().unwrap();

        let reads_one = generate_reads(1, "frag1", 100);
        let reads_two = generate_reads(2, "frag2", 100);

        let file_one = dir.path().join("reads1.fastq.gz");
        let file_two = dir.path().join("reads2.fastq.gz");
        write_reads_to_file(reads_one.clone().into_iter(), &file_one);
        write_reads_to_file(reads_two.clone().into_iter(), &file_two);

        let reader = super::ThreadReader::new(
            vec![file_one, file_two],
            num::NonZeroUsize::new(10).unwrap(),
            1,
        );
        let expected_reads: Vec<OwnedRecord> =
            reads_one.into_iter().chain(reads_two.into_iter()).collect();
        assert_eq!(expected_reads.len(), 200);
        let seen_reads: Vec<OwnedRecord> = reader
            .rx
            .iter()
            .flat_map(|chunk| {
                chunk.into_iter().map(|record| record.to_owned_record()).collect::<Vec<_>>()
            })
            .collect();
        assert_eq!(seen_reads, expected_reads);
    }
}
