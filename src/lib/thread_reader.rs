//! A [`ThreadReader`] allows for pushing the reading of a compressed FASTQ file onto a separate thread.
//!
//! This will spun up pooled decompressor with a default of 4 threads for decompression and will also
//! do some initial parsing of the FASTQ file into N sized chunks of FASTQ records.

use std::{fs::File, io::BufReader, num::NonZeroUsize, path::PathBuf, thread::JoinHandle};

use anyhow::{anyhow, Context, Result};
use flume::{bounded, Receiver};
use gzp::{deflate::Bgzf, par::decompress::ParDecompressBuilder, BUFSIZE};
use seq_io::fastq::{self, RecordSet};

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
    /// Create a new [`ThreadReader`] for a given file that will read `chunksize` records at a time.
    ///
    ///
    /// # Panics
    /// - The reader thread will panic if the input file can't be opened for reading.
    /// - The reader thread will panic if it is unable to send over the channel.
    // Chunksize can't be zero or `seq_io` would try to read the whole file into one record set.
    pub fn new(
        file: PathBuf,
        chunksize: NonZeroUsize,
        decompression_threads_per_reader: usize,
    ) -> Self {
        let (tx, rx) = bounded(READER_CHANNEL_SIZE);
        let handle = std::thread::spawn(move || {
            let reader = BufReader::with_capacity(
                BUFSIZE,
                File::open(&file)
                    .with_context(|| format!("Failed to open {}", file.to_string_lossy()))?,
            );

            let mut reader = fastq::Reader::with_capacity(
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

            // Developer note: read in the first record set so that a well-formatted message is
            // given when the input file is not BGZF.
            let mut record_set = RecordSet::default();
            let mut filled_set = reader
                .read_record_set_exact(&mut record_set, usize::from(chunksize))
                .map_err(|outer_error| {
                    if let seq_io::fastq::ErrorKind::Io(inner_error) = outer_error.kind() {
                        // This implies the inner_error is of kind `std::io::ErrorKind::Other`
                        if let Some(gzp::GzpError::InvalidHeader(_)) =
                            inner_error.get_ref().and_then(|i| i.downcast_ref::<gzp::GzpError>())
                        {
                            let filename = file.to_string_lossy();
                            let message = format!(
                                "
Error reading from: {}
Hint: is the input FASTQ a valid BGZF (Blocked GNU Zipped Format) and not GZIP?
Try re-compressing with bgzip (install with conda install -c bioconda htslib):
    bgzip -d {} | bgzip --stdout --threads {} > {}.bgz",
                                filename, filename, decompression_threads_per_reader, filename,
                            );
                            return anyhow!(message).context(outer_error);
                        }
                    };
                    anyhow!(outer_error)
                })?;

            while filled_set {
                tx.send(record_set).context("Failed to send record set from reader")?;
                record_set = RecordSet::default();
                filled_set = reader
                    .read_record_set_exact(&mut record_set, usize::from(chunksize))
                    .with_context(|| {
                        format!("Failed reading record set from {}", file.to_string_lossy())
                    })?;
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

        let reader = super::ThreadReader::new(file, num::NonZeroUsize::new(chunksize).unwrap(), 4);
        let seen_reads: Vec<OwnedRecord> = reader
            .rx
            .iter()
            .flat_map(|chunk| {
                chunk.into_iter().map(|record| record.to_owned_record()).collect::<Vec<_>>()
            })
            .collect();
        assert_eq!(seen_reads, reads);
    }
}
