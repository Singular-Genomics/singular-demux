//! A small manager struct to orchestrate writing sets of reads (R1, R2, etc) to a set of [`PooledWriter`]s.

use anyhow::{ensure, Result};
use pooled_writer::PooledWriter;
use seq_io::BaseRecord;

use crate::demux::OutputPerSampleReads;

/// A struct that holds onto the [`PooledWriter`]s for a given sample.
#[derive(Debug)]
pub struct PooledSampleWriter {
    pub writers: Vec<PooledWriter>,
    pub reads_seen: usize,
}

impl PooledSampleWriter {
    /// Create a new [`PooledSampleWriter`].
    pub fn new(writers: Vec<PooledWriter>) -> Result<Self> {
        ensure!(!writers.is_empty(), "At least one writer must be provided");
        Ok(Self { writers, reads_seen: 0 })
    }

    /// Write the demultiplexed reads to each of their respective FASTQ files.
    #[allow(clippy::missing_panics_doc)]
    pub fn write_records(&mut self, sample_reads: OutputPerSampleReads) -> Result<()> {
        let mut seen = 0;
        for (i, reads) in sample_reads.per_fastq_reads.into_iter().enumerate() {
            seen += reads.len();
            if let Some(mut writer) = self.writers.get_mut(i) {
                for read in reads {
                    read.write(&mut writer)?;
                }
            } else {
                unreachable!();
            }
        }

        self.reads_seen += seen / self.writers.len();
        Ok(())
    }

    /// Consumes [`Self`]. For each held writer, call finish and drop to flush all writers.
    pub fn finish(self) -> Result<()> {
        for writer in self.writers {
            writer.close()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::{
        fs::File,
        io::{BufReader, BufWriter},
        result::Result,
    };

    use gzp::BgzfSyncReader;
    use itertools::Itertools;
    use pooled_writer::{bgzf::BgzfCompressor, Pool};
    use tempfile::tempdir;

    use crate::{demux::OutputPerSampleReads, utils::test_commons::generate_reads};

    use super::PooledSampleWriter;

    #[test]
    fn test_pooled_sample_writer_simple() {
        let dir = tempdir().unwrap();
        let files = vec![
            dir.path().join("s1_R1.fastq.gz"),
            dir.path().join("s1_R2.fastq.gz"),
            dir.path().join("s1_I1.fastq.gz"),
            dir.path().join("s2_R1.fastq.gz"),
            dir.path().join("s2_R2.fastq.gz"),
            dir.path().join("s2_I1.fastq.gz"),
        ];

        let raw_writers =
            files.iter().map(|name| BufWriter::new(File::create(name).unwrap())).collect();

        let (mut pool, pooled_writers) =
            Pool::new::<_, BgzfCompressor>(1, 1, 2, raw_writers).unwrap();

        let group_size = 3;
        let mut writers = vec![];
        for grouped_writers in &pooled_writers.into_iter().chunks(group_size) {
            writers.push(PooledSampleWriter::new(grouped_writers.collect::<Vec<_>>()).unwrap());
        }

        let mut reads = vec![];
        for (sample_number, read_name) in
            [(1, "R1"), (1, "R2"), (1, "I2"), (2, "R1"), (2, "R2"), (2, "I2")]
        {
            reads.push(generate_reads(sample_number, read_name, 2000));
        }

        for (counter, grouped_reads) in (&reads.iter().chunks(group_size)).into_iter().enumerate() {
            let grouped_reads = grouped_reads.collect::<Vec<_>>();
            let mut reads = OutputPerSampleReads::new(group_size);
            reads.per_fastq_reads[0].extend(grouped_reads[0].iter().cloned());
            reads.per_fastq_reads[1].extend(grouped_reads[1].iter().cloned());
            reads.per_fastq_reads[2].extend(grouped_reads[2].iter().cloned());
            writers[counter].write_records(reads).unwrap();
        }

        let reads_seen: Vec<_> = writers
            .into_iter()
            .map(|w| {
                let reads = w.reads_seen;
                w.finish().unwrap();
                reads
            })
            .collect();
        pool.stop_pool().unwrap();
        for seen in reads_seen {
            assert_eq!(seen, 2000);
        }

        for (i, file) in files.iter().enumerate() {
            let reader = BgzfSyncReader::new(BufReader::new(File::open(file).unwrap()));
            let mut reader = seq_io::fastq::Reader::new(reader);
            let records = reader.records().map(Result::unwrap).collect::<Vec<_>>();
            assert_eq!(&records, &reads[i]);
        }
    }
}
