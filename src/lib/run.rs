use std::{fs::File, io::BufWriter, sync::Arc, vec::Vec};

use anyhow::{bail, ensure, Context, Result};
use gzp::BUFSIZE;
use itertools::Itertools;
use log::{debug, info};
use parking_lot::Mutex;
use pooled_writer::{bgzf::BgzfCompressor, Pool};
use rayon::prelude::*;
use read_structure::ReadStructure;

use crate::{
    demux::{Demultiplex, Demultiplexer, DemuxedGroup, PerFastqRecordSet},
    matcher::{CachedHammingDistanceMatcher, MatcherKind, PreComputeMatcher},
    metrics::{DemuxedGroupMetrics, UnmatchedCounterThread},
    opts::{Opts, LOGO},
    pooled_sample_writer::PooledSampleWriter,
    sample_sheet::{self},
    thread_reader::ThreadReader,
    utils::{check_bgzf, filenames, InputFastq, MultiZip},
};

/// Run demultiplexing.
#[allow(clippy::too_many_lines)]
pub fn run(opts: Opts) -> Result<(), anyhow::Error> {
    eprint!("{}", LOGO);

    let output_types_to_write = opts.output_types_to_write()?;
    let sample_sheet = sample_sheet::SampleSheet::from_path(opts)?;
    let samples = sample_sheet.samples;
    let opts = sample_sheet.opts;

    // If there is only a single sample in the metadata and that sample has no barcode
    // specified the tool is being run to filter/mask/etc. _without_ demultiplexing
    // and so all accepted records will go into file(s) for one sample with no Undetermined
    let is_no_demux = samples.len() == 1 && samples[0].barcode.len() == 0;

    // If a path prefix is given as the single argument to --fastqs then auto-detect all the FASTQs
    // with that path prefix.  Also, if the read-structures is not given, build it from the auto-
    // detected
    let opts = if opts.fastqs.len() == 1 && !opts.fastqs[0].is_file() {
        let input_fastqs: Vec<InputFastq> = InputFastq::slurp(opts.fastqs[0].clone());
        let fastqs = input_fastqs.iter().map(|f| f.path.clone()).collect();
        let read_structures: Vec<ReadStructure> = if opts.read_structures.is_empty() {
            input_fastqs.iter().map(InputFastq::read_structure).collect::<Vec<ReadStructure>>()
        } else {
            opts.read_structures
        };
        Opts { fastqs, read_structures, ..opts }
    } else {
        opts
    };

    // Important: this must be created **after** updating the number of read structures
    let read_filter_config = opts.as_read_filter_config();

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
        opts.quality_mask_threshold.is_empty()
            || opts.quality_mask_threshold.len() == 1
            || opts.quality_mask_threshold.len() == opts.read_structures.len(),
        "Either a single quality mask threshold, or one per fastq must be provided."
    );
    ensure!(
        samples.iter().all(|s| s.barcode.len() == samples[0].barcode.len()),
        "Sample metadata barcodes are unequal lengths."
    );

    // If there is a read structure that's all sample barcode, we need to replace it with the
    // expected length to enable index hopping metrics.  Do so by inspecting the first read in the
    // corresponding FASTQ
    // TODO: move to it's own method to test
    let opts = opts.with_fixed_sample_barcodes()?;

    // All sample barcode read segments should now have a fixed length, so check the sum of their
    // lengths with the sum of length of the sample barcode(s) in the sample sheet.
    ensure!(
        is_no_demux
            || opts
                .read_structures
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
                        // Create the output directory. This is important if the output is being written
                        // to a project-specific directory
                        if let Some(parent) = name.parent() {
                            std::fs::create_dir_all(parent).with_context(|| {
                                format!("Could not create output directory: {:?}", parent)
                            })?;
                        }
                        // Now create the output file
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

#[cfg(test)]
mod test {
    use std::{
        fs::create_dir,
        path::{Path, PathBuf},
        str::FromStr,
    };

    use fgoxide::io::{DelimFile, Io};
    use itertools::Itertools;
    use read_structure::{ReadStructure, SegmentType};
    use rstest::rstest;
    use seq_io::{fastq::OwnedRecord, BaseRecord};

    use crate::{
        fastq_header::FastqHeader,
        matcher::{MatcherKind, UNDETERMINED_NAME},
        metrics::{BarcodeCount, RunMetrics, SampleMetricsProcessed},
        sample_metadata::SampleMetadata,
        sample_sheet::{SampleSheet, SampleSheetError},
        utils::{
            filename,
            test_commons::{
                create_preset_sample_metadata_file, slurp_fastq, write_reads_to_file, Fq,
                SAMPLE_BARCODE_1, SAMPLE_BARCODE_4,
            },
            INPUT_FASTQ_SUFFIX,
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
        masking_threshold_qual: u8,
        threads: usize,
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
            sample_metadata: metadata.clone(),
            read_structures: vec![read_structure],
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            quality_mask_threshold: vec![masking_threshold_qual],
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        run(sample_sheet.opts).unwrap();

        let fastq = output.join(filename(&sample_sheet.samples[0], &SegmentType::Template, 1));
        slurp_fastq(&fastq)
    }

    #[rstest]
    #[case(0, 1, & [b'A'; 39])]
    #[case(11, 1, & [[b'N'; 9].as_slice(), & [b'A'; 30]].concat())]
    #[case(0, 2, & [b'A'; 39])]
    #[case(11, 2, & [[b'N'; 9].as_slice(), & [b'A'; 30]].concat())]
    #[case(40, 2, & [[b'N'; 38].as_slice(), & [b'A'; 1]].concat())]
    #[case(100, 2, & [b'N'; 39])]
    fn test_run_end_to_end_with_quality_threshold(
        #[case] quality: u8,
        #[case] threads: usize,
        #[case] expected: &[u8],
    ) {
        // https://github.com/fulcrumgenomics/fgbio/blob/ed906cd29a7f38c9da0583c8307dc7fbcd94880b/src/test/scala/com/fulcrumgenomics/fastq/DemuxFastqsTest.scala#L470
        // Note that this tool masks anything less than or equal to the quality threshold, but fgbio masks only less than
        // so this will mask one more base than fgbio
        let records = test_end_to_end_with_quality_threshold(quality, threads);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq, expected);
    }

    #[rstest]
    #[should_panic(expected = "Same number of read structures should be given as FASTQs")]
    #[case(vec ! [ReadStructure::from_str("8B100T").unwrap(), ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(expected = "No sample barcodes found in read structures")]
    #[case(vec ! [ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(expected = "No templates found in read structures")]
    #[case(vec ! [ReadStructure::from_str("8B").unwrap()])]
    #[should_panic(expected = "Same number of read structures should be given as FASTQs")]
    #[case(vec ! [ReadStructure::from_str("8B92T").unwrap(), ReadStructure::from_str("100T").unwrap(), ReadStructure::from_str("100T").unwrap()])]
    #[should_panic(
        expected = "The number of sample barcode bases in read structures does not match sample metadata"
    )]
    #[case(vec ! [ReadStructure::from_str("18B100T").unwrap()])]
    fn test_read_structure_failures(#[case] read_structures: Vec<ReadStructure>) {
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

    fn kind_to_char(kind: SegmentType) -> char {
        match kind {
            SegmentType::SampleBarcode => 'I',
            SegmentType::Skip => 'S',
            SegmentType::MolecularBarcode => 'U',
            SegmentType::Template => 'R',
            kind => panic!("Could not determine kind from {:?}", kind),
        }
    }

    fn fastq_file_name(prefix: &str, lane: usize, kind: SegmentType, kind_number: usize) -> String {
        format!(
            "{}_L00{}_{}{}{}",
            prefix,
            lane,
            kind_to_char(kind),
            kind_number,
            INPUT_FASTQ_SUFFIX
        )
    }

    fn fastq_path(dir: impl AsRef<Path>) -> PathBuf {
        let file_name = fastq_file_name("test", 1, SegmentType::Template, 1);
        let path = dir.as_ref().join(file_name);
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
        #[values("T", "B", "TB")] output_types: String,
        #[values(true, false)] use_path_prefix: bool,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let read_structure = ReadStructure::from_str("17B100T").unwrap();
        let input = fastq_path(&dir.path());
        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());
        let fastqs: PathBuf = if use_path_prefix { dir.path().to_path_buf() } else { input };

        let opts = Opts {
            fastqs: vec![fastqs],
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures: vec![read_structure],
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            output_types,
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        let output_types_to_write = sample_sheet.opts.output_types_to_write().unwrap();

        run(sample_sheet.opts).unwrap();

        // Check outputs
        for sample in &sample_sheet.samples {
            for kind in &output_types_to_write {
                let fastq = output.join(filename(sample, kind, 1));
                let records = slurp_fastq(&fastq);
                let (names, barcodes) = names_and_barcodes(&records);

                if sample.sample_id == "Sample1" {
                    assert_eq!(records.len(), 2);
                    assert!(names.contains(&b"frag1".to_vec()));
                    assert!(names.contains(&b"frag2".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGATTACAGA".to_vec()));
                    assert!(barcodes.contains(&b"AAAAAAAAGATTACAGT".to_vec()));
                } else if sample.sample_id == UNDETERMINED_NAME {
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

        for (metric, sample) in per_sample_metrics.into_iter().zip(sample_sheet.samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == sample_sheet.samples[0].sample_id {
                assert_eq!(metric.templates, 2);
                assert_eq!(metric.perfect_matches, 1);
                assert_eq!(metric.one_mismatch_matches, 1);
            } else if metric.barcode_name == sample_sheet.samples[4].sample_id {
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
            filter_control_reads: omit_control_reads,
            filter_failing_quality: omit_failing_reads,
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        run(sample_sheet.opts).unwrap();

        for (i, sample) in sample_sheet.samples.iter().enumerate() {
            let barcodes: Vec<_> = barcodes_per_sample[i].iter().cloned().collect();
            let assignment: Vec<_> =
                assignments_per_sample[i].iter().map(|s| (*s).as_bytes()).collect();
            let r1_fastq = output.join(filename(&sample, &SegmentType::Template, 1));
            let r2_fastq = output.join(filename(&sample, &SegmentType::Template, 2));
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
    fn test_demux_fragment_reads_with_standard_qualities(#[values(1, 2)] threads: usize) {
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
        let sample =
            SampleMetadata::new(String::from("foo"), b"TCGT".as_slice().into(), 0, 2).unwrap();
        let sample_result: Result<SampleMetadata, SampleSheetError> = Ok(sample.clone());
        delim.write_csv(&metadata, sample_result).unwrap();

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
            ..Opts::default()
        };

        run(opts).unwrap();

        let record = &slurp_fastq(output.join(filename(&sample, &SegmentType::Template, 1)))[0];
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
            ..Opts::default()
        };

        run(opts).unwrap();
    }

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_demux_dual_index_paired_end_reads(
        #[values(1, 2)] threads: usize,
        #[values(true, false)] use_path_prefix: bool,
        #[values(true, false)] empty_read_structures: bool,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let fq1_path = dir.path().join(fastq_file_name("test", 1, SegmentType::SampleBarcode, 1));
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: b"AAAAAAAA", ..Fq::default() }.to_owned_record(),
            ),
            &fq1_path,
        );
        let fq2_path = dir.path().join(fastq_file_name("test", 1, SegmentType::Template, 1));
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: &[b'A'; 100], ..Fq::default() }.to_owned_record(),
            ),
            &fq2_path,
        );
        let fq3_path = dir.path().join(fastq_file_name("test", 1, SegmentType::Template, 2));
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: &[b'T'; 100], ..Fq::default() }.to_owned_record(),
            ),
            &fq3_path,
        );
        let fq4_path = dir.path().join(fastq_file_name("test", 1, SegmentType::SampleBarcode, 2));
        write_reads_to_file(
            std::iter::once(
                Fq { name: "frag", bases: b"GATTACAGA", ..Fq::default() }.to_owned_record(),
            ),
            &fq4_path,
        );

        let read_structures = if use_path_prefix && empty_read_structures {
            vec![]
        } else {
            vec![
                ReadStructure::from_str("8B").unwrap(),
                ReadStructure::from_str("100T").unwrap(),
                ReadStructure::from_str("100T").unwrap(),
                ReadStructure::from_str("9B").unwrap(),
            ]
        };

        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = create_preset_sample_metadata_file(&dir.path());

        let fastqs = if use_path_prefix {
            vec![dir.as_ref().to_path_buf()]
        } else {
            vec![fq1_path, fq2_path, fq3_path, fq4_path]
        };

        let opts = Opts {
            fastqs,
            output_dir: output.clone(),
            sample_metadata: metadata,
            read_structures,
            allowed_mismatches: 2,
            min_delta: 3,
            demux_threads: threads,
            compressor_threads: threads,
            writer_threads: threads,
            filter_control_reads: false,
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        run(sample_sheet.opts).unwrap();

        for (i, sample) in sample_sheet.samples.iter().enumerate() {
            let r1_fastq = output.join(filename(sample, &SegmentType::Template, 1));
            let r2_fastq = output.join(filename(sample, &SegmentType::Template, 2));
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
        for (metric, sample) in per_sample_metrics.into_iter().zip(sample_sheet.samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == sample_sheet.samples[0].sample_id {
                assert_eq!(metric.templates, 1);
                assert_eq!(metric.perfect_matches, 1);
                assert_eq!(metric.one_mismatch_matches, 0);
            } else if metric.barcode_name == sample_sheet.samples[4].sample_id {
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
            filter_control_reads: omit_control_reads,
            filter_failing_quality: omit_failing_reads,
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        run(sample_sheet.opts).unwrap();

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

        for (metric, sample) in per_sample_metrics.into_iter().zip(sample_sheet.samples.iter()) {
            assert_eq!(metric.barcode_name, *sample.sample_id);
            assert_eq!(metric.barcode, *sample.barcode.to_string());
            assert_eq!(metric.library_name, *sample.sample_id);

            if metric.barcode_name == sample_sheet.samples[0].sample_id {
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

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_end_to_end_with_all_matchers(
        #[values(None, Some(MatcherKind::PreCompute), Some(MatcherKind::CachedHammingDistance))]
        matcher: Option<MatcherKind>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let fq_path = dir.path().join("fq1.fastq.gz");
        let metadata_path = dir.path().join("metadata.csv");
        let output_path = dir.path().join("output");

        let reads = vec![
            Fq { name: "q1", bases: b"AAAAAAAACGACTCGTCATGA", ..Fq::default() }.to_owned_record(),
            Fq { name: "q2", bases: b"AAAATAAAATATCGCGTCTAT", ..Fq::default() }.to_owned_record(),
            Fq { name: "q3", bases: b"ACAAATAACCGTATCGGCTTA", ..Fq::default() }.to_owned_record(),
        ];

        write_reads_to_file(reads.clone().into_iter(), &fq_path);

        let read_structure = ReadStructure::from_str("8B+T").unwrap();
        let metadata = "Sample_ID,Sample_Barcode\ns1,AAAAAAAA\ns2,ACGTACGT";
        Io::default().write_lines(&metadata_path, vec![metadata]).unwrap();

        create_dir(&output_path).unwrap();

        let opts = Opts {
            fastqs: vec![fq_path],
            output_dir: output_path.clone(),
            sample_metadata: metadata_path.clone(),
            read_structures: vec![read_structure],
            allowed_mismatches: 1,
            min_delta: 2,
            demux_threads: 2,
            compressor_threads: 2,
            writer_threads: 2,
            filter_control_reads: false,
            override_matcher: matcher,
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        let samples = sample_sheet.samples;

        // Run the tool and then read back the data for s1
        run(sample_sheet.opts).unwrap();

        let s1_recs =
            slurp_fastq(&output_path.join(filename(&samples[0], &SegmentType::Template, 1)));
        let s2_recs =
            slurp_fastq(&output_path.join(filename(&samples[1], &SegmentType::Template, 1)));
        let un_recs =
            slurp_fastq(&output_path.join(filename(&samples[2], &SegmentType::Template, 1)));

        assert_eq!(s1_recs.len(), 2);
        assert_eq!(s2_recs.len(), 0);
        assert_eq!(un_recs.len(), 1);
    }

    #[rstest]
    fn test_single_sample_no_demux() {
        let dir = tempfile::tempdir().unwrap();
        let fq_path = dir.path().join("fq1.fastq.gz");
        let metadata_path = dir.path().join("metadata.csv");
        let output_path = dir.path().join("output");

        let reads = vec![
            Fq { name: "q1", bases: b"AAAAAAAACGACTCGTCATGA", ..Fq::default() }.to_owned_record(),
            Fq { name: "q2", bases: b"NNNNNNNNATATCGCGTCTAT", ..Fq::default() }.to_owned_record(),
            Fq { name: "q3", bases: b"ACAAATAACCGTATCGGCTTA", ..Fq::default() }.to_owned_record(),
        ];
        write_reads_to_file(reads.clone().into_iter(), &fq_path);

        let read_structure = ReadStructure::from_str("8B+T").unwrap();
        let metadata = "Sample_ID,Sample_Barcode\ns1,";
        Io::default().write_lines(&metadata_path, vec![metadata]).unwrap();

        create_dir(&output_path).unwrap();

        let opts = Opts {
            fastqs: vec![fq_path],
            output_dir: output_path.clone(),
            sample_metadata: metadata_path.clone(),
            read_structures: vec![read_structure],
            demux_threads: 2,
            compressor_threads: 2,
            writer_threads: 2,
            filter_control_reads: false,
            ..Opts::default()
        };

        // Run the tool and then read back the data for s1
        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        let undetermined_id = sample_sheet.opts.undetermined_sample_name.clone();
        run(sample_sheet.opts).unwrap();

        let s1_recs = slurp_fastq(&output_path.join(filename(
            &sample_sheet.samples[0],
            &SegmentType::Template,
            1,
        )));

        assert_eq!(s1_recs.len(), 3);
        for (idx, rec) in s1_recs.iter().enumerate() {
            let expected = &reads[idx].seq[8..];
            assert_eq!(rec.seq.as_slice(), expected);
        }

        // Check that there is no Undetermined file
        let undetermined = output_path.join(format!("{}_R1.fastq.gz", undetermined_id));
        assert!(!undetermined.exists(), "Undetermined file should not exist.");
    }

    #[rstest]
    #[allow(clippy::too_many_lines)]
    fn test_demux_multiple_q_masks() {
        let dir = tempfile::tempdir().unwrap();
        let fastqs = vec!["r1", "r2", "i1"]
            .into_iter()
            .map(|name| dir.path().join(format!("{}.fastq.gz", name)))
            .collect_vec();

        let read_structures = vec![
            ReadStructure::from_str("6M+T").unwrap(), // 6bp UMI then template
            ReadStructure::from_str("+T").unwrap(),   // just template
            ReadStructure::from_str("8B").unwrap(),   // just an 8bp index
        ];

        // Note we're going to filter R1 <= Q20, R2 <= Q10 and I1 <= @15 _but_ only mask template
        // segments.  Quals used in qual strings:
        //   - q30 = ?  => q30 = ?
        //   - q21 = 6  => q20 = 5
        //   - q20 = 5  => q19 = 4
        //   - q16 = 1  => q15 = 0
        //   - q15 = 0  => q14 = /
        //   - q11 = ,  => q10 = +
        //   - q10 = +  => q9  = *
        let r1 = Fq {
            name: "q1",
            bases: b"ACGTAATTTTAAAA", // u1=ACGTAA, r1=TTTNAANN
            quals: Some(b"*+/045??54?50/"),
            ..Fq::default()
        };
        let r2 = Fq { name: "q2", bases: b"GTAGCTAC", quals: Some(b"++**/045"), ..Fq::default() };
        let i1 = Fq { name: "q1", bases: b"GGTCAGAT", quals: Some(b"4444****"), ..Fq::default() };

        for (fastq, read) in fastqs.iter().zip(vec![r1, r2, i1].iter()) {
            write_reads_to_file(std::iter::once(read.to_owned_record()), fastq);
        }

        let output = dir.path().join("output");
        create_dir(&output).unwrap();

        let metadata = dir.path().join("samples.csv");
        let metadata_txt = "Sample_ID,Sample_Barcode\ns1,GGTCAGAT\n";
        std::fs::write(&metadata, metadata_txt).expect("Failed to write sample metadata.");

        let opts = Opts {
            fastqs,
            output_dir: output.clone(),
            sample_metadata: metadata.clone(),
            read_structures,
            output_types: "TBM".to_owned(),
            allowed_mismatches: 0,
            min_delta: 1,
            demux_threads: 1,
            compressor_threads: 1,
            writer_threads: 1,
            filter_control_reads: false,
            quality_mask_threshold: vec![20, 10, 15],
            ..Opts::default()
        };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        let samples = sample_sheet.samples;
        run(sample_sheet.opts).unwrap();

        // Check the undetermined sample got no reads
        let un_recs = slurp_fastq(&output.join(filename(&samples[1], &SegmentType::Template, 1)));

        assert!(un_recs.is_empty());

        // Check that our reads got masked appropriately
        let r1s = slurp_fastq(&output.join(filename(&samples[0], &SegmentType::Template, 1)));
        let r2s = slurp_fastq(&output.join(filename(&samples[0], &SegmentType::Template, 2)));
        let i1s = slurp_fastq(&output.join(filename(&samples[0], &SegmentType::SampleBarcode, 1)));
        let u1s =
            slurp_fastq(&output.join(filename(&samples[0], &SegmentType::MolecularBarcode, 1)));

        for recs in vec![&r1s, &r2s, &i1s, &u1s].into_iter() {
            assert_eq!(recs.len(), 1);
        }

        assert_eq!(std::str::from_utf8(&r1s[0].seq).unwrap(), "TTTNAANN");
        assert_eq!(std::str::from_utf8(&r2s[0].seq).unwrap(), "GTNNCTAC");
        assert_eq!(std::str::from_utf8(&i1s[0].seq).unwrap(), "GGTCAGAT"); // no masking
        assert_eq!(std::str::from_utf8(&u1s[0].seq).unwrap(), "ACGTAA"); // no masking
    }
}
