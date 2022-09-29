# sgdemux

[![Check](https://github.com/Singular-Genomics/singular-demux/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/Singular-Genomics/singular-demux/actions/workflows/build_and_test.yml)

This repository is home to the `sgdemux` tool for demultiplexing sequencing data generated on Singular Genomics' sequencing instruments.

## Installation

See the separate [installation doc](docs/02_Installation.md) for instructions on how to install `sgdemux`.

## Usage

`sgdemux` performs sample demultiplexing on block-compressed (BGZF) FASTQs such as those produced by the Singular Genomics G4 platform.
The input FASTQs _must_ be block compressed (e.g. with [`bgzip`](http://www.htslib.org/doc/bgzip.html)); uncompressed or non-bgzf gzipped input files are not supported as performance would be significantly degraded.

The primary inputs to the tool are:

1. A set of undemultiplexed FASTQ files (BGZF compressed)
2. A set of read-structures, one per input FASTQ file
3. A file of sample metadata including sample names and barcode sequences
4. A directory into which the demultiplexed FASTQ files should be written

Reads are written to per-sample, per-instrument-read files within the output directory.  An additional `Undetermined` set of files will be written containing those reads that did not match any expected barcodes.

An example invocation follows:

```shell
sgdemux \
  --fastqs R1.fastq.gz r2.fastq.gz I1.fastq.gz I2.fastq.gz \
  --read-structures +T +T 8B 8B \
  --sample-metadata sample-metadata.csv \
  --output-dir demuxed/
```

### Inputs

##### FASTQ Files

The full set of FASTQ files generated for a run, or lane, or sequencing should be provided, including all template and index reads.  For example if a paired-end sequencing run was performed with dual sample index reads, four files should be provided:

```shell
  --fastqs R1.fastq.gz r2.fastq.gz I1.fastq.gz I2.fastq.gz
```

If multiple FASTQ files are available per instrument reads, they should be concatenated prior to running `sgdemux`, e.g.:

```shell
for read in R1 R2 I1 I2; do cat L*/${read}.fastq.gz > ./${read}.fastq.gz; done
```

FASTQ files _must_ be BGZF compressed.

##### Read Structures

Read Structures are short strings that describe the origin and/or purpose of bases within sequencing reads.  They are made up of a sequence of `<number><operator>` pairs.  Four kinds of operators are recognized:

1. **T** identifies template reads/bases
2. **B** identifies sample barcode reads/bases
3. **M** identifies unique molecular index reads/bases
4. **S** identifies a set of bases to be skipped or ignored

The last `<number><operator>` pair in a Read Structure may use `+` instead of a number to denote "all remaining bases".
This is useful if, e.g., FASTQs have been trimmed and/or contain reads of varying length. 

For more dteails on Read Structures, and how to validate them, see [this detailed description](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures).

One Read Structure must be provided for each input FASTQ file, in the same order.  Matching the set of reads specified in the FASTQ files section above one might specify:

```shell
  --read-structures +T +T 8B 8B
```

##### Sample Metadata

The sample-metadata file specified via the `--sample-metadata` option is a simple 2-column CSV file.  An example follows:

```text
Sample_ID,Sample_Barcode
s1,ACTGGTCA
s2,ATACGAAC
```

Both `Sample_ID`s and `Sample_Barcode`s must be unique within the file, and both columns are required for all samples.

If multiple sample barcodes are are present (e.g. dual-indexing runs, additional inline sample indices) then the `Sample_Barcode` field should contain the full set of barcode bases expected to be read for the sample.  Characters other than A, C, G and T may be present in the `Sample_Barcode` field, but will be ignored.  The ordering of the concatenated barcodes is important, and should match the ordering of the FASTQs and Read Structures given.

For example if a dual-indexing run was performed with an additional inline sample barcode in read 1, and `sgdemux` was invoked with the following options:

```shell
--fastqs R1.fastq.gz I1.fastq.gz I2.fastq.gz R2.fastq.gz \
--read-structures 10B+T 8B 8B +T
```

then the `Sample_Barcode` field for each sample should be composed as follows:

```shell
  {10 base inline index}-{8 base I1 index}-{8 base I2 index}
```

##### Full Argument List

|Argument Name|Required|Default Value|Description|
|-------------|--------|-------------|-----------|
|--fastqs                  |Yes|n/a         |Path(s) to the input FASTQs to be demultiplexed.|
|--sample-metadata         |Yes|n/a         |Path to CSV of sample metadata with sample IDs and barcode sequences.|
|--read-structures         |Yes|n/a         |One read structure per input FASTQ describing how to parse the reads.|
|--output-dir              |Yes|n/a         |Path to an output directory to write into. Directory must exist.|
|--filter-control-reads    |No |False       |If true, filter out reads marked as control reads in their FASTQ headers.|
|--filter-failing-quality  |No |False       |If true, filter out reads marked as failing quality control in their FASTQ headers.|
|--allowed-mismatches      |No |1           |The number of mismatches allowed, in total, between expected and observed barcode bases in order to match a read to a sample.|
|--min-delta               |No |2           |The minimum number of mismatches by which the best match for a read is better than the next-best match for a read in order to accept the best match.|
|--free-ns                 |No |1           |The number of observed Ns (no-calls) in the barcode read(s) that are allowed for "free" before treating subsequent Ns as mismatches.|
|--max-no-calls            |No |n/a         |If specified, do not match any reads whose barcode reads contain more than this many Ns.|
|--quality-mask-threshold  |No |0           |Mask to N all bases in all reads whose base quality is at or below the specified value.|
|--output-types            |No |T           |The types of bases/reads for which output files should be generated. A single string containing one or more of `T` (template), `B` (sample barcode), `M` (UMI), and `S` (skipped).|
|--undetermined-sample-name|No |Undetermined|The name used as a prefix to generate FASTQ files for reads that didn't match to any sample.|
|--compressor-threads      |No |12          |The number of threads to use in compressing the output FASTQ files.|
|--deumux-threads          |No |4           |The number of threads to use to perform demultiplexing in memory.|
|--writer-threads          |No |5           |The number of threads to use to write compressed FASTQ data to disk.|
|--override-matcher        |No |n/a         |The algorithm to use for matching, either `CachedHammingDistance` or `PreCompute`.  By default if barcodes are 12bp or shorter `PreCompute` is used which pre-computes all possible matches, or if barcodes are longer than 12bp `CachedHammingDistance` is used which calculates matches when needed then caches the results.|
|--most-unmatched-to-output|No |1000        |Report on the top N most frequently observed unmatched barcode sequences.|


##### Performance Considerations

Various `--*-threads` options are available to control the number of threads used by `sgdemux` for various purposes.  The defaults are intended to fully utilize a 32-core machine.  The defaults to the available options do not add up to 32 as several threads are used to read the input FASTQ files and for ancillary purposes.

For running on larger or smaller instances it is advised to start with the following and tune from there:

- 1/3 of available threads for compression
- 1/6 of available threads for writing
- 1/6-1/3 of available threads for demultiplexing

Currently this tool does not provide a way place a hard limit on the number of threads used.

### Outputs

##### Demultiplexed FASTQs

One or more BGZF compressed FASTQ files will be created per sample in the specified output directory. For
paired end data, the output will have the suffix `_R1.fastq.gz` and `_R2.fastq.gz` for read one and read two
respectively.

Samples barcodes, and unique molecular indices (UMIs), will be inserted into the FASTQ headers if present.  If either multiple sample barcodes or multiple UMIs are present they will be concatenated with `+` between individual barcodes prior to insertion.  For example if a FASTQ record had sample barcodes `ACGT` and `TTGA`, and UMIs of `ACCTAG` and `TCCTGG` the the output header might look like:

```shell
sg001:17:A30ZZ:1:4:1234:4567:ACCTAG+TCCTGG 1:N:1:ACGT+TTGA
```

##### Metrics

Up to four metrics files are generated to help assess run and demultiplexing quality:

##### `per_sample_metrics.tsv` 

This file is always produced and contains the following columns:

|Column|Description|
|------|-----------|
|`barcode_name`|The name for the sample barcode, typically the same name from the SampleSheet.|
|`library_name`|The name of the library, typically the library identifier from the SampleSheet.|
|`barcode`|The sample barcode bases. Dual index barcodes will have two sample barcode sequences delimited by a `+`.|
|`templates`|The total number of templates matching the given barcode.|
|`perfect_matches`|The number of templates that match perfectly the given barcode.|
|`one_mismatch_matches`|The number of pass-filter templates that match the given barcode with exactly one mismatch.|
|`q20_bases`|The number of bases in a template with a quality score 20 or above.|
|`q30_bases`|The number of bases in a template with a quality score 30 or above.|
|`total_number_of_bases`|The total number of bases in the templates combined.|
|`fraction_matches`|The fraction of all templates that match the given barcode.|
|`ratio_this_barcode_to_best_barcode`|The ratio of templates for this barcode to the number of templates of the most prevelant barcode (excluding Undetermined).|
|`frac_q20_bases`|The fraction of bases in a template with a quality score 20 or above.|
|`frac_q30_bases`|The fraction of bases in a template with a quality score 30 or above.|
|`mean_index_base_quality`|The mean quality of index bases.|

##### `run_metrics.tsv`

This file is always produced and contains a small number of summary statistics across the demultiplexing run:

|Column|Description|
|------|-----------|
|`control_reads_omitted`|The number of reads that were omitted for being control reads.|
|`failing_reads_omitted`|The number of reads that were omitted for having failed QC.|
|`total_templates`|The total number of template reads that were output.|

##### `most_frequent_unmatched.tsv`

This file is optional and will only be produced if `--most-unmatched-to-output` is not set to zero. It contains the (approximate) counts of the most prevelant observed barcode sequences that did not match to one of the expected barcodes.

|Column|Description|
|------|-----------|
|`barcode`|The observed barcode sequence.|
|`count`|The approximate number of times that barcode sequences was observed.|

##### `sample_barcode_hop_metrics.tsv`

This file is only output for dual-indexed runs.  It contains frequently observed barcodes that are unexpected combinations of expected barcodes.  For example if two samples are present with barcodes `AA-CC` and `GG-TT`, this file would report on observations of `AA-TT` and `GG-CC` if seen.

|Column|Description|
|------|-----------|
|`barcode`|The observed barcode sequence.|
|`count`|The approximate number of times that barcode sequences was observed.|
