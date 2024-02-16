# sgdemux

[![Check](https://github.com/Singular-Genomics/singular-demux/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/Singular-Genomics/singular-demux/actions/workflows/build_and_test.yml)

This repository is home to the `sgdemux` tool for demultiplexing sequencing data generated on Singular Genomics' sequencing instruments.


* [Installation](#installation)
  * [From Bioconda](#from-bioconda)
  * [From Releases](#from-releases)
  * [From Source](#from-source)
* [Contributing](#contributing)
* [Overview](#overview)
* [Usage](#usage)
  * [Inputs](#inputs)
      * [FASTQ Files](#fastq-files)
        * [Auto-detecting FASTQS from a Path Prefix](#auto-detecting-fastqs-from-a-path-prefix)
      * [Read Structures](#read-structures)
      * [Specifying Sample Information](#specifying-sample-information)
        * [Sample Sheet](#sample-sheet)
        * [Simple Two-column CSV](#simple-two-column-csv)
      * [Full Argument List](#full-argument-list)
      * [Performance Considerations](#performance-considerations)
  * [Outputs](#outputs)
      * [Demultiplexed FASTQs](#demultiplexed-fastqs)
      * [Metrics](#metrics)
        * [per_sample_metrics.tsv](#per_sample_metricstsv)
        * [per_project_metrics.tsv](#per_project_metricstsv)
        * [metrics.tsv](#metricstsv)
        * [most_frequent_unmatched.tsv](#most_frequent_unmatchedtsv)
        * [sample_barcode_hop_metrics.tsv](#sample_barcode_hop_metricstsv)
* [Advance Usage](#advance-usage)
  * [Single Sample](#single-sample)

## Installation

`sgdemux` may be installed from bioconda, downloaded from the releases page, or built from source.

### From Bioconda

Install from [`bioconda`](https://bioconda.github.io/) with:

```console
conda create -n sgdemux -c bioconda sgdemux
conda activate sgdemux
```

### From Releases

Install from pre-built binaries on the [Releases page][releases]

[releases]: https://github.com/Singular-Genomics/singular-demux/releases

### From Source

1. Install [rust and cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html)
2. Install dependencies. For example, cmake and build-essentials are required for Ubuntu 22.04. Install using commands below.

```console
sudo apt-get update
sudo apt-get install build-essential cmake -y
```
Note: cmake for older OS version such as Ubuntu 18.04 is not incompatible.

3. Clone the repo and build the software:

```console
git clone https://github.com/Singular-Genomics/singular-demux.git
cd singular-demux
cargo install --path ../singular-demux --locked
```

## Contributing

Contributions are welcome.  See the [Contributing Guidelines](Contributing.md) for details.

## Overview

`sgdemux` performs sample demultiplexing on block-compressed (BGZF) FASTQs such as those produced by the Singular Genomics G4 platform.
The input FASTQs _must_ be block compressed (e.g. with [`bgzip`](http://www.htslib.org/doc/bgzip.html)); uncompressed or non-bgzf gzipped input files are not supported as performance would be significantly degraded.

The primary options that affect demultiplexing are `--allowed-mismatches` and `--min-delta`.  Together these specify a) how well a sample barcode in a sequencing read must match an expected barcode and b) how much worse the next best match must be.
The default options of `--allowed-mismatches 1 --min-delta 2` will only match a set of FASTQ records to an expected barcode if, across all barcode reads, there is at most one mismatch (allowed mismatches) vs. the expected barcode _and_ the difference (minimmum delta) between the number of mismatches of the best and second best matching barcode is greater than two mismatches.
Note: the allowed mismatches is not used when determining the next-best matching barcode.

For additional examples, consider `--allowed-mismatches 3 --min-delta 1`, with two barcodes `b1` and `b2`:

1. If b1 matches with 2 mismatches, and b2 matches with 3 mismatches, then the delta between the number of mismatches is 1, which _is not_ greater than `--min-delta`, and therefore the read is not assigned to a barcode.
2. If b1 matches with 1 mismatch, and b2 matches with 3 mismatches, then the delta between the number of mismatches is 2, which _is_ greater than `--min-delta`, and therefore the read is assigned to barcode `b1`.
3. If b1 matches with 0 mismatch, and b2 matches with 2 mismatches, then the delta between the number of mismatches is 2, which _is_ greater than `--min-delta`, and therefore the read is assigned to barcode `b1`.
4. If b1 matches with 0 mismatches, and b2 matches with 1 mismatches, then the delta between the number of mismatches is 1, which _is not_ greater than `--min-delta`, and therefore the read is not assigned to a barcode.
5. If b1 matches with 3 mismatches, and b2 matches with 6 mismatches, then the delta between the number of mismatches is 2, which _is_ greater than `--min-delta`, and the number of mismatches for b1 is less than equal to than `--allowed-mismatches`, and thefore read is assigned to barcode `b1`.
6. If b1 matches with 4 mismatches, and b2 matches with 6 mismatches, then the delta between the number of mismatches is 2, which _is_ greater than `--min-delta`, and but the number of mismatches for b1 is greater than `--allowed-mismatches`, and thefore the read is not assigned to a barcode.
7. If b1 matches with 2 mismatch, and b2 matches with 2 mismatches, then the delta between the number of mismatches is 0, which _is not_ greater than `--min-delta`, and therefore the read is not assigned to a
barcode.

Several other options affect how demultiplexing is performed, and for these to be fully understood it is necessary to understand the order in which they are applied in the demultiplexing process.  Operations are ordered as follows:

1. A record is read in from each of the input FASTQ files and broken into read "segments" using the supplied read structures.
2. If `--filter-control-reads` is specified and the reads are marked as controls in the FASTQ header, the reads are discarded (i.e. they do not get written to _any_ output files).
3. If `--filter-failing-quality` is specified and the reads are marked as quality failures in the FASTQ header, the reads are discarded (i.e. they do not get written to _any_ output files).
4. If one or more `--quality-mask-threshold` values are supplied, template bases in all input reads that have base quality below the given threshold value are masked to `N`.
5. Match the reads against the set of expected barcodes; if the sample barcode has more `N` bases in it that specified by `--max-no-calls` or does not match to an expected barcode within defined parameters, the reads will be assigned to the undetermined sample.
6. Write out the subset of the FASTQs/read segments specified by `--output-types` to the FASTQ file(s) for the assigned sample.

## Usage

The primary inputs to the tool are:

1. A set of undemultiplexed FASTQ files (BGZF compressed)
2. A set of read-structures, one per input FASTQ file
3. A file of sample metadata including sample names and barcode sequences
4. A directory into which the demultiplexed FASTQ files should be written

Reads are written to per-sample, per-instrument-read files within the output directory.  An additional `Undetermined` set of files will be written containing those reads that did not match any expected barcodes.

An example invocation follows:

```shell
sgdemux \
  --fastqs R1.fastq.gz R2.fastq.gz I1.fastq.gz I2.fastq.gz \
  --read-structures +T +T 8B 8B \
  --sample-metadata sample-metadata.csv \
  --output-dir demuxed/
```

### Inputs

#### FASTQ Files

The full set of FASTQ files generated for a run, or lane, or sequencing should be provided, including all template and index reads.  For example if a paired-end sequencing run was performed with dual sample index reads, four files should be provided:

```shell
  --fastqs R1.fastq.gz R2.fastq.gz I1.fastq.gz I2.fastq.gz
```

If multiple FASTQ files are available per instrument reads, they should be concatenated prior to running `sgdemux`.
BGZF files, due to their block-compressed nature, can be concatenated simply using standard `cat`, e.g.:

```shell
for read in R1 R2 I1 I2; do cat L*/${read}.fastq.gz > ./${read}.fastq.gz; done
```

FASTQ files _must_ be BGZF compressed.

##### Auto-detecting FASTQS from a Path Prefix

Alternatively, the FASTQS can be auto-detected when a path prefix is given to `--fastqs <dir>/<prefix>`.
The FASTQs must be named `<dir>/<prefix>_L00<lane>_<kind><kind-number>_001.fastq.gz`, where `kind` is
one of R (read/template), I (index/sample barcode), or U (umi/molecular barcode).

The Read Structure will be derived from file names (kind and kind number), with the full read length used for the given kind.
The derived Read Structure and FASTQs will be ordered first by `kind` (I then R then U), second by
read number (e.g. R1 before R2).  This is important for command line options that can be specified once per read kind and number.
 E.g. if the following FASTQs are present with path prefix `/path/to/prefix`:

```
/path/to/prefix_L001_I1_001.fastq.gz
/path/to/prefix_L001_I2_001.fastq.gz
/path/to/prefix_L001_R1_001.fastq.gz
/path/to/prefix_L001_R2_001.fastq.gz
```

then the `+B +B +T +T` read structure will be used.  Since this tool requires all sample barcode
segments to have a fixed length, the first read in any index/sample-barcode FASTQ will be examined
and its length used as the expected sample barcode length.

Furthermore, multiple lanes may be given and will be used for demultiplexing:

```
/path/to/prefix_L001_I1_001.fastq.gz
/path/to/prefix_L002_I1_001.fastq.gz
/path/to/prefix_L001_I2_001.fastq.gz
/path/to/prefix_L002_I2_001.fastq.gz
/path/to/prefix_L001_R1_001.fastq.gz
/path/to/prefix_L002_R1_001.fastq.gz
/path/to/prefix_L001_R2_001.fastq.gz
/path/to/prefix_L002_R2_001.fastq.gz
```

When data for multiple lanes is provided, each lane must have the same number and types of input fastqs.

The auto-detected/derived Read Structure may be overriden on the command line or in the sample sheet
by providing the `--read-structures` argument.  In this case, the new read structure must be given
and will be applied in the same order as described above (e.g. I1, I2, R1, R2 for a dual index paired end run).

#### Read Structures

Read Structures are short strings that describe the origin and/or purpose of bases within sequencing reads.  They are made up of a sequence of `<number><operator>` pairs (segments).  Four kinds of operators are recognized:

1. **T** identifies template reads/bases
2. **B** identifies sample barcode reads/bases
3. **M** identifies unique molecular index reads/bases
4. **S** identifies a set of bases to be skipped or ignored

The last `<number><operator>` pair in a Read Structure may use `+` instead of a number to denote "all remaining bases".
This is useful if, e.g., FASTQs have been trimmed and/or contain reads of varying length.

For more details on Read Structures, and how to validate them, see [this detailed description](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures).

Read Structures are not required to be provided when using a path prefix for the input FASTQs.  In that case,
the read structure will be inferred from the FASTQ name.  See: [Auto-detecting FASTQS from a Path Prefix](#auto-detecting-fastqs-from-a-path-prefix).

When providing the input FASTQs explicitly, one Read Structure must be provided for each input FASTQ file, in the same order.  Matching the set of reads specified in the FASTQ files section above one might specify:

```shell
  --read-structures +T +T 8B 8B
```

All sample barcode segments must be a fixed length.  E.g. `8B+T` is allowed but `10S+B` is not.

#### Specifying Sample Information

The sample metadata file may be a Sample Sheet or a simple two-column CSV file with headers.

##### Sample Sheet

Information about the sample(s) to demultiplex is specified within a Sample Sheet.
Command line options for demultiplexing may also be passed via the Sample Sheet.

The Sample Sheet may have a `[Demux]` section for command line options, and must have a `[Data]`
section for sample information.

The `[Demux]` section must contain a line per command line option.
The first column must contain the option long name _without_ the leading `--` (e.g. `fastqs` or
`read-structures`).
The second column contains the option value, or empty if the option takes no value (i.e. a flag).
If the option accepts multiple values, they must be space separated in the second column.
The command line options specified in the sample sheet override those provided on the command line.
The order of the FASTQs must match the order read structures.

The `[Data]` section must contain a header line.
The `Sample_ID` column must contain a unique, non-empty identifier
for each sample.  One or both of `Index1_Sequence` and `Index2_Sequence` must be present with values for
indexed runs.  For non-indexed runs, a single sample must be given with an empty value for both the
`Index1_Sequence` and `Index2_Sequence` columns.
Both `Sample_ID`s and the `Index1_Sequence`/`Index2_Sequence` combinations must be unique within the file, and both columns are required for all samples.

An example follows:

```text
[Demux]
fastqs,/path/to/i1.fq.gz /path/to/r1.fq.gz
read-structures,+B +T
[Data]
Sample_ID,Index1_Sequence,Index2_Sequence
s1,ACTGGTCA,
s2,ATACGAAC,
```

##### Simple Two-column CSV

For the simple two-column CSV, the `Sample_Barcode` column must contain the unique set of sample barcode bases for the sample(s).
If multiple sample barcodes are are present (e.g. dual-indexing runs, additional inline sample indices) then the `Sample_Barcode` field should contain the full set of barcode bases expected to be read for the sample.
The ordering of the concatenated barcodes is important, and should match the ordering of the FASTQs and Read Structures given.
Both `Sample_ID`s and `Sample_Barcode`s must be unique within the file, and both columns are required for all samples.
An example follows:

```text
Sample_ID,Sample_Barcode
s1,ACTGGTCA
s2,ATACGAAC
```

For example if a dual-indexing run was performed with an additional inline sample barcode in read 1, and `sgdemux` was invoked with the following options:

```console
--fastqs R1.fastq.gz I1.fastq.gz I2.fastq.gz R2.fastq.gz \
--read-structures 10B+T 8B 8B +T
```

then the `Sample_Barcode` field for each sample should be composed as follows:

```shell
  {10 base inline index}-{8 base I1 index}-{8 base I2 index}
```

#### Full Argument List

| Argument Name                    | Required | Default Value | Description                                                                                                                                                                                                                                                                                                                     |
|----------------------------------|----------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --fastqs                         | Yes      | n/a           | Path(s) to the input FASTQs, or path prefix if not a file.                                                                                                                                                                                                                                                                      |
| --sample-metadata                | Yes      | n/a           | Path to CSV of sample metadata with sample IDs and barcode sequences.                                                                                                                                                                                                                                                           |
| --read-structures                | No       | n/a           | Read structures, one per input FASTQ. Do not provide when using a path prefix for FASTQs.                                                                                                                                                                                                                                       |
| --output-dir                     | Yes      | n/a           | Path to an output directory to write into.                                                                                                                                                                                                                                                                                      |
| --allowed-mismatches             | No       | 1             | The number of mismatches allowed, in total, between expected and observed barcode bases in order to match a read to a sample.                                                                                                                                                                                                   |
| --min-delta                      | No       | 2             | The minimum number of mismatches by which the best match for a read is better than the next-best match for a read in order to accept the best match.                                                                                                                                                                            |
| --free-ns                        | No       | 1             | The number of observed Ns (no-calls) in the barcode read(s) that are allowed for "free" before treating subsequent Ns as mismatches.                                                                                                                                                                                            |
| --max-no-calls                   | No       | n/a           | If specified, do not match any reads whose barcode reads contain more than this many Ns.                                                                                                                                                                                                                                        |
| --quality-mask-threshold         | No       | n/a           | Mask to N template bases in all input reads whose base quality is below the specified value(s). A single value may be specified, which is then applied to all input reads/FASTQs. Alternatively one value per input FASTQ may be provided in the same order as the FASTQs. Sample barcode/index and UMI bases are never masked. |
| --filter-control-reads           | No       | False         | If true, filter out reads marked as control reads in their FASTQ headers.                                                                                                                                                                                                                                                       |
| --filter-failing-quality         | No       | False         | If true, filter out reads marked as failing quality control in their FASTQ headers.                                                                                                                                                                                                                                             |
| --output-types                   | No       | T             | The types of bases/reads for which output files should be generated. A single string containing one or more of `T` (template), `B` (sample barcode), `M` (UMI), and `S` (skipped).                                                                                                                                              |
| --undetermined-sample-name       | No       | Undetermined  | The name used as a prefix to generate FASTQ files for reads that didn't match to any sample.                                                                                                                                                                                                                                    |
| --most-unmatched-to-output       | No       | 1000          | Report on the top N most frequently observed unmatched barcode sequences.                                                                                                                                                                                                                                                       |
| --demux-threads                  | No       | 4             | The number of threads to use to perform demultiplexing in memory.                                                                                                                                                                                                                                                               |
| --compressor-threads             | No       | 12            | The number of threads to use in compressing the output FASTQ files.                                                                                                                                                                                                                                                             |
| --writer-threads                 | No       | 5             | The number of threads to use to write compressed FASTQ data to disk.                                                                                                                                                                                                                                                            |
| --override-matcher               | No       | n/a           | The algorithm to use for matching, either `CachedHammingDistance` or `PreCompute`.  By default if barcodes are 12bp or shorter `PreCompute` is used which pre-computes all possible matches, or if barcodes are longer than 12bp `CachedHammingDistance` is used which calculates matches when needed then caches the results.  |
| --skip-read-name-check           | No       | False         | If this is true, then all the read names across FASTQs will not be enforced to be the same. This may be useful when the read names are known to be the same and performance matters. Regardless, the first read name in each FASTQ will always be checked.                                                                      |
| --sample-barcode-in-fastq-header | No       | False         | If this is true, then the sample barcode is expected to be in the FASTQ read header.  For dual indexed data, the barcodes must be `+` (plus) delimited.  Additionally, if true, then neither index FASTQ files nor sample barcode segments in the read structure may be specified.                                              |
| --metric-prefix                  | No       | n/a           | Prepend this prefix to all output metric file names.                                                                                                                                                                                                                                                                            |
| --lane                           | No       | n/a           | Select a subset of lanes to demultiplex.  Will cause only samples and input FASTQs with the given `Lane`(s) to be demultiplexed.  Samples without a lane will be ignored, and FASTQs without lane information will be ignored.                                                                                                  |

#### Performance Considerations

Various `--*-threads` options are available to control the number of threads used by `sgdemux` for various purposes.  The defaults are intended to fully utilize a 32-core machine.  The defaults to the available options do not add up to 32 as several threads are used to read the input FASTQ files and for ancillary purposes.

For running on larger or smaller instances it is advised to start with the following and tune from there:

- 1/3 of available threads for compression
- 1/6 of available threads for writing
- 1/6-1/3 of available threads for demultiplexing

Currently this tool does not provide a way place a hard limit on the number of threads used.

### Outputs

#### Demultiplexed FASTQs

One or more BGZF compressed FASTQ files will be created per sample in the specified output directory. For
paired end data, the output will have the suffix `_R1.fastq.gz` and `_R2.fastq.gz` for read one and read two
respectively.

Samples barcodes, and unique molecular indices (UMIs), will be inserted into the FASTQ headers if present.  If either multiple sample barcodes or multiple UMIs are present they will be concatenated with `+` between individual barcodes prior to insertion.  For example if a FASTQ record had sample barcodes `ACGT` and `TTGA`, and UMIs of `ACCTAG` and `TCCTGG` the the output header might look like:

```shell
sg001:17:A30ZZ:1:4:1234:4567:ACCTAG+TCCTGG 1:N:1:ACGT+TTGA
```

#### Metrics

Up to five metrics files are generated to help assess run and demultiplexing quality:

##### `per_sample_metrics.tsv`

This file always produced and contains the following columns:

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

The `per_sample_metrics.tsv` file produces a row per sample.

##### `per_project_metrics.tsv`

The `per_project_metrics.tsv` file aggregates the metrics by project (aggregates the metrics across
samples with the same project) and has the same columns as [per_sample_metrics.tsv](#per_sample_metricstsv).
In this case, `barcode_name` and `library_name` will contain the project name (or `None` if no
project is given).
THe `barcode` will contain all `N`s.
The undetermined sample will not be aggregated with any other sample.

##### `metrics.tsv`

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


## Advance Usage

### Single Sample

It is possible to run `sgdemux` on a single sample _without_ demultiplexing, in order to make use of the remaining functionality such as filtering control reads, extracting UMIs, etc.  This mode is enabled by providing a sample metadata file that contains a single sample, with no barcode sequence.  For example:

```text
Sample_ID,Sample_Barcode
lone_sample,
```

The `Sample_Barcode` column must still be present, but empty for the sample.  When running in this mode:

* All reads are assigned to the single sample
* No `Undetermined` files are created
* Sample barcodes, if read, will be inserted into the headers of the output FASTQ reads
