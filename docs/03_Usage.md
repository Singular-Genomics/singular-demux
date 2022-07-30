# Usage

<!---toc start-->
   * [Help](#help)
   * [Read Structures](#read-structures)
   * [Examples](#examples)
      * [Dual-Indexed Paired-End](#dual-indexed-paired-end)

<!---toc end-->

## Help

Execute the following to see the full set of command line options

```console
singular-demux --help
```

## Read Structures

A Read Structure refers to a String that describes how the bases in a sequencing run should be allocated into logical reads.
The read structure _per-read_ must be given to assign the bases in each read to logical reads (e.g. Sample Barcode, Template, Unique Molecular Identifiers, Skipped Bases).
See [the following document][read-structure-link] for more information.

## Examples

### Dual-Indexed Paired-End

```console
singular-demux \
	--sample-metadata samples.csv \
	--fastqs \
		Undetermined_S1_L001_I1_001.fastq.gz \
		Undetermined_S1_L001_R1_001.fastq.gz \
		Undetermined_S1_L001_R2_001.fastq.gz \
		Undetermined_S1_L001_I2_001.fastq.gz \
	--read-structures +B +T +T +B \
	--output-dir /path/to/output-dir;
```

[read-structure-link]: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
