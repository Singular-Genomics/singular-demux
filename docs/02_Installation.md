# Installation

<!---toc start-->
   * [From Bioconda](#from-bioconda)
   * [From Releases](#from-releases)
   * [From Source](#from-source)

<!---toc end-->

## From Bioconda

Install from [`bioconda`][bioconda-link] with:

```console
conda install -c bioconda singular-demux
```

## From Releases

Install from pre-built binaries on the [Releases page][releases-link]

## From Source

1. Install [rust and cargo][cargo-install-link]


2. Build the tool

```console
$ cargo build --release
```

Then run:

```console
$ ./target/release/singular-demux --help
```

[bioconda-link]:      https://bioconda.github.io/
[releases-link]:      https://github.com/Singular-Genomics/singular-demux/releases
[cargo-install-link]: https://doc.rust-lang.org/cargo/getting-started/installation.html
