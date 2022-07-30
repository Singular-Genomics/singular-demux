//! A library of functionality for demultiplexing FASTQ files.
//!
//! # Overview
//!
//! The flow of data is as follows:
//!
//! - The [`thread_reader::ThreadReader`] extracts chunks of FASTQs from input files, zipping the chunks
//!   together for processing.
//! - The [`demux::Demultiplexer`] struct performs the demultiplexing, converting the [`demux::PerFastqRecordSet`] into
//!   the [`demux::DemuxedGroup`].
//! - The [`pooled_sample_writer::PooledSampleWriter`] takes each of the [`demux::OutputPerSampleReads`] in the [`demux::DemuxedGroup`]
//!   and writes them to their respective writers.
//! - [`metrics`] are collected on the writer threads, collated when those threads are joined.
#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
pub mod demux;
pub mod fastq_header;
pub mod matcher;
pub mod metrics;
pub mod pooled_sample_writer;
pub mod run;
pub mod sample_metadata;
pub mod thread_reader;
pub mod utils;
