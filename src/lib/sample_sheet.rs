use crate::opts::Opts;
use crate::opts::TOOL_NAME;
use crate::sample_metadata::coelesce_samples;
use crate::sample_metadata::{validate_samples, SampleMetadata};
use clap::Parser;
use csv::{ReaderBuilder, StringRecord, Trim};
use fgoxide::io::Io;
use itertools::Itertools;
use std::fmt::Display;
use thiserror::Error;

/// The optional line number from the [`SampleMetadata`] file where an error ocurred.
#[derive(Debug)]
pub struct ErrorLine(pub Option<usize>);

impl Display for ErrorLine {
    /// Writes the line number if present, nothing if it is not None.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Some(number) => write!(f, "Line {}", number),
            None => Ok(()),
        }
    }
}

/// The reason that a barcode has been deemed invalid.
#[derive(Debug)]
#[non_exhaustive]
pub enum ReasonBarcodeInvalid {
    EmptyString,
}

impl Display for ReasonBarcodeInvalid {
    /// Proper error wording for each reason a barcode is invalid.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyString => write!(f, "Barcode is an empty string"),
        }
    }
}

/// The error that may occur when parsing the [`SampleSheet`].
#[derive(Error, Debug)]
pub enum SampleSheetError {
    #[error("Io error occurred")]
    Io(#[from] std::io::Error),

    #[error("Io error occurred")]
    FgError(#[from] fgoxide::FgError),

    #[error("The sample sheet was empty")]
    Empty,

    #[error("The '[Demux]' section is missing")]
    NoDemuxHeader,

    #[error("Could not parse the demultiplexing options, {kind}: {args}")]
    DemuxOptionsParsing { kind: String, args: String },

    #[error("The '[Data]' section is missing")]
    NoData,

    #[error("The '[Data]' section is missing samples")]
    NoSamples,

    #[error("The header line in the '[Data]' section is missing")]
    NoDataHeader,

    #[error("The sample on line {line_number} had {actual} fields, expected {expected} fields, for line: {line}")]
    SampleInvalidNumberOfColumns {
        actual: usize,
        expected: usize,
        line_number: usize,
        line: String,
    },

    #[error("Unable to parse the sample info line number {line}: {source}")]
    SampleInvalidLine { source: csv::Error, line: usize },

    #[error(
        "{sample_a}:{barcode_a} and {sample_b}:{barcode_b} barcodes are within {distance} hamming distance."
    )]
    BarcodeCollision {
        sample_a: String,
        barcode_a: String,
        sample_b: String,
        barcode_b: String,
        distance: usize,
    },

    #[error("Duplicate Sample_ID found: {id}")]
    DuplicateSampleId { id: String },

    #[error("Invalid barcode sequence for {id} `{barcode}` - {reason}. {line}")]
    InvalidBarcode { barcode: String, id: String, reason: ReasonBarcodeInvalid, line: ErrorLine },

    // #[error("Io error occurred")]
    // Io(#[from] std::io::Error),
    #[error(transparent)]
    Deserialize(#[from] csv::Error),

    #[error("Unable to deserialize line number {line}")]
    DeserializeRecord { source: csv::Error, line: usize },

    #[error("Both 'Sample_Barcode' and '{index_name}' were specified on line number {line}.")]
    BarcodeColumnCollision { index_name: String, line: usize },

    #[error(
        "{sample_a}:{barcode_a} and {sample_b}:{barcode_b} have barcodes with different lengths."
    )]
    UnequalBarcodeLengths {
        sample_a: String,
        barcode_a: String,
        sample_b: String,
        barcode_b: String,
    },

    #[error("Sample metadata must include at least one sample")]
    ZeroSamples,
}

#[derive(Debug, Clone)]
pub struct SampleSheet {
    pub samples: Vec<SampleMetadata>,
    pub opts: Opts,
}

impl SampleSheet {
    /// Builds a `SampleSheet` from the CSV at the given path.
    ///
    /// The CSV may be a Singular Genomics Sample Sheet, with the optional `[Demux]` section first
    /// and required `[Data]` section next, or a simple CSV file with headers.
    ///
    /// # Sample Sheet
    ///
    /// The given default command line options are updated using the first to columns in the
    /// `[Demux]` section.  The keys must specified in the first column and be the long-form of the
    /// corresponding command line option (e.g. `fastqs` for the command line option `--fastqs`).
    /// The value(s) must be specified in the second column.  For options that take multiple
    /// values (e.g. `--read-structures`), the multiple values must be space separated.
    ///
    /// The `[Data]` section must have `Sample_Id`, `Index1_Sequence`, and `Index2_Sequence`
    /// header.  Each subsequent row corresponds to a single sample.
    ///
    /// # Arguments
    /// - `path` - the path to the sample sheet or simple CSV
    /// - `opts` - the default command line options to update using the `[Demux]` section.  The
    ///            options are also used to ensure no barcode collisions occur, and to set the
    ///            undetermined sample name when demultiplexing.
    pub fn from_path(opts: Opts) -> Result<Self, SampleSheetError> {
        // Read in all the lines so we can check if we have a simple CSV file or a full-fledged
        // Sample Sheet.
        let io = Io::default();
        let lines = io.read_lines(&opts.sample_metadata).map_err(SampleSheetError::FgError)?;

        if lines.is_empty() {
            return Err(SampleSheetError::Empty);
        }

        // If either [Demux] or [Data] found at the start of the line, assume a full-fledged
        // sample sheet.  Otherwise, fallback to a simple CSV with a header line.
        let is_sample_sheet: bool =
            lines.iter().any(|l| l.starts_with("[Demux]") || l.starts_with("[Data]"));

        // Use the read in lines again, so we don't have to re-read the input file (could be
        // piped?)
        let data = lines.join("\n");
        let reader = data.as_bytes();
        if is_sample_sheet {
            SampleSheet::from_sample_sheet_reader(reader, opts)
        } else {
            SampleSheet::from_metadata_csv_reader(reader, opts)
        }
    }

    /// Reads sample metadata and command line options from a full-fledged Singular Genomics
    /// Sample Sheet.
    ///
    /// # Arguments
    /// - `reader` - the reader from which to read the sample sheet
    /// - `opts` - the default command line options to update using the `[Demux]` section.  The
    ///            options are also used to ensure no barcode collisions occur, and to set the
    ///            undetermined sample name when demultiplexing.
    fn from_sample_sheet_reader<R: std::io::Read>(
        reader: R,
        opts: Opts,
    ) -> Result<Self, SampleSheetError> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b',')
            .has_headers(false)
            .quoting(true)
            .flexible(true)
            .trim(Trim::All)
            .from_reader(reader);

        // Read in all the records
        let mut records: Vec<StringRecord> = vec![];
        for record in reader.records() {
            let rec = record?;
            records.push(rec);
        }
        SampleSheet::from_sample_sheet_string_records(&records, opts)
    }

    /// Reads sample metadata from a simple CSV.  No command line options are able to be specified
    /// in this format; use a sample sheet in that case.
    ///
    /// # Arguments
    /// - `reader` - the reader from which to read the sample sheet
    /// - `opts` - the command line options used to ensure no barcode collisions occur, and to set
    ///            the undetermined sample name when demultiplexing.
    fn from_metadata_csv_reader<R: std::io::Read>(
        reader: R,
        opts: Opts,
    ) -> Result<Self, SampleSheetError> {
        let mut reader =
            csv::ReaderBuilder::new().has_headers(true).delimiter(b',').from_reader(reader);

        let mut samples = vec![];
        for (ordinal, record) in reader.deserialize().enumerate() {
            // Note that line numbers a +2 to account for the header and convert to 1-based counting
            let mut record: SampleMetadata = record.map_err(|e| {
                SampleSheetError::DeserializeRecord { source: e, line: ordinal + 2 }
            })?;
            record = record.update_with_and_set_demux_barcode(ordinal, ordinal + 2)?;
            samples.push(record);
        }

        let samples = coelesce_samples(samples, &opts.lane);
        let samples = validate_samples(
            samples,
            Some(opts.allowed_mismatches),
            &opts.undetermined_sample_name,
        )?;

        Ok(SampleSheet { samples, opts })
    }

    /// Finds the start and end line index (0-based inclusive) of the section with the given key
    /// (e.g.  "[Header]", "[Demux]", "[Data]"), returning `None` if the key wasn't found.  The key
    /// must occur in the first column.  The end of the section is identified by the first column
    ///  starting with "[".  The line with the header key is not returned.
    ///
    /// Arguments:
    /// - `records` - the `StringRecord`s, one per line
    /// - `section_key` - the value of the first column that identifies the start of the section, e.g. "[Data]"
    fn find_section(records: &[StringRecord], section_key: &str) -> Option<(usize, usize)> {
        let mut start_line_index = 0;
        while start_line_index < records.len() {
            let record = &records[start_line_index];
            start_line_index += 1;
            if !record.is_empty() && &record[0] == section_key {
                break;
            }
        }
        if start_line_index == records.len() {
            return None;
        }
        let mut end_line_index = start_line_index;
        while end_line_index < records.len() {
            let record = &records[end_line_index];
            if record[0].starts_with('[') {
                break;
            }
            end_line_index += 1;
        }
        Some((start_line_index, end_line_index - 1))
    }

    /// Parses the demultiplexing options from given string records and updates the given command
    /// line options.
    ///
    /// The records are one per line.
    /// Updates the given command line options with values from the `[Demux]` section.
    /// Starting at the given line index, read in any demultiplexing command line options and
    /// values.  The option names must be the same name as on the command line, with the leading
    /// `--` ommited.  If the value is empty, the command line option is assumed to be a flag.
    /// Command line options with empty values are not supported.  The given opts are updated.
    ///
    /// Returns the line index after the last demultiplex option.  The end of the demultiplexing
    /// options is identified either when no more lines remain, or when the `[Data]` section is
    /// found.
    ///
    /// # Arguments
    /// - `records` - the full list of `StringRecord`s, one per line
    /// - `opts` - the default command line options to update using the `[Demux]` section
    /// - `line_index` - the index of the first line in the demultiplexing section in the sample sheet
    fn parse_and_update_demux_options(
        records: &[StringRecord],
        mut opts: Opts,
    ) -> Result<Opts, SampleSheetError> {
        if records.is_empty() || records.iter().all(csv::StringRecord::is_empty) {
            return Ok(opts);
        }
        let mut argv: Vec<String> = vec![TOOL_NAME.to_string()];
        for record in records {
            if record.len() >= 2 {
                // Command line arguments and their values are not checked here to be valid.  This
                // happens below when we try to build the opts.
                let key = record[0].to_string();
                let argument = format!("--{}", key);
                argv.push(argument);

                // For boolean options, only the presence of a key and the value is ignored.  Thus,
                // when an empty value is given, assume the argument is a boolean.
                if !record[1].is_empty() {
                    for value in record[1].split(' ') {
                        argv.push(value.to_string());
                    }
                }
            }
        }

        // Build the opts given the arguments.  This is where any unknown command line argument
        // will raise an error, and illegal value too.
        match opts.try_update_from(argv.into_iter()) {
            Ok(()) => Ok(opts),
            Err(err) => {
                let kind = err
                    .kind
                    .as_str()
                    .unwrap_or(&format!("Bug: unknown kind: {:?}", err.kind))
                    .to_string();
                let mut args: Vec<String> = vec![];
                for string in err.info {
                    args.push(string.split(' ').next().unwrap_or(&string).to_string());
                }
                Err(SampleSheetError::DemuxOptionsParsing { kind, args: args.join(", ") })
            }
        }
    }

    /// Converts the given string records from a sample sheet to samples.  Assumes a header, and
    /// all rows have the same number of columns.  No validation is performed across samples.
    fn slurp_samples(
        records: &[StringRecord],
        mut line_index: usize,
    ) -> Result<Vec<SampleMetadata>, SampleSheetError> {
        if records.is_empty() {
            return Err(SampleSheetError::NoDataHeader);
        }

        // The header line with column names
        let header = &records[0];

        // parse the samples
        let mut ordinal = 1;
        let mut samples: Vec<SampleMetadata> = vec![];
        for record in &records[1..] {
            // allow an empty line
            if record.is_empty() {
                line_index += 1;
                continue;
            }
            // make sure we have the correct number of columns
            if header.len() != record.len() {
                return Err(SampleSheetError::SampleInvalidNumberOfColumns {
                    actual: record.len(),
                    expected: header.len(),
                    line_number: line_index + 1,
                    line: record.into_iter().join(","),
                });
            }

            // Build the sample and convert it to SampleMetadata
            let sample: SampleMetadata = record.deserialize(Some(header)).map_err(|e| {
                SampleSheetError::SampleInvalidLine { source: e, line: line_index + 1 }
            })?;
            let sample_metadata: SampleMetadata =
                sample.update_with_and_set_demux_barcode(ordinal, line_index + 1)?;
            samples.push(sample_metadata);

            ordinal += 1;
            line_index += 1;
        }
        if samples.is_empty() {
            return Err(SampleSheetError::NoSamples);
        }

        Ok(samples)
    }

    /// Reads from a Sample Sheet represented as a vector if `StringRecord`s.  Skips over lines
    /// optionally until the `[Demux]` section is found, upating the given command line options,
    /// and then the `[Data]` section, to parse a list of samples.  Validates the samples are
    /// unique and their barcodes are dissimilar enough.
    fn from_sample_sheet_string_records(
        records: &[StringRecord],
        opts: Opts,
    ) -> Result<Self, SampleSheetError> {
        if records.is_empty() {
            return Err(SampleSheetError::Empty);
        }

        // Parse demultiplexing settings ([Demux] section)
        let opts: Opts = {
            match SampleSheet::find_section(records, "[Demux]") {
                None => Ok(opts),
                Some((start, end)) => {
                    SampleSheet::parse_and_update_demux_options(&records[start..=end], opts)
                }
            }?
        };

        // Parse the sample metadata ([Data] section)
        let (start, end) = match SampleSheet::find_section(records, "[Data]") {
            None => return Err(SampleSheetError::NoData),
            Some(tuple) => tuple,
        };

        let samples = SampleSheet::slurp_samples(&records[start..=end], start)?;
        let samples = coelesce_samples(samples, &opts.lane);

        // Validate the samples
        let samples = validate_samples(
            samples,
            Some(opts.allowed_mismatches),
            &opts.undetermined_sample_name,
        )?;

        Ok(SampleSheet { samples, opts })
    }
}

#[cfg(test)]
mod test {
    use std::path::{Path, PathBuf};

    use crate::opts::Opts;
    use crate::sample_sheet::{SampleSheet, SampleSheetError};
    use bstr::BString;
    use clap::error::ErrorKind::{MissingRequiredArgument, UnknownArgument};
    use csv::StringRecord;
    use itertools::Itertools;
    use matches::assert_matches;

    #[test]
    fn test_find_section_not_found_empty_records() {
        assert_eq!(SampleSheet::find_section(&[], "[Hello]"), None);
    }

    #[test]
    fn test_find_section_not_found() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["unknown", ""]),
            StringRecord::from(vec!["[Data]"]),
        ];
        assert_eq!(SampleSheet::find_section(&records, "[Hello]"), None);
    }

    #[test]
    fn test_find_section_not_found_empty_section() {
        let records: Vec<StringRecord> = vec![StringRecord::from(vec!["[Hello]"])];
        assert_eq!(SampleSheet::find_section(&records, "[Hello]"), None);
    }

    #[test]
    fn test_find_section_non_empty_section() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["unknown", ""]),
            StringRecord::from(vec!["[Data]"]),
            StringRecord::from(vec!["unknown", ""]),
            StringRecord::from(vec!["[End]"]),
        ];
        assert_eq!(SampleSheet::find_section(&records, "[Header]"), Some((1, 2)));
        assert_eq!(SampleSheet::find_section(&records, "[Demux]"), Some((4, 6)));
        assert_eq!(SampleSheet::find_section(&records, "[Data]"), Some((8, 8)));
        assert_eq!(SampleSheet::find_section(&records, "[End]"), None);
    }

    // TODO: uncomment when https://github.com/clap-rs/clap/pull/4618 is released and the clap
    // dependency is updated
    /*
    #[test]
    fn test_demux_missing_fastqs_and_read_structures() {
        // It's ok that we do not have any read structures, as it may be a path prefix.
        let records: Vec<StringRecord> = vec![StringRecord::from(vec!["unkonwn"])];

        let result = SampleSheet::parse_and_update_demux_options(&records, Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, MissingRequiredArgument.as_str().unwrap());
            assert_eq!(args, "--fastqs".to_string());
        }
    }

    #[test]
    fn test_demux_missing_fastqs() {
        let records: Vec<StringRecord> = vec![StringRecord::from(vec!["read-structures", "8B +T"])];

        let result = SampleSheet::parse_and_update_demux_options(&records, Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, MissingRequiredArgument.as_str().unwrap());
            assert_eq!(args, "--fastqs".to_string());
        }
    }
    */

    #[test]
    fn test_demux_missing_read_structure() {
        // It's ok that we do not have any read structures, as it may be a path prefix.
        let fastqs = vec![Path::new("/dev/null")];
        let records: Vec<StringRecord> = vec![StringRecord::from(vec!["fastqs", "/dev/null"])];
        let opts = SampleSheet::parse_and_update_demux_options(&records, Opts::default()).unwrap();
        assert_eq!(opts.fastqs, fastqs);
        assert!(opts.read_structures.is_empty());
    }

    #[test]
    fn test_unknown_demux_option_flag() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["unknown", ""]),
        ];

        let result = SampleSheet::parse_and_update_demux_options(&records, Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, UnknownArgument.as_str().unwrap());
            assert_eq!(args, "--unknown".to_string());
        }
    }

    #[test]
    fn test_unknown_demux_option_with_value() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["unknown", "value"]),
        ];

        let result = SampleSheet::parse_and_update_demux_options(&records, Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, UnknownArgument.as_str().unwrap());
            assert_eq!(args, "--unknown".to_string());
        }
    }

    #[test]
    fn test_parse_and_update_demux_options_ok() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["allowed-mismatches", "123"]),
        ];

        let opts = SampleSheet::parse_and_update_demux_options(&records, Opts::default()).unwrap();
        assert_eq!(opts.fastqs, vec![PathBuf::from("/dev/null")]);
        assert_eq!(opts.read_structures.iter().map(|r| format!("{}", r)).join(" "), "8B +T");
        assert_eq!(opts.allowed_mismatches, 123);
    }

    #[test]
    fn test_demux_existing_fastqs() {
        // FASTQS is given on the command line, so already exists in Opts.  The sample sheet
        // should parse just fine without it.
        let opts = Opts { fastqs: vec![PathBuf::from("/dev/null")], ..Opts::default() };
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["allowed-mismatches", "123"]),
        ];
        let opts = SampleSheet::parse_and_update_demux_options(&records, opts).unwrap();
        assert_eq!(opts.fastqs, vec![PathBuf::from("/dev/null")]);
        assert_eq!(opts.read_structures.iter().map(|r| format!("{}", r)).join(" "), "8B +T");
        assert_eq!(opts.allowed_mismatches, 123);
    }

    #[test]
    fn test_slurp_samples_no_header() {
        let result = SampleSheet::slurp_samples(&[], 12);
        assert_matches!(result, Err(SampleSheetError::NoDataHeader));
    }

    #[test]
    fn test_slurp_samples_ok() {
        let records = vec![
            StringRecord::from(vec!["Sample_ID", "Index1_Sequence", "Index2_Sequence"]),
            StringRecord::from(vec!["S1", "AAAA", "CCCC"]),
            StringRecord::from(vec!["S2", "GGGG", "TTTT"]),
        ];
        let samples = SampleSheet::slurp_samples(&records, 12).unwrap();
        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0].sample_id, "S1");
        assert_eq!(samples[0].index1, Some(BString::from("AAAA")));
        assert_eq!(samples[0].index2, Some(BString::from("CCCC")));
        assert_eq!(samples[0].ordinal, 1);
        assert_eq!(samples[1].sample_id, "S2");
        assert_eq!(samples[1].index1, Some(BString::from("GGGG")));
        assert_eq!(samples[1].index2, Some(BString::from("TTTT")));
        assert_eq!(samples[1].ordinal, 2);
    }

    #[test]
    fn test_slurp_samples_too_few_lines() {
        let records = vec![
            StringRecord::from(vec!["Sample_ID", "Index1_Sequence", "Index2_Sequence"]),
            StringRecord::from(vec!["S1", "AAAA", "CCCC"]),
            StringRecord::from(vec!["S2", "TTTT"]),
            StringRecord::from(vec!["S3", "AGAG", "TCTC"]),
        ];
        let result = SampleSheet::slurp_samples(&records, 12);
        assert_matches!(
            result,
            Err(SampleSheetError::SampleInvalidNumberOfColumns {
                actual: _,
                expected: _,
                line_number: _,
                line: _,
            })
        );
        if let Err(SampleSheetError::SampleInvalidNumberOfColumns {
            actual,
            expected,
            line_number,
            line,
        }) = result
        {
            assert_eq!(actual, 2);
            assert_eq!(expected, 3);
            assert_eq!(line_number, 14);
            assert_eq!(line, "S2,TTTT");
        }
    }

    #[test]
    fn test_slurp_samples_too_many_lines() {
        let records = vec![
            StringRecord::from(vec!["Sample_ID", "Index1_Sequence", "Index2_Sequence"]),
            StringRecord::from(vec!["S1", "AAAA", "CCCC"]),
            StringRecord::from(vec!["S2", "GGGG", "TTTT", "ACGT"]),
            StringRecord::from(vec!["S3", "AGAG", "TCTC"]),
        ];
        let result = SampleSheet::slurp_samples(&records, 12);
        assert_matches!(
            result,
            Err(SampleSheetError::SampleInvalidNumberOfColumns {
                actual: _,
                expected: _,
                line_number: _,
                line: _,
            })
        );
        if let Err(SampleSheetError::SampleInvalidNumberOfColumns {
            actual,
            expected,
            line_number,
            line,
        }) = result
        {
            assert_eq!(actual, 4);
            assert_eq!(expected, 3);
            assert_eq!(line_number, 14);
            assert_eq!(line, "S2,GGGG,TTTT,ACGT");
        }
    }

    #[test]
    fn test_slurp_samples_deserialize_error() {
        let records = vec![
            StringRecord::from(vec!["Sample_ID", "Index1_Sequence", "Index2_Sequence", "Lane"]),
            StringRecord::from(vec!["S1", "AAAA", "CCCC", "1"]),
            StringRecord::from(vec!["S2", "GGGG", "TTTT", "ACCT"]),
        ];
        let result = SampleSheet::slurp_samples(&records, 12);
        assert_matches!(result, Err(SampleSheetError::SampleInvalidLine { source: _, line: _ }));
        if let Err(SampleSheetError::SampleInvalidLine { source: _, line }) = result {
            assert_eq!(line, 14);
        }
    }

    #[test]
    fn test_slurp_samples_both_sample_barcode_and_index_sequences() {
        let records = vec![
            StringRecord::from(vec![
                "Sample_ID",
                "Index1_Sequence",
                "Index2_Sequence",
                "Sample_Barcode",
            ]),
            StringRecord::from(vec!["S1", "AAAA", "CCCC", "AAAACCCC"]),
            StringRecord::from(vec!["S2", "GGGG", "TTTT", "GGGGTTTT"]),
            StringRecord::from(vec!["S3", "AGAG", "TCTC", "AGAGTCTC"]),
        ];
        let result = SampleSheet::slurp_samples(&records, 12);
        assert_matches!(
            result,
            Err(SampleSheetError::BarcodeColumnCollision { index_name: _, line: _ })
        );
        if let Err(SampleSheetError::BarcodeColumnCollision { index_name, line }) = result {
            assert_eq!(index_name, "Index1_Sequence and Index2_Sequence");
            assert_eq!(line, 13);
        }
    }

    // TODO:
    // - tests for from_sample_sheet_string_records
    //   - empty sample sheet
    //   - no Data section
    //   - ok sample sheet (csv, sample sheet, no demux)

    #[test]
    fn test_sample_sheet_empty() {
        let records: Vec<StringRecord> = vec![];
        assert_matches!(
            SampleSheet::from_sample_sheet_string_records(&records, Opts::default()),
            Err(SampleSheetError::Empty)
        );
    }

    #[test]
    fn test_sample_sheet_no_data() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
        ];
        assert_matches!(
            SampleSheet::from_sample_sheet_string_records(&records, Opts::default()),
            Err(SampleSheetError::NoData)
        );
    }

    #[test]
    fn test_ok_sample_sheet() {
        let file_contents ="[Header]\n\
            Date,Today\n\
            Run Name,Foo\n\
            [Demux]\n\
            Bar\n\
            fastqs,/dev/null\n\
            read-structures,8B +T\n\
            [Data]\n\
            Sample_ID,Index1_Name,Index1_Sequence,Index2_Name,Index2_SequenceLane,Lane,Lane_Name,Project,Loading_Concentration,Application,Notes,Reference\n\
            S1,I1,AAAA,I2,CCCC,1,Lane1,S1_Project,,,,\n\
            S2,I1,GGGG,I2,TTTT,1,Lane1,S2_Project,,,,\n\
            S3,I1,ATTA,I2,GCCG,1,Lane1,S3_Project,,,,\n\
            ".to_string();

        let dir = tempfile::tempdir().unwrap();
        let output = dir.as_ref().join("sample_metadata.csv");
        std::fs::write(&output, file_contents).expect("Failed to write sample metadata to file.");
        let opts = Opts { sample_metadata: output, ..Opts::default() };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        // read structure
        assert_eq!(
            sample_sheet.opts.read_structures.iter().map(|r| format!("{}", r)).join(""),
            "8B+T"
        );

        // input path
        assert_eq!(
            sample_sheet.opts.fastqs.iter().map(|r| format!("{:?}", r)).collect::<Vec<String>>(),
            vec!["\"/dev/null\""]
        );

        // samples (one extra for the undetermined)
        assert_eq!(sample_sheet.samples.len(), 4);
    }

    #[test]
    fn test_ok_metadata_csv() {
        let file_contents = "Sample_ID,Index1_Name,Index1_Sequence,Index2_Name,Index2_SequenceLane,Lane,Lane_Name,Project,Loading_Concentration,Application,Notes,Reference\n\
            S1,I1,AAAA,I2,CCCC,1,Lane1,S1_Project,,,,\n\
            S2,I1,GGGG,I2,TTTT,1,Lane1,S2_Project,,,,\n\
            S3,I1,ATTA,I2,GCCG,1,Lane1,S3_Project,,,,\n\
            ".to_string();
        let dir = tempfile::tempdir().unwrap();
        let output = dir.as_ref().join("sample_metadata.csv");
        std::fs::write(&output, file_contents).expect("Failed to write sample metadata to file.");
        let opts = Opts { sample_metadata: output, ..Opts::default() };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        // samples (one extra for the undetermined)
        assert_eq!(sample_sheet.samples.len(), 4);
    }

    #[test]
    fn test_ok_sample_sheet_no_demux() {
        let file_contents ="[Header]\n\
            Date,Today\n\
            Run Name,Foo\n\
            [Demux]\n\
            Bar\n\
            fastqs,/dev/null\n\
            read-structures,+T\n\
            [Data]\n\
            Sample_ID,Index1_Name,Index1_Sequence,Index2_Name,Index2_SequenceLane,Lane,Lane_Name,Project,Loading_Concentration,Application,Notes,Reference\n\
            S1,I1,,I2,,1,Lane1,S1_Project,,,,\n\
            ".to_string();

        let dir = tempfile::tempdir().unwrap();
        let output = dir.as_ref().join("sample_metadata.csv");
        std::fs::write(&output, file_contents).expect("Failed to write sample metadata to file.");
        let opts = Opts { sample_metadata: output, ..Opts::default() };

        let sample_sheet = SampleSheet::from_path(opts).unwrap();

        // read structure
        assert_eq!(
            sample_sheet.opts.read_structures.iter().map(|r| format!("{}", r)).join(""),
            "+T"
        );

        // input path
        assert_eq!(
            sample_sheet.opts.fastqs.iter().map(|r| format!("{:?}", r)).collect::<Vec<String>>(),
            vec!["\"/dev/null\""]
        );

        // samples (one extra for the undetermined)
        assert_eq!(sample_sheet.samples.len(), 1);
        assert_eq!(sample_sheet.samples[0].barcode.to_string(), "");
    }

    // TODO test same sample, multiple lanes
}
