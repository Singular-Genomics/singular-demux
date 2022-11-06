use crate::opts::Opts;
use crate::opts::TOOL_NAME;
use crate::sample_metadata::{validate_samples, SampleMetadata};
use clap::Parser;
use csv::{ReaderBuilder, StringRecord, Trim};
use fgoxide::io::Io;
use itertools::Itertools;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::fmt::Display;
use std::path::Path;
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

    #[error("Both 'Sample_Barcode' and '{index_name}' were specified on line nmber {line}.")]
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
    /// The CSV may be a Singular Genomics Sample Sheet, with the optional `[Demux]` section and
    /// optiona `[Data]` section, or a simple CSV file with
    pub fn from_path<P: AsRef<Path>>(path: &P, opts: &Opts) -> Result<Self, SampleSheetError> {
        // Build the input CSV reader
        let io = Io::default();
        let lines = io.read_lines(path).unwrap();

        if lines.is_empty() {
            return Err(SampleSheetError::Empty);
        }

        // If either [Demux] or [Data] found at the start of the line, assume a full-fledged
        // sample sheet.  Otherwise, fallback to a simple CSV with a header line.
        let is_sample_sheet: bool = {
            let mut found = false;
            for line in &lines {
                if line.starts_with("[Demux]") || line.starts_with("[Data]") {
                    found = true;
                    break;
                }
            }
            found
        };

        let data = lines.join("\n");
        let reader = data.as_bytes();
        if is_sample_sheet {
            SampleSheet::from_sample_sheet_reader(reader, opts)
        } else {
            SampleSheet::from_metadata_csv_reader(reader, opts)
        }
    }

    fn from_sample_sheet_reader<R: std::io::Read>(
        reader: R,
        opts: &Opts,
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
            let rec = record.unwrap();
            records.push(rec);
        }
        SampleSheet::from_string_records(&records, opts)
    }

    fn from_metadata_csv_reader<R: std::io::Read>(
        reader: R,
        opts: &Opts,
    ) -> Result<Self, SampleSheetError> {
        let mut reader =
            csv::ReaderBuilder::new().has_headers(true).delimiter(b',').from_reader(reader);

        let mut samples = vec![];
        for (ordinal, record) in reader.deserialize().enumerate() {
            // Note that line numbers a +2 to account for the header and convert to 1-based counting
            let mut record: SampleMetadata = record.map_err(|e| {
                SampleSheetError::DeserializeRecord { source: e, line: ordinal + 2 }
            })?;
            record = record.update_with(ordinal, ordinal + 2)?;
            println!("record {:?}", record);
            samples.push(record);
        }

        let samples = validate_samples(
            samples,
            Some(opts.allowed_mismatches),
            &opts.undetermined_sample_name,
        )?;

        Ok(SampleSheet { samples, opts: opts.clone() })
    }

    fn slurp_demux_opts(
        records: &Vec<StringRecord>,
        start_line_index: usize,
        opts: &mut Opts,
    ) -> Result<usize, SampleSheetError> {
        let mut line_index = start_line_index;
        if records.len() == line_index {
            // no "[Demux]" found
            return Err(SampleSheetError::NoDemuxHeader);
        }
        line_index += 1; // skip over "[Demux]"
        let mut argv: Vec<String> = vec![TOOL_NAME.to_string()];
        while line_index < records.len() {
            let record = &records[line_index];
            if !record.is_empty() && &record[0] == "[Data]" {
                break;
            }
            if record.len() >= 2 {
                // Command line arguments and their values are not checked here to be valid.  This
                // happens below when we try to build the opts.
                let key = record[0].to_string();
                let argument = format!("--{}", key);
                argv.push(argument);

                // For boolean options, only the presence of a key and the value is ignored.  Thus,
                // when an empty value is given, assume the argument is a boolean.
                if !record[1].is_empty() {
                    argv.push(record[1].to_string());
                }
            }
            line_index += 1;
        }

        // Build the opts given the arguments.  This is where any unknown command line argument
        // will raise an error, and illegal value too.
        match opts.try_update_from(argv.into_iter()) {
            Ok(()) => Ok(line_index),
            Err(err) => {
                let kind = err.kind.as_str().unwrap_or("Unknown error").to_string();
                let mut args: Vec<String> = vec![];
                for string in err.info {
                    args.push(string.split(' ').next().unwrap_or(&string).to_string());
                }
                Err(SampleSheetError::DemuxOptionsParsing { kind, args: args.join(", ") })
            }
        }
    }

    fn from_string_records(
        records: &Vec<StringRecord>,
        opts: &Opts,
    ) -> Result<Self, SampleSheetError> {
        if records.is_empty() {
            return Err(SampleSheetError::Empty);
        }

        let mut line_index = 0;

        // ignore data pre-"[Demux]"
        while line_index < records.len()
            && !records[line_index].is_empty()
            && &records[line_index][0] != "[Demux]"
        {
            line_index += 1;
        }

        // Parse demultiplexing settings
        let mut opts = opts.clone();
        let mut line_index = SampleSheet::slurp_demux_opts(records, line_index, &mut opts)?;
        if records.len() == line_index {
            // no "[Data]" found
            return Err(SampleSheetError::NoData);
        }

        // Parse data
        line_index += 1; // skip over the "[Data]" line
        if records.len() == line_index {
            // no header line in the "[Data]" section
            return Err(SampleSheetError::NoSamples);
        }
        let data_header = &records[line_index];
        //let data_header: Vec<&str> = records[line_index].into_iter().collect();
        line_index += 1; // skip over the header
        let mut ordinal = 1;
        let mut samples: Vec<SampleMetadata> = vec![];
        while line_index < records.len() {
            // allow an empty line
            // TODO: allow a line with the correct # of columns but all empty lines
            if records[line_index].is_empty() {
                line_index += 1;
                continue;
            }
            // make sure we have the correct number of columns
            if data_header.len() != records[line_index].len() {
                return Err(SampleSheetError::SampleInvalidNumberOfColumns {
                    actual: records[line_index].len(),
                    expected: data_header.len(),
                    line_number: line_index + 1,
                    line: records[line_index].into_iter().join(","),
                });
            }

            // Build the sample and convert it to SampleMetadata
            let sample: SampleMetadata =
                records[line_index].deserialize(Some(data_header)).map_err(|e| {
                    SampleSheetError::SampleInvalidLine { source: e, line: line_index + 1 }
                })?;
            let sample_metadata: SampleMetadata = sample.update_with(ordinal, line_index + 1)?;
            samples.push(sample_metadata);

            ordinal += 1;
            line_index += 1;
        }

        let samples = validate_samples(
            samples,
            Some(opts.allowed_mismatches),
            &opts.undetermined_sample_name,
        )?;

        Ok(SampleSheet { opts: opts.clone(), samples })
    }
}

#[cfg(test)]
mod test {
    use crate::opts::Opts;
    use crate::sample_sheet::{SampleSheet, SampleSheetError};
    use clap::error::ErrorKind::{MissingRequiredArgument, UnknownArgument};
    use csv::StringRecord;
    use itertools::Itertools;
    use matches::assert_matches;

    #[test]
    fn test_empty_sample_sheet() {
        let records: Vec<StringRecord> = vec![];
        assert_matches!(
            SampleSheet::from_string_records(&records, &Opts::default()),
            Err(SampleSheetError::Empty)
        );
    }

    #[test]
    fn test_unknown_demux_option_flag() {
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

        let result = SampleSheet::from_string_records(&records, &Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, UnknownArgument.as_str().unwrap());
            assert_eq!(args, "--unknown".to_string());
        }
    }

    #[test]
    fn test_unknown_demux_option_with_value() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["unknown", "value"]),
            StringRecord::from(vec!["[Data]"]),
        ];

        let result = SampleSheet::from_string_records(&records, &Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, UnknownArgument.as_str().unwrap());
            assert_eq!(args, "--unknown".to_string());
        }
    }

    #[test]
    fn test_demux_missing_fastqs_and_read_structures() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["Bar"]),
        ];

        let result = SampleSheet::from_string_records(&records, &Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, MissingRequiredArgument.as_str().unwrap());
            assert_eq!(args, "--fastqs, --read-structures".to_string());
        }
    }

    #[test]
    fn test_demux_missing_fastqs() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
        ];

        let result = SampleSheet::from_string_records(&records, &Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, MissingRequiredArgument.as_str().unwrap());
            assert_eq!(args, "--fastqs".to_string());
        }
    }

    #[test]
    fn test_demux_missing_read_structure() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
        ];

        let result = SampleSheet::from_string_records(&records, &Opts::default());
        assert_matches!(result, Err(SampleSheetError::DemuxOptionsParsing { kind: _, args: _ }));
        if let Err(SampleSheetError::DemuxOptionsParsing { kind, args }) = result {
            assert_eq!(kind, MissingRequiredArgument.as_str().unwrap());
            assert_eq!(args, "--read-structures".to_string());
        }
    }

    #[test]
    fn test_no_data() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
        ];
        assert_matches!(
            SampleSheet::from_string_records(&records, &Opts::default()),
            Err(SampleSheetError::NoData)
        );
    }

    #[test]
    fn test_no_samples() {
        let records: Vec<StringRecord> = vec![
            StringRecord::from(vec!["[Header]"]),
            StringRecord::from(vec!["Date", "Today"]),
            StringRecord::from(vec!["Run Name", "Foo"]),
            StringRecord::from(vec!["[Demux]"]),
            StringRecord::from(vec!["read-structures", "8B +T"]),
            StringRecord::from(vec!["fastqs", "/dev/null"]),
            StringRecord::from(vec!["[Data]"]),
        ];
        assert_matches!(
            SampleSheet::from_string_records(&records, &Opts::default()),
            Err(SampleSheetError::NoSamples)
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

        let sample_sheet = SampleSheet::from_path(&output, &Opts::default()).unwrap();

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

        let sample_sheet = SampleSheet::from_path(&output, &Opts::default()).unwrap();

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

        let sample_sheet = SampleSheet::from_path(&output, &Opts::default()).unwrap();

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
