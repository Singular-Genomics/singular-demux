#![forbid(unsafe_code)]
#![allow(clippy::must_use_candidate)]
//! [`SampleMetadata`] is a small utility for parsing simplified sample sheets of the following format:
//!
//! * **Column Header**: Sample_Barcode
//!     * **Description**: The sample barcode sequence. <br>
//!                    If the sample barcode is present across multiple reads (ex. dual-index, or inline in both reads of a pair),
//!                    then the expected barcode bases from each read should be concatenated in the same order as the order of the reads' FASTQs and read structures given to this tool.
//!                    Non-ACGT characters will be removed.
//!     * **Required**: True
//! * **Column Header**: Sample_ID
//!     * **Description**: The unique identifier for the sample
//!     * **Required**: True

use std::{collections::HashSet, fmt::Display, path::Path};

use bstr::{BStr, BString};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// The bases that are allowed in the [`SampleMetadata::barcode`].
const ALLOWED_BASES: &[u8] = &[b'A', b'C', b'T', b'G'];

/// The optional line number from the [`SampleMetadata`] file where an error ocurred.
#[derive(Debug)]
pub struct ErrorLine(Option<usize>);

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

/// The error that may occur when parsing the [`SampleMetadata`].
#[derive(Error, Debug)]
pub enum SampleMetadataError {
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

    #[error("Io error occurred")]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Deserialize(#[from] csv::Error),

    #[error("Unable to deserialize line number {line}")]
    DeserializeRecord { source: csv::Error, line: usize },

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

/// A single row in a sample metadata file.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Hash, Eq)]
pub struct SampleMetadata {
    /// The unique identifier for the sample.
    #[serde(alias = "Sample_ID", rename(serialize = "Sample_ID"))]
    pub sample_id: String,
    /// The sample barcode sequence.
    ///
    /// If the sample barcode is present across multiple reads (ex. dual-index, or inline in both reads of a pair),
    /// then the expected barcode bases from each read should be concatenated in the same order as the order of
    /// the reads' FASTQs and read structures given. Non-ACTG characters will be removed.
    #[serde(alias = "Sample_Barcode", rename(serialize = "Sample_Barcode"))]
    pub raw_barcode: BString,

    /// The sanitized `raw_barcode`
    #[serde(skip)]
    pub barcode: BString,

    /// The number of the sample in the samplesheet (corresponds to the row number), starts at 0.
    #[serde(skip)]
    pub ordinal: usize,
}

impl SampleMetadata {
    /// Create a new [`SampleMetadata`] object.
    ///
    /// # Errors
    ///
    /// - [`SampleMetadataError::InvalidBarcode`] if the barcode is invalid
    pub fn new(
        sample_id: String,
        barcode: BString,
        number: usize,
    ) -> Result<Self, SampleMetadataError> {
        let fixed = Self::sanitize_barcode(barcode.as_ref());
        Self::validate_barcode(fixed.as_ref(), &sample_id, None)?;
        Ok(Self { sample_id, raw_barcode: barcode, barcode: fixed, ordinal: number })
    }

    /// Create a new [`SampleMetadata`] object without sanitizing the bases first.
    ///
    /// This is useful for creating barcodes of all `N`s. Other barcode rules still apply,
    /// such as barcodes being greater than length 0.
    ///
    /// # Errors
    ///
    /// - [`SampleMetadataError::InvalidBarcode`] if the barcode is invalid
    pub fn new_allow_invalid_bases(
        sample_id: String,
        barcode: BString,
        number: usize,
    ) -> Result<Self, SampleMetadataError> {
        Ok(Self { sample_id, raw_barcode: barcode.clone(), barcode, ordinal: number })
    }

    /// Run a set of validations on a barcode to ensure that it is well formed.
    ///
    /// # Errors
    ///
    /// - [`SampleMetadataError::InvalidBarcode`] if the barcode is invalid
    pub fn validate_barcode(
        barcode: &BStr,
        id: &str,
        line_number: Option<usize>,
    ) -> Result<(), SampleMetadataError> {
        if barcode.len() == 0 {
            Err(SampleMetadataError::InvalidBarcode {
                id: id.to_owned(),
                barcode: barcode.to_string(),
                reason: ReasonBarcodeInvalid::EmptyString,
                line: ErrorLine(line_number),
            })
        } else {
            Ok(())
        }
    }
    /// Sanitize a barcode sequence to remove any symbols other than `ACTG`.
    pub fn sanitize_barcode(raw_barcode: &BStr) -> BString {
        raw_barcode
            .to_ascii_uppercase()
            .iter()
            .filter(|&b| ALLOWED_BASES.contains(b))
            .copied()
            .collect()
    }

    /// Validate a collection of [`SampleMetadata`] barcodes by checking that all sample barcodes
    /// are >= `min_mismatches` away from each other by hamming distance.
    ///
    /// # Errors
    ///
    /// - [`SampleMetadataError::BarcodeCollision`] if two barcodes are with the allowed distance
    /// - [`SampleMetadataError::UnequalBarcodeLengths`] if two barcodes have unequal length
    pub fn validate_barcode_pairs(
        samples: &[Self],
        min_mismatches: usize,
    ) -> Result<(), SampleMetadataError> {
        for (i, sample) in samples.iter().enumerate() {
            let barcode = &sample.barcode;
            for other in samples.iter().skip(i + 1) {
                // Ensure barcodes have equal length
                if barcode.len() != other.barcode.len() {
                    return Err(SampleMetadataError::UnequalBarcodeLengths {
                        sample_a: sample.sample_id.clone(),
                        barcode_a: sample.raw_barcode.to_string(),
                        sample_b: other.sample_id.clone(),
                        barcode_b: other.raw_barcode.to_string(),
                    });
                }

                let distance = Self::hamming_distance(barcode, &other.barcode);

                // Ensure barcodes don't collide
                if distance < min_mismatches {
                    return Err(SampleMetadataError::BarcodeCollision {
                        sample_a: sample.sample_id.clone(),
                        barcode_a: sample.raw_barcode.to_string(),
                        sample_b: other.sample_id.clone(),
                        barcode_b: other.raw_barcode.to_string(),
                        distance,
                    });
                }
            }
        }

        Ok(())
    }

    /// Hamming distance on slice of bytes.
    ///
    /// Skips length check and will stop comparing after alpha is exhausted.
    fn hamming_distance(alpha: &[u8], beta: &[u8]) -> usize {
        alpha.iter().zip(beta.iter()).map(|(a, b)| if a == b { 0 } else { 1 }).sum()
    }
}

/// Implementation of AsRef trait for [`SampleMetadata`] to convert `SampleMetadata` to `&SampleMetadata`.
impl AsRef<SampleMetadata> for SampleMetadata {
    fn as_ref(&self) -> &SampleMetadata {
        self
    }
}

/// Deserialize a [`SampleMetadata`] file.
///
/// If the `min_mismatch` is provided, the barcodes will be validated to ensure that there are no collisions.
///
/// # Errors
///
/// - [`SampleMetadataError::DuplicateSampleId`]
/// - [`SampleMetadataError::InvalidBarcode`]
/// - [`SampleMetadataError::BarcodeCollision`]
pub fn from_path<P: AsRef<Path>>(
    path: P,
    min_mismatch: Option<usize>,
    undetermined_name: &str,
) -> Result<Vec<SampleMetadata>, SampleMetadataError> {
    let mut reader = csv::ReaderBuilder::new().has_headers(true).delimiter(b',').from_path(path)?;
    let mut ids = HashSet::new();

    let mut records = vec![];
    for (number, record) in reader.deserialize().enumerate() {
        // Note that line numbers a +2 to account for the header and convert to 1-based counting
        let mut record: SampleMetadata = record
            .map_err(|e| SampleMetadataError::DeserializeRecord { source: e, line: number + 2 })?;
        let fixed = SampleMetadata::sanitize_barcode(record.raw_barcode.as_ref());
        if !ids.insert(record.sample_id.clone()) {
            return Err(SampleMetadataError::DuplicateSampleId { id: record.sample_id });
        }
        record.barcode = fixed;
        record.ordinal = number;
        records.push(record);
    }

    if records.is_empty() {
        return Err(SampleMetadataError::ZeroSamples);
    }
    else if records.len() > 1 {
        // Only validate barcodes if we have more than one sample, since we allow an empty
        // barcode on a single sample.
        for sample in records.iter() {
            SampleMetadata::validate_barcode(sample.barcode.as_ref(), &sample.sample_id, Some(sample.ordinal + 2))?;
        }

        if let Some(min_mismatch) = min_mismatch {
            SampleMetadata::validate_barcode_pairs(&records, min_mismatch)?;
        }
    }

    // If we have more than one sample, or we have one sample with an actual barcode
    // then add in the undetermined sample.
    if records.len() > 1 || records[0].barcode.len() > 0 {
        records.push(SampleMetadata::new_allow_invalid_bases(
            String::from(undetermined_name),
            BString::from(vec![b'N'; records[0].barcode.len()]),
            records.len(),
        )?);
    }

    Ok(records)
}

/// Serialize a collection of [`SampleMetadata`] into a file.
///
/// # Errors
///
/// - [`SampleMetadataError::Io`]
pub fn to_path<P: AsRef<Path>, S: AsRef<SampleMetadata> + Serialize, I: IntoIterator<Item = S>>(
    path: P,
    samples: I,
) -> Result<(), SampleMetadataError> {
    let mut writer = csv::WriterBuilder::new().has_headers(true).delimiter(b',').from_path(path)?;
    for s in samples {
        writer.serialize(s)?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {

    use matches::assert_matches;
    use tempfile::tempdir;

    use crate::matcher::UNDETERMINED_NAME;

    use super::*;
    use std::fs;
    use itertools::Itertools;

    #[test]
    fn test_all_valid_data() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        let found_samples = from_path(&input_path, None, UNDETERMINED_NAME)
            .unwrap().into_iter().filter(|s| s.sample_id != UNDETERMINED_NAME).collect_vec();

        assert_eq!(found_samples.len(), 4, "Found wrong number of serialized samples.");
        for (found, expected) in found_samples.iter().zip(samples.iter()) {
            assert_eq!(found, expected);
        }
    }

    #[test]
    fn test_invalid_barcode() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACNN"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        let r = from_path(&input_path, Some(1), UNDETERMINED_NAME);
        assert_matches!(r, Err(SampleMetadataError::UnequalBarcodeLengths { .. }));
    }

    #[test]
    fn test_invalid_barcode_with_line_number() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");

        let mut invalid =
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2).unwrap();
        // overwrite the "good" barcode that will be serialized
        invalid.raw_barcode = BString::from("");

        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1).unwrap(),
            invalid,
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();
        let r = from_path(&input_path, Some(1), UNDETERMINED_NAME);

        assert_matches!(r, Err(SampleMetadataError::InvalidBarcode { .. }));
        if let Err(SampleMetadataError::InvalidBarcode { barcode, id, reason, line }) = r {
            assert_eq!(line.0, Some(4));
            assert_eq!(barcode, String::from(""));
            assert_eq!(id, String::from("Sample3"));
            assert_matches!(reason, ReasonBarcodeInvalid::EmptyString);
        } else {
            panic!("Wrong error returned");
        }
    }

    #[test]
    fn test_duplicate_sample_id() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACGT"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample1"), BString::from("CTGA"), 1).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        assert_matches!(
            from_path(&input_path, None, UNDETERMINED_NAME),
            Err(SampleMetadataError::DuplicateSampleId { .. })
        );
    }

    #[test]
    fn test_barcode_collision() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3).unwrap(),
            SampleMetadata::new(String::from("Sample5"), BString::from("GGGC"), 4).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        assert_matches!(
            from_path(&input_path, Some(2), UNDETERMINED_NAME),
            Err(SampleMetadataError::BarcodeCollision { .. })
        );
    }

    #[test]
    fn test_blank_lines_skipped() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");

        let bytes = "\
Sample_ID,Sample_Barcode
Sample1,ACTG

Sample2,GGGG



";
        fs::write(&input_path, bytes).unwrap();

        let expected = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("GGGG"), 1).unwrap(),
            SampleMetadata::new_allow_invalid_bases(UNDETERMINED_NAME.to_string(), BString::from("NNNN"), 2).unwrap()
        ];

        assert_eq!(from_path(&input_path, Some(1), UNDETERMINED_NAME).unwrap(), expected);
    }

    #[test]
    fn test_blank_records_fail() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");

        let bytes = "\
Sample_ID,Sample_Barcode
Sample1,ACTG

Sample2,GGGG
,

";
        fs::write(&input_path, bytes).unwrap();
        // Failes with barcode empty error
        assert_matches!(
            from_path(&input_path, Some(1), UNDETERMINED_NAME),
            Err(SampleMetadataError::InvalidBarcode { .. })
        );
    }

    #[test]
    fn test_empty_barcode() {
        let r = SampleMetadata::new(String::from("Sample1"), BString::from(""), 0);
        assert_matches!(r, Err(SampleMetadataError::InvalidBarcode { .. }));
        if let Err(SampleMetadataError::InvalidBarcode { barcode, id, reason, line }) = r {
            assert!(line.0.is_none());
            assert_eq!(barcode, String::from(""));
            assert_eq!(id, String::from("Sample1"));
            assert_matches!(reason, ReasonBarcodeInvalid::EmptyString);
        } else {
            panic!("Wrong error returned");
        }
    }
}
