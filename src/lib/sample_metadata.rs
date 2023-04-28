#![forbid(unsafe_code)]
#![allow(clippy::must_use_candidate)]

use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use bstr::{BStr, BString};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::sample_sheet::{ErrorLine, ReasonBarcodeInvalid, SampleSheetError};

/// The bases that are allowed in the [`SampleMetadata::barcode`].
const ALLOWED_BASES: &[u8] = &[b'A', b'C', b'T', b'G'];

/// Metadata about a sample.
///
/// Metadata may be derived from a simple CSV file or a Singular Genomics Sample Sheet.
///
/// See the `update_with()` method to derive the barcode for demultiplexing from either the
/// `raw_barcode` or the `index1`/`index2` fields.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Hash, Eq)]
pub struct SampleMetadata {
    /// The unique identifier for the sample.
    #[serde(alias = "Sample_ID", rename(serialize = "Sample_ID"))]
    pub sample_id: String,

    /// The sample barcode sequence concatenated across all reads in the order the reads are
    /// given.
    #[serde(alias = "Sample_Barcode", rename(serialize = "Sample_Barcode"))]
    pub raw_barcode: Option<BString>,

    /// The lane number.
    #[serde(alias = "Lane", rename(serialize = "Lane"))]
    pub lane: Option<usize>,

    /// The sample barcode in the index1 read
    #[serde(alias = "Index1_Sequence", rename(serialize = "Index1_Sequence"))]
    pub index1: Option<BString>,

    /// The sample barcode in the index2 read
    #[serde(alias = "Index2_Sequence", rename(serialize = "Index2_Sequence"))]
    pub index2: Option<BString>,

    /// The project associated with this the sample
    #[serde(alias = "Project", rename(serialize = "Project"))]
    pub project: Option<String>,

    /// The sanitized `raw_barcode` used for demultiplexing.
    #[serde(skip)]
    pub barcode: BString,

    /// The number of the sample in the samplesheet (corresponds to the row number), starts at 0.
    #[serde(skip)]
    pub ordinal: usize,

    /// The line number in the input in which this sample was defined
    pub line_number: Option<usize>,
}

impl SampleMetadata {
    /// Create a new [`SampleMetadata`] object.
    ///
    /// # Errors
    ///
    /// - [`SampleSheetError::InvalidBarcode`] if the barcode is invalid
    pub fn new(
        sample_id: String,
        barcode: BString,
        number: usize,
        line_number: usize,
    ) -> Result<Self, SampleSheetError> {
        let fixed = Self::sanitize_barcode(barcode.as_ref());
        Self::validate_barcode(fixed.as_ref(), &sample_id, None)?;
        Ok(Self {
            sample_id,
            raw_barcode: Some(barcode),
            index1: None,
            index2: None,
            lane: None,
            project: None,
            barcode: fixed,
            ordinal: number,
            line_number: Some(line_number),
        })
    }

    /// Create a new [`SampleMetadata`] object without sanitizing the bases first.
    ///
    /// This is useful for creating barcodes of all `N`s. Other barcode rules still apply,
    /// such as barcodes being greater than length 0.
    ///
    /// # Errors
    ///
    /// - [`SampleSheetError::InvalidBarcode`] if the barcode is invalid
    pub fn new_allow_invalid_bases(
        sample_id: String,
        barcode: BString,
        number: usize,
        lane: Option<usize>,
    ) -> Result<Self, SampleSheetError> {
        Ok(Self {
            sample_id,
            raw_barcode: Some(barcode.clone()),
            index1: None,
            index2: None,
            lane,
            project: None,
            barcode,
            ordinal: number,
            line_number: None,
        })
    }

    /// Updates the sample metadata.
    ///
    /// 1. Sets the ordinal
    /// 2. Sets the sample barcode to use for demultiplexing.
    ///
    /// If `raw_barcode` is not `None`, then both `index1` and `index2` must be `None`.
    /// Furthermore, if the sample barcode is present across multiple reads (ex. dual-index, or
    /// inline in both reads of a pair), then the barcode bases from each read should be
    /// concatenated in the same order as the order of the reads' FASTQs and read structures given.
    ///
    /// If either `index1` or `index2` is not `None`, then `raw_barcode` must be `None`.  If both
    /// `index1` and `index2` are given, then `barcode` should be the concatenation of the two,
    /// otherwise `barcode` will be the index that is present.
    ///
    /// The final `barcode` will have non-ACTG characters removed.
    pub fn update_with_and_set_demux_barcode(
        mut self,
        ordinal: usize,
        line_number: usize,
    ) -> Result<Self, SampleSheetError> {
        // Update the ordinal and line number
        self.ordinal = ordinal;
        self.line_number = Some(line_number);

        // Update the sample barcode to use for demultiplexing
        let raw_barcode = match (self.raw_barcode.clone(), self.index1.clone(), self.index2.clone())
        {
            // Use the Sample_Barcode
            (Some(raw), None, None) => raw,
            // Concatenate i1 and i2
            (None, Some(i1), Some(i2)) => BString::from([i1, i2].iter().join("")),
            // Just i1
            (None, Some(i1), None) => i1,
            // Just i2
            (None, None, Some(i2)) => i2,
            // No sample barcode
            (None, None, None) => BString::from(""),
            // Both Sample_Barcode and i1/i2, so an error
            (_, i1, i2) => {
                let index_name = match (i1, i2) {
                    (Some(_), Some(_)) => "Index1_Sequence and Index2_Sequence",
                    (Some(_), None) => "Index1_Sequence",
                    (None, Some(_)) => "Index2_Sequence",
                    (None, None) => "Bug",
                };
                return Err(SampleSheetError::BarcodeColumnCollision {
                    index_name: index_name.to_string(),
                    line: line_number,
                });
            }
        };
        self.barcode = SampleMetadata::sanitize_barcode(raw_barcode.as_ref());
        Ok(self)
    }

    /// Run a set of validations on a barcode to ensure that it is well formed.
    ///
    /// # Errors
    ///
    /// - [`SampleSheetError::InvalidBarcode`] if the barcode is invalid
    pub fn validate_barcode(
        barcode: &BStr,
        id: &str,
        line_number: Option<usize>,
    ) -> Result<(), SampleSheetError> {
        if barcode.len() == 0 {
            Err(SampleSheetError::InvalidBarcode {
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
    /// - [`SampleSheetError::BarcodeCollision`] if two barcodes are with the allowed distance
    /// - [`SampleSheetError::UnequalBarcodeLengths`] if two barcodes have unequal length
    pub fn validate_barcode_pairs(
        samples: &[Self],
        min_mismatches: usize,
    ) -> Result<(), SampleSheetError> {
        for (i, sample) in samples.iter().enumerate() {
            let barcode = &sample.barcode;
            for other in samples.iter().skip(i + 1) {
                // Ensure barcodes have equal length
                if barcode.len() != other.barcode.len() {
                    return Err(SampleSheetError::UnequalBarcodeLengths {
                        sample_a: sample.sample_id.clone(),
                        barcode_a: sample.barcode.to_string(),
                        sample_b: other.sample_id.clone(),
                        barcode_b: other.barcode.to_string(),
                    });
                }

                let distance = Self::hamming_distance(barcode, &other.barcode);

                // Ensure barcodes don't collide
                if distance < min_mismatches {
                    return Err(SampleSheetError::BarcodeCollision {
                        sample_a: sample.sample_id.clone(),
                        barcode_a: sample.barcode.to_string(),
                        sample_b: other.sample_id.clone(),
                        barcode_b: other.barcode.to_string(),
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

/// Validates a set of samples ([`SampleMetadata`] objects).
///
/// If the `min_mismatch` is provided, the barcodes will be validated to ensure that there are no
/// collisions.
///
/// If there is more than one sample, or we have one sample with an actual barcode then add in
/// an undetermined sample.
///
/// Also validates that the sample ordinals are in strictly ascending order.
///
/// # Errors
///
/// - [`SampleSheetError::ZeroSamples`]
/// - [`SampleSheetError::ZeroSamplesForLane`]
/// - [`SampleSheetError::DuplicateSampleId`]
/// - [`SampleSheetError::InvalidBarcode`]
/// - [`SampleSheetError::BarcodeCollision`]
/// - [`SampleSheetError::SampleOrdinalsOutOfOrder`]
pub fn validate_samples(
    mut samples: Vec<SampleMetadata>,
    min_mismatch: Option<usize>,
    undetermined_name: &str,
    lanes: &[usize],
) -> Result<Vec<SampleMetadata>, SampleSheetError> {
    if samples.is_empty() {
        if lanes.is_empty() {
            return Err(SampleSheetError::ZeroSamples);
        }
        let lanes_str = lanes.iter().map(|lane| format!("{}", lane)).join(",");
        return Err(SampleSheetError::ZeroSamplesForLane { lanes: lanes_str });
    }

    // Check for duplicate sample identifiers
    // TODO: same sample id, but across lanes
    let mut ids = HashSet::new();
    let mut prev_ordinal = 0;
    for record in &samples {
        // Developer note: leverages `ids` to determine if we are on the first sample (ids is empty)
        if !ids.is_empty() && prev_ordinal >= record.ordinal {
            return Err(SampleSheetError::SampleOrdinalsOutOfOrder);
        }
        prev_ordinal = record.ordinal;
        if !ids.insert(record.sample_id.clone()) {
            return Err(SampleSheetError::DuplicateSampleId { id: record.sample_id.clone() });
        }
    }

    // Support a special case where a single-sample samplesheet with no barcode sequence for the
    // sample indicates a non-demultiplexing run.  When this is the case we will not validate
    // sample barcodes or insert an Undetermined record into the metadata.
    if samples.len() > 1 || samples[0].barcode.len() > 0 {
        // Only validate barcodes if we have more than one sample, since we allow an empty
        // barcode on a single sample.
        for sample in &samples {
            SampleMetadata::validate_barcode(
                sample.barcode.as_ref(),
                &sample.sample_id,
                sample.line_number,
            )?;
        }

        if let Some(min_mismatch) = min_mismatch {
            SampleMetadata::validate_barcode_pairs(&samples, min_mismatch)?;
        }

        // If we have more than one sample, or we have one sample with an actual barcode
        // then add in the undetermined sample.  The ordinal must be larger than the maximum
        // ordinal of an existing sample.
        let undetermined_ordinal = prev_ordinal + 1;
        let lane: Option<usize> = if lanes.len() == 1 { Some(lanes[0]) } else { None };
        samples.push(SampleMetadata::new_allow_invalid_bases(
            String::from(undetermined_name),
            BString::from(vec![b'N'; samples[0].barcode.len()]),
            undetermined_ordinal,
            lane,
        )?);
    }

    Ok(samples)
}

/// Subset the samples to those with the specific set of lanes.  If no lanes are provided, all
/// samples are considered.  Next, aggregate samples with the same `Sample_ID` and barcode
/// combination.
pub fn coelesce_samples(samples: Vec<SampleMetadata>, lanes: &[usize]) -> Vec<SampleMetadata> {
    // subset to just the samples for the given lane
    let samples = if lanes.is_empty() {
        samples
    } else {
        // subset to just the samples for the given lane
        samples
            .iter()
            .filter(|sample| sample.lane.filter(|lane| lanes.contains(lane)).is_some())
            .cloned()
            .collect()
    };

    // Aggregate samples with the same Sample_ID and barcode
    let mut sample_groups: HashMap<(String, BString), Vec<SampleMetadata>> = HashMap::new();
    for sample in &samples {
        let key = (sample.sample_id.clone(), sample.barcode.clone());
        sample_groups.entry(key).or_insert_with(|| Vec::with_capacity(1)).push(sample.clone());
    }

    // Create a new sample per group (per unique Sample_ID/barcode combination)
    let mut samples: Vec<SampleMetadata> = sample_groups
        .into_iter()
        .map(|(_, group)| {
            if group.len() > 1 {
                // Pick the sample with the lowest line number to preserve ordering from the input file
                let sample = group.iter().min_by_key(|s| s.line_number).unwrap().clone();
                // only remove lane information if there was > 1 sample found
                SampleMetadata { lane: None, ..sample }
            } else {
                group[0].clone()
            }
        })
        .collect();

    // Sort by line number to keep the order of samples in the outputs (e.g. metrics) the same
    // as the order in the input (e.g. sample sheet).
    samples.sort_by_key(|sample| sample.line_number);

    // Update the ordinals
    samples = samples
        .iter()
        .enumerate()
        .map(|(index, sample)| SampleMetadata { ordinal: index, ..sample.clone() })
        .collect();

    samples
}

/// Serialize a collection of [`SampleMetadata`] into a file.
///
/// # Errors
///
/// - [`SampleSheetError::Io`]
pub fn to_path<P: AsRef<Path>, S: AsRef<SampleMetadata> + Serialize, I: IntoIterator<Item = S>>(
    path: P,
    samples: I,
) -> Result<(), SampleSheetError> {
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

    use super::*;
    use crate::opts::Opts;
    use crate::sample_sheet::SampleSheet;
    use itertools::Itertools;
    use std::fs;

    #[test]
    fn test_all_valid_data() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1, 3).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2, 4).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3, 5).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        let opts = Opts { sample_metadata: input_path, ..Opts::default() };
        let sample_sheet = SampleSheet::from_path(opts).unwrap();
        let found_samples = sample_sheet
            .samples
            .into_iter()
            .filter(|s| s.sample_id != sample_sheet.opts.undetermined_sample_name)
            .collect_vec();

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
            SampleMetadata::new(String::from("Sample1"), BString::from("ACNN"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1, 3).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2, 4).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3, 5).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        let opts = Opts { sample_metadata: input_path, ..Opts::default() };
        let r = SampleSheet::from_path(opts);

        assert_matches!(r, Err(SampleSheetError::UnequalBarcodeLengths { .. }));
    }

    #[test]
    fn test_invalid_barcode_with_line_number() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");

        let mut invalid =
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2, 4).unwrap();
        // overwrite the "good" barcode that will be serialized
        invalid.raw_barcode = Some(BString::from(""));

        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1, 3).unwrap(),
            invalid,
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3, 5).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();
        let opts = Opts { sample_metadata: input_path, allowed_mismatches: 1, ..Opts::default() };
        let r = SampleSheet::from_path(opts);

        assert_matches!(r, Err(SampleSheetError::InvalidBarcode { .. }));
        if let Err(SampleSheetError::InvalidBarcode { barcode, id, reason, line }) = r {
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
            SampleMetadata::new(String::from("Sample1"), BString::from("ACGT"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample1"), BString::from("CTGA"), 1, 3).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2, 4).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3, 5).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();
        let opts = Opts { sample_metadata: input_path, ..Opts::default() };

        assert_matches!(
            SampleSheet::from_path(opts),
            Err(SampleSheetError::DuplicateSampleId { .. })
        );
    }

    #[test]
    fn test_barcode_collision() {
        let tempdir = tempdir().unwrap();
        let input_path = tempdir.path().join("input.csv");
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("CTGA"), 1, 3).unwrap(),
            SampleMetadata::new(String::from("Sample3"), BString::from("AAAA"), 2, 4).unwrap(),
            SampleMetadata::new(String::from("Sample4"), BString::from("GCGC"), 3, 5).unwrap(),
            SampleMetadata::new(String::from("Sample5"), BString::from("GGGC"), 4, 6).unwrap(),
        ];

        to_path(&input_path, samples.iter()).unwrap();

        let opts = Opts { sample_metadata: input_path, allowed_mismatches: 2, ..Opts::default() };

        assert_matches!(
            SampleSheet::from_path(opts),
            Err(SampleSheetError::BarcodeCollision { .. })
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

        let opts = Opts::default();
        let expected = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("GGGG"), 1, 3).unwrap(),
            SampleMetadata::new_allow_invalid_bases(
                opts.undetermined_sample_name.to_string(),
                BString::from("NNNN"),
                2,
                None,
            )
            .unwrap(),
        ];

        let opts = Opts { sample_metadata: input_path, allowed_mismatches: 1, ..Opts::default() };

        assert_eq!(SampleSheet::from_path(opts).unwrap().samples, expected);
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

        let opts = Opts { sample_metadata: input_path, allowed_mismatches: 1, ..Opts::default() };
        let result = SampleSheet::from_path(opts);
        assert_matches!(
            result,
            Err(SampleSheetError::InvalidBarcode { barcode: _, id: _, reason: _, line: _ })
        );
        // do not validate the line number because blank lines are not counted
    }

    #[test]
    fn test_empty_barcode() {
        let r = SampleMetadata::new(String::from("Sample1"), BString::from(""), 0, 2);
        assert_matches!(r, Err(SampleSheetError::InvalidBarcode { .. }));
        if let Err(SampleSheetError::InvalidBarcode { barcode, id, reason, line }) = r {
            assert!(line.0.is_none());
            assert_eq!(barcode, String::from(""));
            assert_eq!(id, String::from("Sample1"));
            assert_matches!(reason, ReasonBarcodeInvalid::EmptyString);
        } else {
            panic!("Wrong error returned");
        }
    }

    /// Create a new [`SampleMetadata`] object.
    ///
    /// # Errors
    ///
    /// - [`SampleSheetError::InvalidBarcode`] if the barcode is invalid
    pub fn new_sample_meta(
        sample_id: String,
        barcode: BString,
        number: usize,
        line_number: usize,
        lane: usize,
    ) -> Result<SampleMetadata, SampleSheetError> {
        SampleMetadata::new(sample_id, barcode, number, line_number)
            .map(|sample| SampleMetadata { lane: Some(lane), ..sample })
    }

    #[test]
    fn test_coalesce_no_lanes() {
        let samples = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            SampleMetadata::new(String::from("Sample2"), BString::from("GGGG"), 1, 3).unwrap(),
        ];
        let actual = coelesce_samples(samples.clone(), &[]);
        assert_eq!(actual, samples);
    }

    #[test]
    fn test_coalesce_lanes_on_sample_but_no_lanes_given() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 1).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        let actual = coelesce_samples(samples.clone(), &[]);
        assert_eq!(actual, samples);
    }

    #[test]
    fn test_coalesce_different_lanes_on_sample_but_no_lanes_given() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 1).unwrap(),
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        let expected = vec![
            SampleMetadata::new(String::from("Sample1"), BString::from("ACTG"), 0, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        let actual = coelesce_samples(samples, &[]);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_coalesce_different_lanes_on_sample_and_subset_given_lane() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 1).unwrap(),
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        // the lane value is kept on both samples
        let expected = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        let actual = coelesce_samples(samples, &[2]);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_coalesce_different_lanes_on_sample_and_subset_given_lane_and_no_samples_remain() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 1).unwrap(),
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 0, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 1, 3, 2).unwrap(),
        ];
        let actual = coelesce_samples(samples, &[3]);
        assert_eq!(actual, []);
    }

    #[test]
    fn test_coalesce_pick_lowest_ordinal() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 2, 2, 1).unwrap(),
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 1, 1, 1).unwrap(),
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 3, 3, 1).unwrap(),
        ];
        let actual = coelesce_samples(samples, &[]);
        let expected = vec![SampleMetadata::new(
            String::from("Sample1"),
            BString::from("ACTG"),
            0,
            1,
        )
        .unwrap()];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_coalesce_sort_by_ordinal_after_filtering() {
        let samples = vec![
            new_sample_meta(String::from("Sample1"), BString::from("ACTG"), 1, 1, 3).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 4, 4, 2).unwrap(),
            new_sample_meta(String::from("Sample3"), BString::from("TTTT"), 3, 3, 2).unwrap(),
            new_sample_meta(String::from("Sample4"), BString::from("CCCC"), 2, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample5"), BString::from("AAAA"), 5, 5, 1).unwrap(),
        ];
        let actual = coelesce_samples(samples, &[2]);
        let expected = vec![
            new_sample_meta(String::from("Sample4"), BString::from("CCCC"), 0, 2, 2).unwrap(),
            new_sample_meta(String::from("Sample3"), BString::from("TTTT"), 1, 3, 2).unwrap(),
            new_sample_meta(String::from("Sample2"), BString::from("GGGG"), 2, 4, 2).unwrap(),
        ];
        assert_eq!(actual, expected);
    }
}
