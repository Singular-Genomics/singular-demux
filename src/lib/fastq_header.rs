#![forbid(unsafe_code)]
#![allow(clippy::must_use_candidate)]
//! Convert a slice of bytes into a [`FastqHeader`].
//!
//! The [`FastqHeader`] internally holds all data as [`Cow`] types
//! which allow it to be very lightweight.
//!
//! **Note**: The number of fields and their structure is validated when parsing the bytes,
//! but the contents of the fields is not validated.
//!
//! # Grammar
//!
//! The FASTQ headers are expected to conform to the following grammar:
//!
//! ```text
//! <header>         ::= <read-name> | “<read-name> <comment>”
//! <read-name>      ::= <instrument>:<run-number>:<flowcell-id>:<lane>:<tile>:<x-pos>:<y-pos>
//! <instrument>     ::= [a-zA-Z0-9_]+
//! <run-number>     ::= <non-zero-digit>
//! <lane>           ::= <non-zero-digit>
//! <tile>           ::= <non-zero-digit>
//! <x-pos>          ::= <non-zero-digit>
//! <y-pox>          ::= <non-zero-digit>
//! <flowcell-id>    ::= [a-zA-Z0-9]+
//! <comment>        ::= “<info>” | “<info> <other>”
//! <info>           ::= <read-number>:<is-filtered>:<control-number>:<sample-barcode>
//! <read-number>    ::= <digit>
//! <is-filtered>    ::= “Y” | “N”
//! <control-number> ::= “0” | <non-zero-digit>
//! <sample-barcode> ::= [ACGTN]+ | “<sample-barcode>+<sample-barcode>”
//! <other>          ::= .*
//! ```
//!
//! # Usage
//!
//! The primary usage of this library is via the [`FastqHeader::try_from`] method:
//!
//! ```rust
//! use sgdemux_lib::fastq_header::FastqHeader;
//! let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000:NNNN 1:N:0:GCATAAGCTT+GGCGACGGAA";
//! let mut header = FastqHeader::try_from(header.as_bytes()).unwrap();
//! header.set_umi("AAAA".as_bytes());
//! assert_eq!(header.read_name.umi.as_ref().unwrap().as_ref(), "AAAA".as_bytes());
//! ```

use std::{
    borrow::{Borrow, Cow},
    error::Error,
    fmt::Display,
};

use bstr::ByteSlice;

/// A Set of errors that are returned when parsing a FASTQ headers bytes
#[derive(Debug)]
pub enum FastqHeaderError {
    MissingReadName,
    // Errors from parsing ReadName
    MissingInstrument,
    MissingRunNumber,
    MissingFlowCellId,
    MissingLane,
    MissingTile,
    MissingXPos,
    MissingYPos,
    MissingComment,
    // Errors from parsingComment
    MissingInfo,
    // Errors from parsing Info
    MissingReadNumber,
    MissingIsFiltered,
    MissingControlNumber,
    MissingSampleBarcode,
}

impl Error for FastqHeaderError {}
impl Display for FastqHeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// A struct representing a parsed FASTQ header.
#[derive(Default, Debug)]
pub struct FastqHeader<'a> {
    /// The required [`ReadName`] information.
    pub read_name: ReadName<'a>,
    /// The optional set of [`Comment`] fields.
    pub comment: Option<Comment<'a>>,
}

impl<'a> FastqHeader<'a> {
    /// Set the `sample_barcode` on the inner [`Info`] field.
    ///
    /// # Errors
    /// - [`FastqHeaderError::MissingInfo`]
    pub fn set_sample_barcode(&mut self, barcode: &'a [u8]) -> Result<(), FastqHeaderError> {
        if let Some(comment) = self.comment.as_mut() {
            comment.info.sample_barcode.to_mut().resize(barcode.len(), 0);
            comment.info.sample_barcode.to_mut().copy_from_slice(barcode);
        } else {
            return Err(FastqHeaderError::MissingComment);
        }
        Ok(())
    }

    /// Get the [`Info::sample_barcode`] field.
    ///
    /// If the optional [`Comment`] field is missing, `None` is returned.
    pub fn sample_barcode(&self) -> Option<&[u8]> {
        if let Some(comment) = self.comment.as_ref() {
            return Some(comment.info.sample_barcode.borrow());
        }
        None
    }

    /// Get the [`ReadName::umi`] field.
    pub fn umi(&self) -> Option<&[u8]> {
        self.read_name.umi.as_ref().map(std::borrow::Borrow::borrow)
    }

    /// Set the [`ReadName::umi`] field.
    pub fn set_umi(&mut self, umi: &'a [u8]) {
        if let Some(umi_field) = self.read_name.umi.as_mut() {
            umi_field.to_mut().resize(umi.len(), 0);
            umi_field.to_mut().copy_from_slice(umi);
        } else {
            self.read_name.umi = Some(Cow::Owned(umi.to_owned()));
        }
    }

    /// Set the [`Info::read_number`] field.
    ///
    /// The number should be the ASCII value for the number being set. i.e. pass in `b'1'` not `1u8`.
    ///
    /// # Errors
    ///
    /// - [`FastqHeaderError::MissingComment`] if the the comment field is `None`
    pub fn set_read_number(&mut self, number: u8) -> Result<(), FastqHeaderError> {
        if let Some(comment) = self.comment.as_mut() {
            comment.info.read_number.to_mut()[0] = number;
        } else {
            return Err(FastqHeaderError::MissingComment);
        }

        Ok(())
    }

    /// Get the [`Info::read_number`] field.
    ///
    /// If the optional [`Comment`] field is missing, `None` is returned.
    pub fn read_number(&self) -> Option<u8> {
        if let Some(comment) = self.comment.as_ref() {
            return comment.info.read_number.first().copied();
        }
        None
    }

    /// Returns true if the [`Info::control_number`] field is both present and set to a value other than 0.
    pub fn is_control(&self) -> bool {
        if let Some(comment) = &self.comment {
            let control_number: &[u8] = comment.info.control_number.borrow();
            if control_number != [b'0'] {
                return true;
            }
        }
        false
    }

    /// Returns true if the [`Info::is_filtered`] field is both present and set
    pub fn is_failed_quality_filter(&self) -> bool {
        if let Some(comment) = &self.comment {
            let quality_filter_value: &[u8] = comment.info.is_filtered.borrow();
            if quality_filter_value == [b'Y'] {
                return true;
            }
        }
        false
    }

    /// Copy the [`FastqHeader`] to a [`Vec<u8>`].
    pub fn copy_to_vec(&self, dest: &mut Vec<u8>) {
        self.read_name.copy_to_vec(dest);
        if let Some(comment) = &self.comment {
            dest.push(b' ');
            comment.copy_to_vec(dest);
        }
    }
}

/// The [`ReadName`] is a required portion of the FASTQ header and contains a sequence identifier
/// along with information about the run and cluster.
#[derive(Default, Debug)]
pub struct ReadName<'a>
where
    [u8]: ToOwned<Owned = Vec<u8>>,
{
    /// The instrument ID.
    pub instrument: Cow<'a, [u8]>,
    /// The run number on the instrument.
    pub run_number: Cow<'a, [u8]>,
    /// The flow cell ID.
    pub flowcell_id: Cow<'a, [u8]>,
    /// The flow cell lane number.
    pub lane: Cow<'a, [u8]>,
    /// The flow cell tile number.
    pub tile: Cow<'a, [u8]>,
    /// The X coordinate of the cluster.
    pub x_pos: Cow<'a, [u8]>,
    /// The Y coordinate of the cluster.
    pub y_pos: Cow<'a, [u8]>,
    /// The optional UMI sequence. When the sample sheet specifies UMIs, a plus sign separates
    /// each non-contiguous UMI sequence.
    pub umi: Option<Cow<'a, [u8]>>,
}

impl<'a> ReadName<'a> {
    /// Copy the [`Self`] to a vec of bytes delimited by the `:` character.
    pub fn copy_to_vec(&self, dest: &mut Vec<u8>) {
        dest.extend(self.instrument.as_ref());
        dest.push(b':');
        dest.extend(self.run_number.as_ref());
        dest.push(b':');
        dest.extend(self.flowcell_id.as_ref());
        dest.push(b':');
        dest.extend(self.lane.as_ref());
        dest.push(b':');
        dest.extend(self.tile.as_ref());
        dest.push(b':');
        dest.extend(self.x_pos.as_ref());
        dest.push(b':');
        dest.extend(self.y_pos.as_ref());
        if let Some(umi) = &self.umi {
            dest.push(b':');
            dest.extend(umi.as_ref());
        }
    }
}

/// A [`Comment`] represents all the optional fields after the first space character in the
/// FASTQ header.
#[derive(Debug)]
pub struct Comment<'a>
where
    [u8]: ToOwned<Owned = Vec<u8>>,
{
    /// The optional [`Info`] field contains sample specific information.
    pub info: Info<'a>,
    /// The optional other field containing all bytes after the second space character (if present).
    /// These bytes are left unparsed.
    pub other: Option<Cow<'a, [u8]>>,
}

impl<'a> Comment<'a> {
    /// Create a new [`Comment`] object.
    pub fn new(info: Info<'a>, other: Option<Vec<u8>>) -> Self {
        Self { info, other: other.map(Cow::Owned) }
    }
    /// Copy a [`Comment`] into a vec of bytes delimited byte the ` ` character.
    pub fn copy_to_vec(&self, dest: &mut Vec<u8>) {
        self.info.copy_to_vec(dest);
        if let Some(other) = &self.other {
            dest.push(b' ');
            dest.extend(other.as_ref());
        }
    }
}

/// The [`Info`] field contains sample specific information.
#[derive(Debug)]
pub struct Info<'a>
where
    [u8]: ToOwned<Owned = Vec<u8>>,
{
    /// What number read this is.
    pub read_number: Cow<'a, [u8]>,
    /// Indicator of whether the read is filtered.
    pub is_filtered: Cow<'a, [u8]>,
    /// Indicator of whether or not this sample is a control
    pub control_number: Cow<'a, [u8]>,
    /// The Index Read sequence
    ///
    /// - If the sample sheet indicates indexing, the index sequence is appended to the end
    ///   of the read identifier.
    /// - If indexing is not indicated (one sample per lane), the sample number is appended to the
    ///   read identifier.
    pub sample_barcode: Cow<'a, [u8]>,
}

impl<'a> Info<'a> {
    /// Create a new [`Info`] object.
    pub fn new(
        read_number: u8,
        is_filtered: u8,
        control_number: u8,
        sample_barcode: Vec<u8>,
    ) -> Self {
        Self {
            read_number: Cow::Owned(vec![read_number]),
            is_filtered: Cow::Owned(vec![is_filtered]),
            control_number: Cow::Owned(vec![control_number]),
            sample_barcode: Cow::Owned(sample_barcode),
        }
    }

    /// Copy an [`Info`] into a vec of bytes delimited byte the `:` character.
    pub fn copy_to_vec(&self, dest: &mut Vec<u8>) {
        dest.extend(self.read_number.as_ref());
        dest.push(b':');
        dest.extend(self.is_filtered.as_ref());
        dest.push(b':');
        dest.extend(self.control_number.as_ref());
        dest.push(b':');
        dest.extend(self.sample_barcode.as_ref());
    }
}

impl<'a> TryFrom<&'a [u8]> for FastqHeader<'a> {
    type Error = FastqHeaderError;

    /// Convert a slice of bytes into a [`FastqHeader`].
    fn try_from(bytes: &'a [u8]) -> Result<Self, Self::Error> {
        let bytes = bytes.trim_end();
        let mut parts = bytes.splitn(2, |c| *c == b' ');
        let read_name = if let Some(read_name_part) = parts.next() {
            let mut parts = read_name_part.splitn(8, |c| *c == b':');

            let instrument = Cow::from(parts.next().ok_or(FastqHeaderError::MissingInstrument)?);
            let run_number = Cow::from(parts.next().ok_or(FastqHeaderError::MissingRunNumber)?);
            let flowcell_id = Cow::from(parts.next().ok_or(FastqHeaderError::MissingFlowCellId)?);
            let lane = Cow::from(parts.next().ok_or(FastqHeaderError::MissingLane)?);
            let tile = Cow::from(parts.next().ok_or(FastqHeaderError::MissingTile)?);
            let x_pos = Cow::from(parts.next().ok_or(FastqHeaderError::MissingXPos)?);
            let y_pos = Cow::from(parts.next().ok_or(FastqHeaderError::MissingYPos)?);
            let umi = parts.next().map(Cow::from);
            ReadName { instrument, run_number, flowcell_id, lane, tile, x_pos, y_pos, umi }
        } else {
            return Err(FastqHeaderError::MissingReadName);
        };
        let comment = if let Some(comment_parts) = parts.next() {
            let mut parts = comment_parts.splitn(2, |c| *c == b' ');

            let info = if let Some(info_parts) = parts.next() {
                let mut parts = info_parts.splitn(4, |c| *c == b':');

                let read_number =
                    Cow::from(parts.next().ok_or(FastqHeaderError::MissingReadNumber)?);
                let is_filtered =
                    Cow::from(parts.next().ok_or(FastqHeaderError::MissingIsFiltered)?);
                let control_number =
                    Cow::from(parts.next().ok_or(FastqHeaderError::MissingControlNumber)?);
                let sample_barcode =
                    Cow::from(parts.next().ok_or(FastqHeaderError::MissingSampleBarcode)?);
                Info { read_number, is_filtered, control_number, sample_barcode }
            } else {
                return Err(FastqHeaderError::MissingInfo);
            };

            let other = parts.next().map(Cow::from);
            Some(Comment { info, other })
        } else {
            None
        };
        Ok(Self { read_name, comment })
    }
}

#[cfg(test)]
mod tests {
    use super::FastqHeader;

    #[test]
    fn test_all_fields() {
        let instrument = "aass235_9";
        let run_number = "1";
        let flowcell_id = "aa30ZZ";
        let lane = "2";
        let tile = "3";
        let x_pos = "4";
        let y_pos = "5";
        let umi = "ACTG+GCTA";
        let read_number = "0";
        let is_filtered = "N";
        let control_number = "1";
        let sample_barcode = "ACTG+GCTA";
        let other = "This is a long comment section for: other.";
        let header = format!(
            "{}:{}:{}:{}:{}:{}:{}:{} {}:{}:{}:{} {}",
            instrument,
            run_number,
            flowcell_id,
            lane,
            tile,
            x_pos,
            y_pos,
            umi,
            read_number,
            is_filtered,
            control_number,
            sample_barcode,
            other
        );

        let header: FastqHeader = header.as_bytes().try_into().unwrap();
        assert_eq!(instrument.as_bytes(), header.read_name.instrument.as_ref());
        assert_eq!(run_number.as_bytes(), header.read_name.run_number.as_ref());
        assert_eq!(flowcell_id.as_bytes(), header.read_name.flowcell_id.as_ref());
        assert_eq!(lane.as_bytes(), header.read_name.lane.as_ref());
        assert_eq!(tile.as_bytes(), header.read_name.tile.as_ref());
        assert_eq!(x_pos.as_bytes(), header.read_name.x_pos.as_ref());
        assert_eq!(y_pos.as_bytes(), header.read_name.y_pos.as_ref());
        assert_eq!(umi.as_bytes(), header.read_name.umi.unwrap().as_ref());
        assert_eq!(
            read_number.as_bytes(),
            header.comment.as_ref().unwrap().info.read_number.as_ref()
        );
        assert_eq!(
            is_filtered.as_bytes(),
            header.comment.as_ref().unwrap().info.is_filtered.as_ref()
        );
        assert_eq!(
            control_number.as_bytes(),
            header.comment.as_ref().unwrap().info.control_number.as_ref()
        );
        assert_eq!(
            sample_barcode.as_bytes(),
            header.comment.as_ref().unwrap().info.sample_barcode.as_ref()
        );
        assert_eq!(
            other.as_bytes(),
            header.comment.as_ref().unwrap().other.as_ref().unwrap().as_ref()
        );
    }
    #[test]
    fn test_all_fields_without_umi() {
        let instrument = "aass235_9";
        let run_number = "1";
        let flowcell_id = "aa30ZZ";
        let lane = "2";
        let tile = "3";
        let x_pos = "4";
        let y_pos = "5";
        let read_number = "0";
        let is_filtered = "N";
        let control_number = "1";
        let sample_barcode = "ACTG+GCTA";
        let other = "This is a long comment section for: other.";
        let header = format!(
            "{}:{}:{}:{}:{}:{}:{} {}:{}:{}:{} {}",
            instrument,
            run_number,
            flowcell_id,
            lane,
            tile,
            x_pos,
            y_pos,
            read_number,
            is_filtered,
            control_number,
            sample_barcode,
            other
        );

        let header: FastqHeader = header.as_bytes().try_into().unwrap();
        assert_eq!(instrument.as_bytes(), header.read_name.instrument.as_ref());
        assert_eq!(run_number.as_bytes(), header.read_name.run_number.as_ref());
        assert_eq!(flowcell_id.as_bytes(), header.read_name.flowcell_id.as_ref());
        assert_eq!(lane.as_bytes(), header.read_name.lane.as_ref());
        assert_eq!(tile.as_bytes(), header.read_name.tile.as_ref());
        assert_eq!(x_pos.as_bytes(), header.read_name.x_pos.as_ref());
        assert_eq!(y_pos.as_bytes(), header.read_name.y_pos.as_ref());
        assert!(header.read_name.umi.is_none());
        assert_eq!(
            read_number.as_bytes(),
            header.comment.as_ref().unwrap().info.read_number.as_ref()
        );
        assert_eq!(
            is_filtered.as_bytes(),
            header.comment.as_ref().unwrap().info.is_filtered.as_ref()
        );
        assert_eq!(
            control_number.as_bytes(),
            header.comment.as_ref().unwrap().info.control_number.as_ref()
        );
        assert_eq!(
            sample_barcode.as_bytes(),
            header.comment.as_ref().unwrap().info.sample_barcode.as_ref()
        );
        assert_eq!(
            other.as_bytes(),
            header.comment.as_ref().unwrap().other.as_ref().unwrap().as_ref()
        );
    }

    #[test]
    fn test_good_fastq_header() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 1:N:0:GCATAAGCTT+GGCGACGGAA Random extra:fields!!.  ";
        assert!(FastqHeader::try_from(header.as_bytes()).is_ok());
    }

    #[test]
    fn test_good_fastq_header_missing_other() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 1:N:0:GCATAAGCTT+GGCGACGGAA";
        assert!(FastqHeader::try_from(header.as_bytes()).is_ok());
    }

    #[test]
    fn test_good_fastq_header_missing_comment() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000";
        assert!(FastqHeader::try_from(header.as_bytes()).is_ok());
    }

    #[test]
    fn test_bad_fastq_header() {
        let header = "Hello World";
        assert!(FastqHeader::try_from(header.as_bytes()).is_err());
    }

    #[test]
    fn test_bad_fastq_header_invalid_read_name() {
        let header =
            "H00233:4:AAAFGW3HV:1:1101:1000 1:N:0:GCATAAGCTT+GGCGACGGAA Random extra:fields!!.  ";
        assert!(FastqHeader::try_from(header.as_bytes()).is_err());
    }

    #[test]
    fn test_bad_fastq_header_invalid_info() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 1:0:GCATAAGCTT+GGCGACGGAA Random extra:fields!!.  ";
        assert!(FastqHeader::try_from(header.as_bytes()).is_err());
    }

    #[test]
    fn test_bad_fastq_header_no_info_with_other() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 Random extra:fields!!.  ";
        assert!(FastqHeader::try_from(header.as_bytes()).is_err());
    }

    #[test]
    fn test_update_umi() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000:NNNN 1:N:0:GCATAAGCTT+GGCGACGGAA";
        let mut header = FastqHeader::try_from(header.as_bytes()).unwrap();
        header.set_umi("AAAA".as_bytes());
        assert_eq!(header.read_name.umi.as_ref().unwrap().as_ref(), "AAAA".as_bytes());
    }

    #[test]
    fn test_update_sample_barcode() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 1:N:0:GCATAAGCTT+GGCGACGGAA";
        let mut header = FastqHeader::try_from(header.as_bytes()).unwrap();
        header.set_sample_barcode("AAAA+CCCC".as_bytes()).unwrap();
        assert_eq!(
            header.comment.as_ref().unwrap().info.sample_barcode.as_ref(),
            "AAAA+CCCC".as_bytes()
        );
    }

    #[test]
    fn test_round_trip_serialize() {
        let header = "H00233:4:AAAFGW3HV:1:1101:59586:1000 1:N:0:GCATAAGCTT+GGCGACGGAA";
        let new_header = FastqHeader::try_from(header.as_bytes()).unwrap();
        let mut dest = vec![];
        new_header.copy_to_vec(&mut dest);
        assert_eq!(header.as_bytes(), dest.as_slice());
    }
}
