//! Create a lookup table that contains all valid matches that are possible given the expected
//! barcodes, the maximum hamming distance allowed from an observed barcode to the expected, and
//! the minimum hamming distance allowed from an observed barcode to any other expected barcode.

use std::cell::RefCell;

use ahash::{AHashMap, AHashSet};
use anyhow::anyhow;
use cached::SizedCache;
use clap::{ArgEnum, PossibleValue};
use itertools::Itertools;

use crate::sample_metadata::SampleMetadata;

/// The bases that are allowed in the [`SampleMetadata::barcode`].
const ALLOWED_BASES: &[u8] = &[b'A', b'C', b'T', b'G', b'N'];

/// The name given to the "undetermined" sample
pub const UNDETERMINED_NAME: &str = "Undetermined";

thread_local! (
    /// The barcode cache used by the `CachedHammingDistanceMatcher`.
    static CACHE: RefCell<SizedCache<Vec<u8>, MatchResult>> =
    RefCell::new(SizedCache::with_size(100_000))
);

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub enum MatchResult {
    Match { sample_index: usize, hamming_dist: usize, barcode: Vec<u8> },
    NoMatch { barcode: Vec<u8> },
}

impl MatchResult {
    pub fn barcode(&self) -> &Vec<u8> {
        match self {
            MatchResult::Match { barcode, .. } | MatchResult::NoMatch { barcode } => barcode,
        }
    }
}

#[derive(ArgEnum, Debug, Clone, Copy, PartialEq)]
pub enum MatcherKind {
    CachedHammingDistance,
    PreCompute,
}

impl MatcherKind {
    pub fn possible_values<'a>() -> impl Iterator<Item = PossibleValue<'a>> {
        MatcherKind::value_variants().iter().filter_map(ArgEnum::to_possible_value)
    }
}

impl std::str::FromStr for MatcherKind {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        for variant in Self::value_variants() {
            if variant.to_possible_value().unwrap().matches(s, false) {
                return Ok(*variant);
            }
        }
        Err(anyhow!("Invalid variant: {}", s))
    }
}

/// The base trait for all matching algorithms.
pub trait Matcher {
    fn find(&self, barcode: Vec<u8>) -> MatchResult;
}

// Matches based on hamming distance with a cache
pub struct CachedHammingDistanceMatcher<'a> {
    pub samples: &'a [SampleMetadata],
    pub max_mismatches: usize,
    pub min_delta: usize,
    pub free_ns: usize,
}

impl<'a> CachedHammingDistanceMatcher<'a> {
    pub fn new(
        samples: &'a [SampleMetadata],
        max_mismatches: usize,
        min_delta: usize,
        free_ns: usize,
    ) -> Self {
        Self { samples, max_mismatches, min_delta, free_ns }
    }
}

impl<'a> Matcher for CachedHammingDistanceMatcher<'a> {
    fn find(&self, barcode: Vec<u8>) -> MatchResult {
        find_independently_cached(
            barcode,
            self.samples,
            self.max_mismatches,
            self.min_delta,
            self.free_ns,
            &hamming_distance,
        )
    }
}

/// Helper struct to track the hamming distance from the expected sample barcode to the permuted observed barcode.
#[derive(Debug, Hash, PartialEq, Eq, Copy, Clone)]
pub struct PrecomputedMatch {
    /// The index of the sample the permuted barcode best matches
    sample: usize,
    /// The hamming distance from the permuted to expected barcode
    hamming_dist: usize,
}

// Precompute allowable matches
pub struct PreComputeMatcher<'a> {
    pub samples: &'a [SampleMetadata],
    pub max_mismatches: usize,
    pub min_delta: usize,
    pub samples_lookup: AHashMap<Vec<u8>, PrecomputedMatch>,
    pub free_ns: usize,
}

impl<'a> Matcher for PreComputeMatcher<'a> {
    fn find(&self, barcode: Vec<u8>) -> MatchResult {
        return match self.samples_lookup.get(&barcode) {
            Some(sample) => MatchResult::Match {
                sample_index: sample.sample,
                hamming_dist: sample.hamming_dist,
                barcode,
            },
            None => MatchResult::NoMatch { barcode },
        };
    }
}

impl<'a> PreComputeMatcher<'a> {
    pub fn new(
        samples: &'a [SampleMetadata],
        max_mismatches: usize,
        min_delta: usize,
        free_ns: usize,
    ) -> Self {
        let samples_lookup = Self::build_map(samples, max_mismatches, min_delta);
        Self { samples, max_mismatches, min_delta, samples_lookup, free_ns }
    }

    /// Build a map of all observable barcodes that are possible given the `max_mismatches` and `min_delta` for the expected
    /// set of barcodes in the [`SampleMetadata`].
    pub fn build_map(
        samples: &[SampleMetadata],
        max_mismatches: usize,
        min_delta: usize,
    ) -> AHashMap<Vec<u8>, PrecomputedMatch> {
        let mut map = AHashMap::new();

        // Fill in all observable barcodes that could generate a match
        for sample in samples.iter() {
            for (read, hamming_dist) in Self::all_permutations(&sample.barcode, 0, max_mismatches) {
                let matches = map.entry(read).or_insert_with(AHashSet::new);
                matches.insert(PrecomputedMatch { sample: sample.ordinal, hamming_dist });
            }
        }

        // Then add higher mismatch "matches" that could disallow a match
        for sample in samples.iter() {
            for (read, hamming_dist) in Self::all_permutations(
                &sample.barcode,
                max_mismatches + 1,
                max_mismatches + min_delta - 1,
            ) {
                if let Some(matches) = map.get_mut(&read) {
                    matches.insert(PrecomputedMatch { sample: sample.ordinal, hamming_dist });
                }
            }
        }

        // Then remove from the map entries that are not ok
        let mut result = AHashMap::new();
        for (observed_barcode, matches) in map {
            if let Some(sample) = Self::pick_matching_sample_barcode(&matches, min_delta) {
                result.insert(observed_barcode, sample);
            }
        }
        result
    }

    /// Generate all permutations of a string between min and max mismatches inclusive
    fn all_permutations(
        barcode: &'_ [u8],
        min_mismatches: usize,
        max_mismatches: usize,
    ) -> impl Iterator<Item = (Vec<u8>, usize)> + '_ {
        let free_ns = 0;
        (0..barcode.len())
            .combinations(max_mismatches)
            .flat_map(move |locations| {
                let mut this_barcode = barcode.iter().map(|c| vec![*c]).collect::<Vec<Vec<u8>>>();
                for location in locations {
                    this_barcode[location] = ALLOWED_BASES.to_vec();
                }
                this_barcode.into_iter().multi_cartesian_product().filter_map(move |bc| {
                    let dist = hamming_distance(barcode, &bc, 1);
                    if dist >= min_mismatches {
                        Some((bc, dist))
                    } else {
                        None
                    }
                })
            })
            .sorted()
            .dedup()
            .flat_map(move |(permuted_bc, _orig_dist)| {
                // permute free Ns
                (0..permuted_bc.len()).combinations(free_ns).flat_map(move |locations| {
                    let mut this_barcode =
                        permuted_bc.iter().map(|c| vec![*c]).collect::<Vec<Vec<u8>>>();
                    for location in locations {
                        this_barcode[location] = vec![b'N'];
                    }
                    this_barcode.into_iter().multi_cartesian_product().map(move |bc| {
                        let dist = hamming_distance(barcode, &bc, free_ns);
                        // an N could fill a permuted spot so the dists won't line up
                        (bc, dist)
                    })
                })
            })
            .sorted()
            .dedup()
    }

    /// Is this match ok, and if so which sample barcode matches best.
    fn pick_matching_sample_barcode(
        matches: &AHashSet<PrecomputedMatch>,
        min_delta: usize,
    ) -> Option<PrecomputedMatch> {
        let best_by_distance = matches.iter().map(|m| m.hamming_dist).min();
        if let Some(best_by_distance) = best_by_distance {
            let best_hits =
                matches.iter().filter(|m| m.hamming_dist == best_by_distance).collect::<Vec<_>>();
            if best_hits.len() == 1 {
                let best = best_hits[0];
                if matches
                    .iter()
                    .filter(|m| m.hamming_dist <= best.hamming_dist + min_delta)
                    .count()
                    == 1
                {
                    return Some(*best);
                }
            }
        }
        None
    }
}

/// Finds the best barcode barcode-by-barcode
fn find_independently(
    barcode: Vec<u8>,
    samples: &[SampleMetadata],
    max_mismatches: usize,
    min_delta: usize,
    free_ns: usize,
    distance: impl Fn(&[u8], &[u8], usize) -> usize,
) -> MatchResult {
    let mut best_sample: usize = samples.len();
    let mut best_dist: usize = usize::MAX;
    let mut next_best_dist: usize = usize::MAX;
    for (i, sample) in samples.iter().enumerate() {
        let sample_barcode = &sample.barcode;
        let dist = distance(&barcode, sample_barcode, free_ns);
        if dist < best_dist {
            best_sample = i;
            next_best_dist = best_dist;
            best_dist = dist;
        } else if dist < next_best_dist {
            next_best_dist = dist;
        }
    }
    if best_sample == samples.len()
        || best_dist > max_mismatches
        || next_best_dist - best_dist < min_delta
    {
        MatchResult::NoMatch { barcode }
    } else {
        MatchResult::Match { sample_index: best_sample, hamming_dist: best_dist, barcode }
    }
}

fn find_independently_cached(
    barcode: Vec<u8>,
    samples: &[SampleMetadata],
    max_mismatches: usize,
    min_delta: usize,
    free_ns: usize,
    distance: &dyn Fn(&[u8], &[u8], usize) -> usize,
) -> MatchResult {
    CACHE.with(|cache| {
        // check the cache
        let c = &mut *cache.borrow_mut();
        let res: Option<&MatchResult> = cached::Cached::cache_get(c, &barcode);
        if let Some(res) = res {
            return res.clone();
        }

        let result =
            find_independently(barcode, samples, max_mismatches, min_delta, free_ns, distance);
        cached::Cached::cache_set(c, result.barcode().clone(), result.clone());
        result
    })
}

/// Hamming distance on slice of bytes.
///
/// Skips length check and will stop comparing after alpha is exhausted.
/// N's are not treated specially given this is a private method and is only expected to compare
/// "validated" barcodes that contain the [`ALLOWED_BASES`].
fn hamming_distance(alpha: &[u8], beta: &[u8], free_ns: usize) -> usize {
    let mut ns_seen = 0;
    let mut dist = 0;
    for (a, b) in alpha.iter().zip(beta.iter()) {
        if (*a == b'N' || *b == b'N') && ns_seen < free_ns {
            // Ignore the first `free_ns`
            ns_seen += 1;
        } else if a != b || *a == b'N' || *b == b'N' {
            // N-N counts as mismatch
            dist += 1;
        }
    }
    dist
}

#[cfg(test)]
mod test {
    use std::array::IntoIter;

    use ahash::AHashMap;
    use bstr::{BString, B};

    use crate::{
        matcher::{hamming_distance, PreComputeMatcher, PrecomputedMatch},
        sample_metadata::SampleMetadata,
    };

    fn create_sample(barcode: BString) -> SampleMetadata {
        SampleMetadata::new(String::from("Sample"), barcode, 0, 2).unwrap()
    }

    #[ignore]
    #[test]
    #[allow(clippy::from_iter_instead_of_collect, clippy::similar_names)]
    fn test_permutations_large() {
        let samples = vec![
            create_sample(BString::from("AAAAAAAAAAAATTTTTTTTTTTT")),
            create_sample(BString::from("TTTTTTTTTTTTAAAAAAAAAAAA")),
            create_sample(BString::from("TTTTTTTTTTTTTTTTTTTTTTTT")),
            create_sample(BString::from("GGGGGGGGGGGGAAAAAAAAAAAA")),
            create_sample(BString::from("CCCCCCCCCCCCAAAAAAAAAAAA")),
            create_sample(BString::from("AAAAAAAAAAAAGGGGGGGGGGGG")),
            create_sample(BString::from("AAAAAAAAAAAACCCCCCCCCCCC")),
            create_sample(BString::from("GGGGGGGGGGGGCCCCCCCCCCCC")),
            // create_sample(BString::from("CCCCCCCCCCCCGGGGGGGGGGGG")),
            // create_sample(BString::from("GGGGGGGGGGGGTTTTTTTTTTTT")),
            // create_sample(BString::from("CCCCCCCCCCCCTTTTTTTTTTTT")),
            // create_sample(BString::from("TTTTTTTTTTTTGGGGGGGGGGGG")),
            // create_sample(BString::from("TTTTTTTTTTTTCCCCCCCCCCCC")),
            // create_sample(BString::from("ACTGACTGACTGTTTTTTTTTTTT")),
            // create_sample(BString::from("TTTTTTTTTTTTACTGACTGACTG")),
            // create_sample(BString::from("ACTGACTGACTGAAAAGGGGCCCC")),
        ];

        let max_mismatches = 3;
        let min_delta = 2;
        let found = PreComputeMatcher::build_map(&samples, max_mismatches, min_delta);
        eprintln!("{} possible barcodes", found.len());

        // let expected = AHashMap::from_iter(IntoIter::new([
        //     (B("AAG").to_vec(), &sample1),
        //     (B("AGA").to_vec(), &sample1),
        //     (B("GAA").to_vec(), &sample1),
        //     (B("AAC").to_vec(), &sample1),
        //     (B("ACA").to_vec(), &sample1),
        //     (B("CAA").to_vec(), &sample1),
        //     (B("AAA").to_vec(), &sample1),
        //     (B("TTG").to_vec(), &sample2),
        //     (B("TGT").to_vec(), &sample2),
        //     (B("GTT").to_vec(), &sample2),
        //     (B("TTC").to_vec(), &sample2),
        //     (B("TCT").to_vec(), &sample2),
        //     (B("CTT").to_vec(), &sample2),
        //     (B("TTT").to_vec(), &sample2),
        // ]));

        // assert_eq!(expected, found);
    }

    #[test]
    #[allow(clippy::from_iter_instead_of_collect, clippy::similar_names)]
    fn test_permutations_small() {
        let mut sample1 = create_sample(BString::from("AAA"));
        sample1.ordinal = 0;
        let mut sample2 = create_sample(BString::from("TTT"));
        sample2.ordinal = 1;

        let samples = vec![sample1.clone(), sample2.clone()];
        let max_mismatches = 1;
        let min_delta = 2;

        let expected = AHashMap::from_iter(IntoIter::new([
            (B("AAG").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("AGA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("GAA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("AAC").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("ACA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("CAA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("AAN").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("ANA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("NAA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 1 }),
            (B("AAA").to_vec(), PrecomputedMatch { sample: 0, hamming_dist: 0 }),
            (B("TTG").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TGT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("GTT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TTC").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TCT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("CTT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TTN").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TNT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("NTT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 1 }),
            (B("TTT").to_vec(), PrecomputedMatch { sample: 1, hamming_dist: 0 }),
        ]));

        let found = PreComputeMatcher::build_map(&samples, max_mismatches, min_delta);
        assert_eq!(expected, found);
    }

    #[test]
    fn test_hamming_dist_no_mismatches() {
        let dist = hamming_distance(b"GATTACA", b"GATTACA", 0);
        assert_eq!(dist, 0);
    }

    #[test]
    fn test_hamming_dist_two_mismatches() {
        let dist = hamming_distance(b"GATTACA", b"GACCACA", 0);
        assert_eq!(dist, 2);
    }

    #[test]
    fn test_hamming_dist_no_calls() {
        let dist = hamming_distance(b"GATTACA", b"GANNACA", 0);
        assert_eq!(dist, 2);
    }

    #[test]
    fn test_hamming_dist_all_mismatches() {
        let dist = hamming_distance(b"GATTACA", b"CTAATGT", 0);
        assert_eq!(dist, 7);
    }
}
