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

    pub fn is_match(&self) -> bool {
        matches!(self, Self::Match { .. })
    }

    pub fn is_no_match(&self) -> bool {
        !self.is_match()
    }
}

#[derive(ArgEnum, Debug, Clone, Copy, PartialEq, Eq)]
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
        barcode: &[u8],
        min_mismatches: usize,
        max_mismatches: usize,
    ) -> impl Iterator<Item = (Vec<u8>, usize)> {
        let free_ns = 0;
        // xs.iter().combinations(n) returns [] if n > xs.len();
        let max_mismatches = max_mismatches.min(barcode.len());
        (0..barcode.len())
            .combinations(max_mismatches) // generate positions to insert mismatches
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
                    .filter(|m| should_reject_delta(m.hamming_dist - best.hamming_dist, min_delta))
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

#[inline]
pub fn should_reject_delta(delta: usize, min_delta: usize) -> bool {
    // Delta is the difference best distance and next best distance.
    // Accept delta > min_delta and reject delta <=  min_delta
    delta <= min_delta
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
        || should_reject_delta(next_best_dist - best_dist, min_delta)
    {
        MatchResult::NoMatch { barcode }
    } else {
        MatchResult::Match { sample_index: best_sample, hamming_dist: best_dist, barcode }
    }
}

type DistanceFunction = dyn Fn(&[u8], &[u8], usize) -> usize;

fn find_independently_cached(
    barcode: Vec<u8>,
    samples: &[SampleMetadata],
    max_mismatches: usize,
    min_delta: usize,
    free_ns: usize,
    distance: &DistanceFunction,
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
    use ahash::AHashMap;
    use bstr::{BString, B};

    use super::{hamming_distance, PreComputeMatcher, PrecomputedMatch};
    use crate::sample_metadata::SampleMetadata;

    pub(super) fn create_sample(barcode: BString) -> SampleMetadata {
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

        let expected = AHashMap::from_iter(IntoIterator::into_iter([
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

        // Check degenerate case (possible in tests) where permutations work with max_matches > barcode len
        let l = PreComputeMatcher::build_map(&samples, 3, 0);
        let r = PreComputeMatcher::build_map(&samples, 10, 0);
        assert_eq!(l, r);
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

    #[test]
    fn test_hamming_dist_free_ns() {
        let dist = hamming_distance(b"GATTACN", b"GATTACN", 0);
        assert_eq!(dist, 1);

        let dist = hamming_distance(b"GATTACA", b"GATTACN", 0);
        assert_eq!(dist, 1);

        let dist = hamming_distance(b"GATTACN", b"GATTACA", 0);
        assert_eq!(dist, 1);

        let dist = hamming_distance(b"GATTACN", b"GATTACN", 1);
        assert_eq!(dist, 0);

        let dist = hamming_distance(b"GATTACA", b"GATTACN", 1);
        assert_eq!(dist, 0);

        let dist = hamming_distance(b"GATTACN", b"GATTACA", 1);
        assert_eq!(dist, 0);
    }
}

#[cfg(test)]
mod test_matches {
    use crate::sample_metadata::SampleMetadata;

    use super::test::create_sample;
    use super::{
        find_independently, hamming_distance, should_reject_delta, MatchResult, Matcher,
        PreComputeMatcher,
    };

    pub(super) fn create_samples(barcodes: &[&str]) -> Vec<SampleMetadata> {
        // need unique ordinals for appropriate test output from PreComputeMatcher::build_map
        barcodes
            .into_iter()
            .enumerate()
            .map(|(i, b)| {
                let mut s = create_sample(bstr::BString::from(b.clone()));
                s.ordinal = i;
                s
            })
            .collect()
    }

    // Helper fun to make find_independent call more concise
    pub(super) fn find_ind(
        barcode: Vec<u8>,
        samples: &[SampleMetadata],
        max_mismatches: usize,
        min_delta: usize,
    ) -> MatchResult {
        find_independently(barcode, samples, max_mismatches, min_delta, 0, hamming_distance)
    }

    // A really pedantic test because prior impls had a bug.
    #[test]
    fn test_reject_delta() {
        // e.g.) min_delta == 0 is the most permissive setting.
        // When min_delta == 0:
        //  if delta <= 0, reject
        //  if delta > 0 (i.e. >= 1), accept
        assert!(should_reject_delta(0, 0)); // reject
        assert!(!should_reject_delta(1, 0)); // accept
        assert!(!should_reject_delta(2, 0));

        // min_delta == 3
        assert!(should_reject_delta(1, 3)); // reject
        assert!(should_reject_delta(3, 3)); // reject
        assert!(!should_reject_delta(4, 3)); // accept
    }

    #[test]
    fn test_find_ind_match() {
        // Small example
        let sample_b1 = "AAAA";
        let sample_b2 = "AACC";
        let samples = create_samples(&[sample_b1, sample_b2]);

        // Observed barcode with best and next distance 0 and 2; delta is 2;
        let observed_barcode = b"AAAA".to_vec();
        assert_eq!(hamming_distance(&observed_barcode, sample_b1.as_bytes(), 0), 0);
        assert_eq!(hamming_distance(&observed_barcode, sample_b2.as_bytes(), 0), 2);

        // Find matches
        let max_mismatches = 2;
        let min_delta = 1;
        let ind_result = find_ind(observed_barcode.clone(), &samples, max_mismatches, min_delta);
        assert!(ind_result.is_match());
    }

    #[test]
    fn test_precompute_match() {
        // Small example
        let sample_b1 = "AAAA";
        let sample_b2 = "AACC";
        let samples = create_samples(&[sample_b1, sample_b2]);

        // Observed barcode with best and next distance 0 and 2; delta is 2;
        let observed_barcode = b"AAAA".to_vec();
        // Find matches
        let max_mismatches = 2;
        let min_delta = 1;
        let matcher = PreComputeMatcher::new(&samples, max_mismatches, min_delta, 0);
        let pc_result = matcher.find(observed_barcode.clone());
        assert!(pc_result.is_match());
    }

    #[test]
    fn test_find_ind_no_match() {
        let sample_b1 = "AAAA";
        let sample_b2 = "AACC";
        let samples: Vec<SampleMetadata> = create_samples(&[sample_b1, sample_b2]);

        let observed_barcode = b"AAAG".to_vec();
        // best and next bast hamming distance; delta is 1.
        assert_eq!(hamming_distance(&observed_barcode, sample_b1.as_bytes(), 0), 1);
        assert_eq!(hamming_distance(&observed_barcode, sample_b2.as_bytes(), 0), 2);

        let max_mismatches = 2;
        let min_delta = 1;
        let ind_result = find_ind(observed_barcode.clone(), &samples, max_mismatches, min_delta);
        assert!(ind_result.is_no_match());

        let matcher = PreComputeMatcher::new(&samples, max_mismatches, min_delta, 0);
        let pc_result = matcher.find(observed_barcode.clone());
        assert!(pc_result.is_no_match());

        assert_eq!(ind_result, pc_result);
    }

    #[test]
    fn test_precompute_no_match() {
        let sample_b1 = "AAAA";
        let sample_b2 = "AACC";
        let samples = create_samples(&[sample_b1, sample_b2]);

        let observed_barcode = b"AAAG".to_vec();

        let max_mismatches = 2;
        let min_delta = 1;
        let matcher = PreComputeMatcher::new(&samples, max_mismatches, min_delta, 0);
        let pc_result = matcher.find(observed_barcode.clone());
        assert!(pc_result.is_no_match());
    }
}

#[cfg(test)]
mod test_readme_examples {
    use crate::sample_metadata::SampleMetadata;

    use super::test_matches::{create_samples, find_ind};
    use super::{hamming_distance, Matcher, PreComputeMatcher};

    struct SimpleTest {
        should_match: bool,
        max_mismatches: usize,
        min_delta: usize,
        b1_dist: usize,
        b2_dist: usize,
    }

    impl SimpleTest {
        fn run(self) {
            assert!(self.b2_dist >= self.b1_dist);

            // generate some barcodes
            let sample_b1 = "A".repeat(self.b2_dist);
            let sample_b2 = "T".repeat(self.b2_dist);

            let samples: Vec<SampleMetadata> = create_samples(&[&sample_b1, &sample_b2]);

            let mut observed_barcode = vec![b'G'; self.b1_dist];
            observed_barcode.append(&mut vec![b'A'; sample_b1.len() - self.b1_dist]);

            let ind_result =
                find_ind(observed_barcode.clone(), &samples, self.max_mismatches, self.min_delta);
            let matcher = PreComputeMatcher::new(&samples, self.max_mismatches, self.min_delta, 0);
            let pc_result = matcher.find(observed_barcode.clone());

            // check that the generated observed barcode has given distances.
            assert_eq!(hamming_distance(&observed_barcode, sample_b1.as_bytes(), 0), self.b1_dist);
            assert_eq!(hamming_distance(&observed_barcode, sample_b2.as_bytes(), 0), self.b2_dist);

            dbg!(sample_b1);
            dbg!(sample_b2);
            dbg!(String::from_utf8(observed_barcode).unwrap());

            if self.should_match {
                assert!(ind_result.is_match(), "independent matcher failed to find match");
                assert!(
                    pc_result.is_match(),
                    "independent matcher found match but precompute matcher failed to find match"
                );
            } else {
                assert!(ind_result.is_no_match(), "independent matcher failed to reject match");
                assert!(
                pc_result.is_no_match(),
                "independent matcher rejected match but precompute matcher failed to reject match"
            );
            }
        }
    }

    fn test_readme_ex(should_match: bool, best_dist: usize, next_dist: usize) {
        SimpleTest {
            should_match,
            max_mismatches: 3,
            min_delta: 1,
            b1_dist: best_dist,
            b2_dist: next_dist,
        }
        .run();
    }

    #[test]
    fn ex1() {
        test_readme_ex(false, 2, 3)
    }

    #[test]
    fn ex2() {
        test_readme_ex(true, 1, 3)
    }

    #[test]
    fn ex3() {
        test_readme_ex(true, 0, 2)
    }

    #[test]
    fn ex4() {
        test_readme_ex(false, 0, 1)
    }

    #[test]
    fn ex5() {
        test_readme_ex(true, 3, 6)
    }

    #[test]
    fn ex6() {
        test_readme_ex(false, 4, 6)
    }

    #[test]
    fn ex7() {
        test_readme_ex(false, 2, 2)
    }
}
