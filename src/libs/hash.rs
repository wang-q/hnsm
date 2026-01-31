use itertools::Itertools;
use minimizer_iter::MinimizerBuilder;
use std::iter::FromIterator;

// These codes were adapted from https://curiouscoding.nl/posts/fast-minimizers/
pub trait Hasher: Clone {
    fn hash(&self, t: &[u8]) -> u64;
    fn hash_kmers(&mut self, k: usize, t: &[u8]) -> Vec<u64> {
        t.windows(k).map(|kmer| self.hash(kmer)).collect::<Vec<_>>()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FxHash;
impl Hasher for FxHash {
    fn hash(&self, t: &[u8]) -> u64 {
        fxhash::hash64(t)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MurmurHash3;
impl Hasher for MurmurHash3 {
    fn hash(&self, t: &[u8]) -> u64 {
        murmurhash3::murmurhash3_x64_128(t, 42).0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RapidHash;
impl Hasher for RapidHash {
    fn hash(&self, t: &[u8]) -> u64 {
        rapidhash::rapidhash(t)
    }
}

pub trait Minimizer {
    /// The absolute positions of all minimizers in the text.
    fn minimizer(&mut self, text: &[u8]) -> Vec<(u64, usize)>;
    fn mins(&mut self, text: &[u8]) -> Vec<u64>;
}

pub struct JumpingMinimizer<H = FxHash> {
    pub w: usize,
    pub k: usize,
    pub hasher: H,
}

impl<H: Hasher> Minimizer for JumpingMinimizer<H> {
    fn minimizer(&mut self, text: &[u8]) -> Vec<(u64, usize)> {
        let mut minimizers = Vec::new();

        // Precompute hashes of all k-mers.
        let hashes = self.hasher.hash_kmers(self.k, text);

        let mut start = 0;
        while start < hashes.len() - self.w {
            // Position_min returns the position of the leftmost minimal hash.
            let min_pos = start
                + hashes[start..start + self.w]
                    .iter()
                    .position_min()
                    .expect("w > 0");
            minimizers.push(min_pos);
            start = min_pos + 1;
        }
        // Possibly add one last minimizer.
        let start = hashes.len() - self.w;
        let min_pos = start + hashes[start..].iter().position_min().expect("w > 0");
        if minimizers.last() != Some(&min_pos) {
            minimizers.push(min_pos);
        }
        minimizers.iter().map(|e| (hashes[*e], *e)).collect()
    }

    fn mins(&mut self, text: &[u8]) -> Vec<u64> {
        self.minimizer(text).iter().map(|(min, _)| *min).collect()
    }
}

pub fn seq_mins(
    seq: &[u8],
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
) -> anyhow::Result<rapidhash::RapidHashSet<u64>> {
    let minimizers: Vec<u64> = match opt_hasher {
        "rapid" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: RapidHash,
        }
        .mins(&seq[..]),
        "fx" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: FxHash,
        }
        .mins(&seq[..]),
        "murmur" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: MurmurHash3,
        }
        .mins(&seq[..]),
        "mod" => {
            let min_iter = minimizer_iter::MinimizerBuilder::<u64, _>::new_mod()
                .canonical()
                .minimizer_size(opt_kmer)
                .width(opt_window as u16)
                .iter(&seq[..]);

            min_iter.map(|(min, _, _)| min).collect()
        }
        _ => unreachable!(),
    };
    let hashset: rapidhash::RapidHashSet<u64> = rapidhash::RapidHashSet::from_iter(minimizers);

    Ok(hashset)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MinimizerInfo {
    pub hash: u64,
    pub seq_id: u32,
    pub pos: u32,
    pub strand: bool, // true: +, false: -
}

// Wrapper for Filter logic
struct FilterBuildHasher<'a, F> {
    filter: &'a F,
}

impl<'a, F> Clone for FilterBuildHasher<'a, F> {
    fn clone(&self) -> Self {
        FilterBuildHasher {
            filter: self.filter,
        }
    }
}

impl<'a, F> Copy for FilterBuildHasher<'a, F> {}

impl<'a, F> std::hash::BuildHasher for FilterBuildHasher<'a, F>
where
    F: Fn(u64) -> bool,
{
    type Hasher = FilterHasher<'a, F>;
    fn build_hasher(&self) -> Self::Hasher {
        FilterHasher {
            filter: self.filter,
            state: 0,
        }
    }
}

struct FilterHasher<'a, F> {
    filter: &'a F,
    state: u64,
}

impl<'a, F> std::hash::Hasher for FilterHasher<'a, F>
where
    F: Fn(u64) -> bool,
{
    fn write(&mut self, bytes: &[u8]) {
        // Use RapidHash logic directly
        self.state = rapidhash::rapidhash(bytes);
    }

    fn finish(&self) -> u64 {
        let h = self.state;
        // If filter returns false, we want to reject this hash.
        // minimizer_iter selects the MINIMUM hash.
        // If we return u64::MAX, it will be ignored unless all hashes in window are MAX.
        if (self.filter)(h) {
            h
        } else {
            u64::MAX
        }
    }
}

/// Sketch a sequence to find minimizers.
///
/// # Arguments
/// * `seq` - The DNA sequence
/// * `seq_id` - ID of the sequence
/// * `k` - K-mer size
/// * `w` - Window size
/// * `soft_mask` - If true, ignore k-mers containing lowercase bases
/// * `filter` - A predicate that returns true if a hash should be KEPT.
pub fn seq_sketch<F>(
    seq: &[u8],
    seq_id: u32,
    k: usize,
    w: usize,
    soft_mask: bool,
    filter: F,
) -> Vec<MinimizerInfo>
where
    F: Fn(u64) -> bool,
{
    // Use minimizer_iter with our custom FilterBuildHasher
    let build_hasher = FilterBuildHasher { filter: &filter };

    let builder = MinimizerBuilder::<u64, _>::new()
        .minimizer_size(k)
        .width(w as u16)
        .canonical() // Ensure canonical minimizers (min of fwd/rev)
        .hasher(build_hasher);

    // If soft_mask is enabled, we filter out minimizers that overlap with lowercase regions.
    // We check the original sequence at the minimizer's position.
    // Note: This effectively drops windows where the minimizer falls in a masked region.
    // It avoids allocation and complex iterator mapping.

    builder
        .iter(seq)
        .map(|(hash, pos, is_rc)| {
            let strand = !is_rc;
            MinimizerInfo {
                hash,
                seq_id,
                pos: pos as u32,
                strand,
            }
        })
        .filter(|m| {
            if m.hash == u64::MAX {
                return false;
            } // Should be filtered by FilterHasher if used, but explicit check is fine.

            if soft_mask {
                let start = m.pos as usize;
                let end = start + k;
                if end > seq.len() {
                    return false;
                } // Should not happen

                // Check if any byte in seq[start..end] is lowercase
                !seq[start..end].iter().any(|&b| b.is_ascii_lowercase())
            } else {
                true
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_sketch_basic() {
        let seq = b"ACGTACGT";
        let k = 3;
        let w = 3; // minimizer_iter requires odd window size?
        let mins = seq_sketch(seq, 1, k, w, false, |_| true);

        assert!(!mins.is_empty());
        for m in &mins {
            assert_eq!(m.seq_id, 1);
            assert!(m.pos < seq.len() as u32);
        }
    }

    #[test]
    fn test_seq_sketch_soft_mask() {
        // "acgt" is lowercase. If soft_mask is true, it should be ignored.
        let seq = b"ACGTacgtACGT";
        let k = 4;
        let w = 1; // Small window to test individual k-mers

        // Without soft mask: "acgt" (lowercase) is a valid k-mer (different hash from ACGT)
        let mins_no_mask = seq_sketch(seq, 1, k, w, false, |_| true);
        // ACGT (0), CGTa (1), GTac (2), Tacg (3), acgt (4), cgtA (5), gtAC (6), tACG (7), ACGT (8)
        // We expect some minimizers.
        assert!(mins_no_mask.len() >= 2);

        // With soft mask: "acgt" and any k-mer containing lowercase should be ignored.
        // k-mers containing lowercase:
        // CGTa, GTac, Tacg, acgt, cgtA, gtAC, tACG
        // Only ACGT (0) and ACGT (8) are purely uppercase.
        let mins_mask = seq_sketch(seq, 1, k, w, true, |_| true);

        // Should only find the two uppercase blocks
        // ACGT at 0
        // ACGT at 8
        // Note: minimizer_iter with w=1 returns all valid k-mers.
        // But if we return u64::MAX, they are filtered out.
        assert_eq!(mins_mask.len(), 2);
        assert_eq!(mins_mask[0].pos, 0);
        assert_eq!(mins_mask[1].pos, 8);
    }

    #[test]
    fn test_seq_sketch_strand() {
        // AAAA (fwd) vs TTTT (rev)
        // If canonical is working, it should pick the smaller hash.
        // Let's rely on consistency.
        let seq = b"ACGT";
        let k = 4;
        let w = 1;
        let mins = seq_sketch(seq, 1, k, w, false, |_| true);
        assert_eq!(mins.len(), 1);
    }
}
