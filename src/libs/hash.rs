use itertools::Itertools;

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
