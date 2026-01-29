use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

pub struct BloomFilter {
    bit_vec: Vec<u64>,
    num_bits: u64,
    num_hashes: u32,
}

impl BloomFilter {
    pub fn new(expected_elements: usize, false_positive_rate: f64) -> Self {
        let num_bits =
            (-(expected_elements as f64) * false_positive_rate.ln() / (2.0f64.ln().powi(2))) as u64;
        let num_hashes = ((num_bits as f64 / expected_elements as f64) * 2.0f64.ln()) as u32;

        let num_u64 = (num_bits + 63) / 64;
        let bit_vec = vec![0; num_u64 as usize];

        Self {
            bit_vec,
            num_bits,
            num_hashes,
        }
    }

    fn get_hashes(&self, item: u64) -> (u64, u64) {
        // Use double hashing with the item itself as the source of entropy
        // item is already a hash (minimizer hash)
        // We can use it directly or hash it again
        let mut hasher1 = DefaultHasher::new();
        item.hash(&mut hasher1);
        let h1 = hasher1.finish();

        let mut hasher2 = DefaultHasher::new();
        item.hash(&mut hasher2);
        // Rotate to get a different hash
        let h2 = hasher2.finish().rotate_left(32);

        (h1, h2)
    }

    pub fn insert(&mut self, item: u64) {
        let (h1, h2) = self.get_hashes(item);
        for i in 0..self.num_hashes {
            let hash = h1.wrapping_add((i as u64).wrapping_mul(h2));
            let bit_idx = hash % self.num_bits;
            let vec_idx = (bit_idx / 64) as usize;
            let bit_offset = (bit_idx % 64) as usize;
            self.bit_vec[vec_idx] |= 1 << bit_offset;
        }
    }

    pub fn contains(&self, item: u64) -> bool {
        let (h1, h2) = self.get_hashes(item);
        for i in 0..self.num_hashes {
            let hash = h1.wrapping_add((i as u64).wrapping_mul(h2));
            let bit_idx = hash % self.num_bits;
            let vec_idx = (bit_idx / 64) as usize;
            let bit_offset = (bit_idx % 64) as usize;
            if (self.bit_vec[vec_idx] & (1 << bit_offset)) == 0 {
                return false;
            }
        }
        true
    }
}
