use rand::{RngCore, SeedableRng};
use rapidhash::RapidRng;
use std::simd::prelude::*;

#[allow(dead_code)]
// The original implementation is i16
fn hash_hv_serial(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Initialize the hypervector (HV)
    // The initial value of HV is -(seed_vec.len() as i32), which is the negative of the k-mer hash set size,
    //   corresponding to the `-1` part of the equation
    // hv_d is the dimension of the hypervector
    let mut hv = vec![-(seed_vec.len() as i32); hv_d];

    // Iterate over each k-mer hash value
    for hash in seed_vec {
        // Initialize the random number generator (RapidRng) using the k-mer hash as the seed
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Generate a random stream of bits of sufficient length
        for i in 0..(hv_d / 32) {
            // Generate a 64-bit random number
            let rnd_bits = rng.next_u32();

            // Iterate over each bit in the random number
            // Formula: hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i32
            // - (rnd_bits >> j) & 1: Extract the j-th bit (0 or 1)
            // - << 1: Shift the bit left by 1 (resulting in 0 or 2, i.e. multiply by 2)
            // - Add the result to the hypervector at position i * 64 + j
            for j in 0..32 {
                hv[i * 32 + j] += (((rnd_bits >> j) & 1) << 1) as i32;
            }
        }
    }

    hv
}

/// Generates a hypervector (HV) from a set of k-mer hash values using a SIMD-optimized implementation.
///
/// # Arguments
/// * `kmer_hash_set` - A set of k-mer hash values.
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// A hypervector of dimension `hv_d` represented as a `Vec<i32>`.
///
/// # Formula
/// The hypervector is generated as:
/// \[
/// \mathbf{H} = \sum_{i=1}^{N} (hv^{i} \times 2 - 1)
/// \]
/// where \(N\) is the number of k-mer hash values, and \(hv^{i}\) is a binary hypervector derived from the k-mer hash.
///
/// # Notes
/// This function uses SIMD instructions to process 8 bits at a time, improving performance over the serial implementation.
pub fn hash_hv_bit(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    let num_seed = seed_vec.len();
    let num_chunk = hv_d / 32;
    let mut hv = vec![-(num_seed as i32); hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // SIMD-based HV encoding
        for i in 0..num_chunk {
            // 32 * 8 can be fit into an AVX2 register
            let rnd_bits = rng.next_u32();

            // Use SIMD to process 8 bits at a time
            for j in (0..32).step_by(8) {
                let bit_mask = u32x8::splat(1);
                let shift = Simd::from_array([
                    j as u32,
                    (j + 1) as u32,
                    (j + 2) as u32,
                    (j + 3) as u32,
                    (j + 4) as u32,
                    (j + 5) as u32,
                    (j + 6) as u32,
                    (j + 7) as u32,
                ]);
                let bits = (u32x8::splat(rnd_bits) >> shift) & bit_mask;

                // Convert bits to i32 and shift left by 1
                let bits_i32 = bits.cast::<i32>() << Simd::splat(1);

                // Load the target HV values
                let mut hv_simd = i32x8::from_slice(&hv[i * 32 + j..i * 32 + j + 8]);

                // Accumulate the bits
                hv_simd += bits_i32;

                // Store the updated HV values
                hv_simd.copy_to_slice(&mut hv[i * 32 + j..i * 32 + j + 8]);
            }
        }
    }

    hv
}

/// Generates a hypervector (HV) using i8 values (-128..=127) as the basic unit.
///
/// # Arguments
/// * `seed_vec` - A set of k-mer hash values (seeds).
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// A hypervector of dimension `hv_d` represented as a `Vec<i32>`.
///
/// # Notes
/// This implementation avoids bit manipulation overhead by using `i8` directly,
/// but requires more RNG calls (1 u64 per 8 dimensions) compared to the bit-based approach.
/// It uses SIMD to process 8 dimensions at a time.
pub fn hash_hv_i8(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Initialize HV with 0.
    // We accumulate random i8 values (-128..=127) directly.
    let mut hv = vec![0i32; hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Process 8 dimensions per chunk (1 u64 = 8 bytes)
        let num_chunk = hv_d / 8;

        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();

            // Load 8 bytes into SIMD vector
            let vec_u8 = u8x8::from_array(bytes);

            // Cast u8 to i8 (0..255 -> 0..127, -128..-1)
            // Then cast to i32 for accumulation
            let vec_vals = vec_u8.cast::<i8>().cast::<i32>();

            // Load current HV values
            let mut hv_simd = i32x8::from_slice(&hv[i * 8..(i + 1) * 8]);

            // Accumulate
            hv_simd += vec_vals;

            // Store back
            hv_simd.copy_to_slice(&mut hv[i * 8..(i + 1) * 8]);
        }
    }

    hv
}

#[allow(dead_code)]
fn hash_hv_i8_serial(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    // Initialize HV with 0.
    let mut hv = vec![0i32; hv_d];

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // Process dimensions in chunks of 8 (matching SIMD behavior for RNG alignment)
        let num_chunk = hv_d / 8;

        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let bytes = rnd_bits.to_ne_bytes();

            // Iterate over each byte in the u64
            for j in 0..8 {
                let val_i8 = bytes[j] as i8;
                hv[i * 8 + j] += val_i8 as i32;
            }
        }
    }

    hv
}

#[allow(dead_code)]
fn hv_norm_l2_sq_serial(hv: &[i32]) -> f32 {
    let norm_sq = hv
        .iter()
        .fold(0., |sum: f32, &num| sum + (num as f32 * num as f32));
    norm_sq
}

/// Computes the squared L2 norm of a hypervector using a SIMD-optimized implementation.
///
/// # Arguments
/// * `hv` - The hypervector represented as a slice of `i32`.
///
/// # Returns
/// The squared L2 norm of the hypervector as an `f32`.
pub fn hv_norm_l2_sq(hv: &[i32]) -> f32 {
    let a_f32: Vec<f32> = hv.iter().map(|&x| x as f32).collect();
    crate::norm_l2_sq(&a_f32)
}

/// Computes the cardinality of a set represented by a hypervector.
///
/// # Arguments
/// * `hv` - The hypervector represented as a slice of `i32`.
/// * `hv_d` - The dimension of the hypervector.
///
/// # Returns
/// The cardinality of the set as a `usize`.
///
/// # Formula
/// The cardinality is computed as:
/// \[
/// |\mathcal{S}_k(A)| = \frac{\|\mathbf{H}_A\|_2^2}{D}
/// \]
/// where \(\|\mathbf{H}_A\|_2^2\) is the squared L2 norm of the hypervector, and \(D\) is the dimension of the hypervector.
pub fn hv_cardinality(hv: &[i32]) -> usize {
    let norm_sq = hv_norm_l2_sq(hv);
    (norm_sq / hv.len() as f32) as usize
}

/// Computes the dot product of two hypervectors.
///
/// # Arguments
/// * `a` - The first hypervector represented as a slice of `i32`.
/// * `b` - The second hypervector represented as a slice of `i32`.
///
/// # Returns
/// The dot product of the two hypervectors as an `f32`.
pub fn hv_dot(a: &[i32], b: &[i32]) -> f32 {
    let hv_d_sqrt = (a.len() as f32).sqrt();
    let a_f32: Vec<_> = a.iter().map(|&x| x as f32 / hv_d_sqrt).collect();
    let b_f32: Vec<_> = b.iter().map(|&x| x as f32 / hv_d_sqrt).collect();

    crate::dot_product(&a_f32, &b_f32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rapidhash::RapidHashSet;

    #[test]
    fn test_hash_hv() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run the SIMD version
        let hv = hash_hv_bit(&seed_vec, hv_d);

        // Check the dimension of the hypervector
        assert_eq!(hv.len(), hv_d, "Hypervector dimension mismatch!");

        // Check that the hypervector is not all zeros
        assert!(
            hv.iter().any(|&x| x != 0),
            "Hypervector should not be all zeros!"
        );
    }

    #[test]
    fn test_hash_hv_serial_vs_simd() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run normal version
        let result_serial = hash_hv_serial(&seed_vec, hv_d);

        // Run SIMD version
        let result_simd = hash_hv_bit(&seed_vec, hv_d);

        println!(
            "Size of seed_vec: {} bytes",
            seed_vec.capacity() * size_of::<u64>()
        );
        println!(
            "Size of result_rapid: {} bytes",
            result_serial.capacity() * size_of::<i16>()
        );

        let cardinality = hv_cardinality(&result_serial);
        println!("HV cardinality: {}", cardinality);

        // Compare results
        assert_eq!(
            result_serial, result_simd,
            "SIMD version does not match normal version!"
        );
    }

    #[test]
    fn test_hash_hv_i8() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run the i8 SIMD version
        let hv = hash_hv_i8(&seed_vec, hv_d);

        // Check the dimension of the hypervector
        assert_eq!(hv.len(), hv_d, "Hypervector dimension mismatch!");

        // Check that the hypervector is not all zeros
        assert!(
            hv.iter().any(|&x| x != 0),
            "Hypervector should not be all zeros!"
        );
    }

    #[test]
    fn test_hash_hv_i8_serial_vs_simd() {
        // Generate random input data
        let mut rng = rand::rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.random::<u64>()).collect();
        let seed_vec: Vec<u64> = kmer_hash_set.into_iter().collect();
        let hv_d = 4096;

        // Run normal version
        let result_serial = hash_hv_i8_serial(&seed_vec, hv_d);

        // Run SIMD version
        let result_simd = hash_hv_i8(&seed_vec, hv_d);

        // Compare results
        assert_eq!(
            result_serial, result_simd,
            "SIMD version does not match serial version for i8 implementation!"
        );
    }

    #[test]
    fn test_hv_norm_l2_sq() {
        // Create a simple hypervector
        let hv = vec![1, 2, 3, 4, 5];

        // Compute the squared L2 norm
        let norm_sq = hv_norm_l2_sq(&hv);

        // Expected result: 1^2 + 2^2 + 3^2 + 4^2 + 5^2 = 55
        assert_eq!(norm_sq, 55.0, "Squared L2 norm calculation is incorrect!");
    }

    #[test]
    fn test_hv_norm_l2_sq_serial_vs_simd() {
        let hv: Vec<_> = (1..=32).map(|x| x as i32).collect();

        let result_scalar = hv_norm_l2_sq_serial(&hv);
        let result_simd = hv_norm_l2_sq(&hv);

        println!("Scalar result: {}", result_scalar);
        println!("SIMD result: {}", result_simd);

        assert_eq!(result_scalar, result_simd, "Results do not match!");
    }

    #[test]
    fn test_hv_cardinality() {
        // Create a simple hypervector
        let hv = vec![1, 2, 3, 4, 5];

        // Compute the cardinality
        let cardinality = hv_cardinality(&hv);

        // Expected result: (1^2 + 2^2 + 3^2 + 4^2 + 5^2) / 5 = 55 / 5 = 11
        assert_eq!(cardinality, 11, "Cardinality calculation is incorrect!");
    }

    #[test]
    fn test_hv_dot() {
        // Create two simple hypervectors
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![2, 3, 4, 5, 6];

        // Compute the dot product
        let dot = hv_dot(&a, &b);

        // Expected result: (1*2 + 2*3 + 3*4 + 4*5 + 5*6) / 5 = 14
        assert_eq!(dot, 14.0, "Dot product calculation is incorrect!");
    }

    #[test]
    fn test_hv_dot_orthogonal() {
        // Create two orthogonal hypervectors
        let a = vec![1, 0, 0];
        let b = vec![0, 1, 0];

        // Compute the dot product
        let dot = hv_dot(&a, &b);

        // Expected result: (1*0 + 0*1 + 0*0) / 3 = 0
        assert_eq!(
            dot, 0.0,
            "Dot product of orthogonal vectors should be zero!"
        );
    }
}
