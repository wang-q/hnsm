use std::simd::prelude::*;

use rand::{RngCore, SeedableRng};
use rapidhash::{RapidHashSet, RapidRng};

#[allow(dead_code)]
// The original implementation is i16
fn hash_hv_serial(kmer_hash_set: &RapidHashSet<u64>, hv_d: usize) -> Vec<i32> {
    // Convert the k-mer hash set to a Vec for faster iteration
    let seed_vec = Vec::from_iter(kmer_hash_set.clone());

    // Initialize the hypervector (HV)
    // The initial value of HV is -(kmer_hash_set.len() as i32), which is the negative of the k-mer hash set size,
    //   corresponding to the `-1` part of the equation
    // hv_d is the dimension of the hypervector
    let mut hv = vec![-(kmer_hash_set.len() as i32); hv_d];

    // Iterate over each k-mer hash value
    for hash in seed_vec {
        // Initialize the random number generator (RapidRng) using the k-mer hash as the seed
        let mut rng = RapidRng::seed_from_u64(hash);

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
/// This function uses SIMD instructions to process 4 bits at a time, improving performance over the serial implementation.
pub fn hash_hv(kmer_hash_set: &RapidHashSet<u64>, hv_d: usize) -> Vec<i32> {
    let num_seed = kmer_hash_set.len();
    let mut hv = vec![-(num_seed as i32); hv_d];

    let num_chunk = hv_d / 32;

    // Convert HashSet to Vec
    let seed_vec: Vec<u64> = kmer_hash_set.iter().cloned().collect();

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(hash);

        // SIMD-based HV encoding
        for i in 0..num_chunk {
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

#[allow(dead_code)]
fn hv_norm_l2_sq_serial(hv: &Vec<i32>) -> f32 {
    let norm_sq = hv
        .iter()
        .fold(0., |sum: f32, &num| sum + (num as f32 * num as f32));
    norm_sq
}

/// Computes the squared L2 norm of a hypervector using a SIMD-optimized implementation.
///
/// # Arguments
/// * `a` - The hypervector represented as a slice of `i32`.
///
/// # Returns
/// The squared L2 norm of the hypervector as an `f32`.
pub fn hv_norm_l2_sq(a: &[i32]) -> f32 {
    let a_f32: Vec<f32> = a.iter().map(|&x| x as f32).collect();
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
pub fn hv_cardinality(hv: &[i32], hv_d: usize) -> usize {
    let norm_sq = hv_norm_l2_sq(hv);
    (norm_sq / hv_d as f32) as usize
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
    let a_f32: Vec<_> = a.iter().map(|&x| x as f32).collect();
    let b_f32: Vec<_> = b.iter().map(|&x| x as f32).collect();

    crate::dot_product(&a_f32, &b_f32)
}

// pub fn compute_pairwise_ani(
//     r: &Vec<i16>,
//     norm2_r: i32,
//     q: &Vec<i16>,
//     norm2_q: i32,
//     ksize: u8,
// ) -> f32 {
//     // Scalar-based inner product
//     let dot_r_q: i32 = r
//         .iter()
//         .zip(q.iter())
//         .map(|(x, y)| (*x as i32) * (*y as i32))
//         .sum();
//
//     let jaccard: f32 = dot_r_q as f32 / (norm2_r + norm2_q - dot_r_q) as f32;
//     let ani: f32 = 1.0 + (2.0 / (1.0 / jaccard + 1.0)).ln() / (ksize as f32);
//
//     if ani.is_nan() {
//         0.0
//     } else {
//         ani.min(1.0).max(0.0) * 100.0
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_hash_hv() {
        // Generate random input data
        let mut rng = rand::thread_rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.gen::<u64>()).collect();
        let hv_d = 4096;

        // Run the SIMD version
        let hv = hash_hv(&kmer_hash_set, hv_d);

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
        let mut rng = rand::thread_rng();
        let kmer_hash_set: RapidHashSet<u64> = (0..1000).map(|_| rng.gen::<u64>()).collect();
        let hv_d = 4096;

        // Run normal version
        let result_serial = hash_hv_serial(&kmer_hash_set, hv_d);

        // Run SIMD version
        let result_simd = hash_hv(&kmer_hash_set, hv_d);

        println!(
            "Size of kmer_hash_set: {} bytes",
            kmer_hash_set.capacity() * size_of::<u64>()
        );
        println!(
            "Size of result_rapid: {} bytes",
            result_serial.capacity() * size_of::<i16>()
        );

        let cardinality = hv_cardinality(&result_serial, hv_d);
        println!("HV cardinality: {}", cardinality);

        // Compare results
        assert_eq!(
            result_serial, result_simd,
            "SIMD version does not match normal version!"
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
        let hv_d = 5;

        // Compute the cardinality
        let cardinality = hv_cardinality(&hv, hv_d);

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

        // Expected result: 1*2 + 2*3 + 3*4 + 4*5 + 5*6 = 2 + 6 + 12 + 20 + 30 = 70
        assert_eq!(dot, 70.0, "Dot product calculation is incorrect!");
    }

    #[test]
    fn test_hv_dot_orthogonal() {
        // Create two orthogonal hypervectors
        let a = vec![1, 0, 0];
        let b = vec![0, 1, 0];

        // Compute the dot product
        let dot = hv_dot(&a, &b);

        // Expected result: 1*0 + 0*1 + 0*0 = 0
        assert_eq!(
            dot, 0.0,
            "Dot product of orthogonal vectors should be zero!"
        );
    }

    #[test]
    fn test_hv_cardinality_zero() {
        // Create a hypervector with all zeros
        let hv = vec![0, 0, 0, 0, 0];
        let hv_d = 5;

        // Compute the cardinality
        let cardinality = hv_cardinality(&hv, hv_d);

        // Expected result: 0
        assert_eq!(
            cardinality, 0,
            "Cardinality of a zero vector should be zero!"
        );
    }
}
