#![feature(array_chunks)]
#![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]
use std::simd::prelude::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::rngs::SmallRng;
use rand::rngs::StdRng;
use rand::{Rng, RngCore, SeedableRng};
use rapidhash::RapidRng;
use std::collections::HashSet;
use wyhash::WyRng;

pub fn encode_hash_hd_simd(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i16> {
    let num_seed = kmer_hash_set.len();
    let mut hv = vec![-(num_seed as i16); hv_d];

    let num_chunk = hv_d / 64;

    // Convert HashSet to Vec
    let seed_vec: Vec<u64> = kmer_hash_set.iter().cloned().collect();

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(hash);

        // SIMD-based HV encoding
        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();

            // Use SIMD to process 4 bits at a time
            for j in (0..64).step_by(4) {
                // Create a SIMD vector of 4 bits
                let bit_mask = u64x4::splat(1);
                let shift =
                    Simd::from_array([j as u64, (j + 1) as u64, (j + 2) as u64, (j + 3) as u64]);
                let bits = (u64x4::splat(rnd_bits) >> shift) & bit_mask;

                // Convert bits to i16 and shift left by 1
                let bits_i16 = bits.cast::<i16>() << Simd::splat(1);

                // Load the target HV values
                let mut hv_simd = i16x4::from_slice(&hv[i * 64 + j..i * 64 + j + 4]);

                // Accumulate the bits
                hv_simd += bits_i16;

                // Store the updated HV values
                hv_simd.copy_to_slice(&mut hv[i * 64 + j..i * 64 + j + 4]);
            }
        }
    }

    hv
}

pub fn encode_hash_hd_simd2(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i32> {
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

pub fn encode_hash_hd_rapid(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i16> {
    let seed_vec = Vec::from_iter(kmer_hash_set.clone());
    let mut hv = vec![-(kmer_hash_set.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_wy(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i16> {
    let seed_vec = Vec::from_iter(kmer_hash_set.clone());
    let mut hv = vec![-(kmer_hash_set.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = WyRng::seed_from_u64(hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_std(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i16> {
    let seed_vec = Vec::from_iter(kmer_hash_set.clone());
    let mut hv = vec![-(kmer_hash_set.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = StdRng::seed_from_u64(hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_small(kmer_hash_set: &HashSet<u64>, hv_d: usize) -> Vec<i16> {
    let seed_vec = Vec::from_iter(kmer_hash_set.clone());
    let mut hv = vec![-(kmer_hash_set.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = SmallRng::seed_from_u64(hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

// Generate a random k-mer hash set
fn generate_kmer_hash_set(size: usize) -> HashSet<u64> {
    let mut rng = StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut kmer_hash_set = HashSet::with_capacity(size);

    for _ in 0..size {
        kmer_hash_set.insert(rng.gen::<u64>());
    }

    kmer_hash_set
}

// Benchmark function
fn bench_encode_hash_hd(c: &mut Criterion) {
    // Create test datasets of different sizes
    let kmer_hash_set_small = generate_kmer_hash_set(1000); // Small dataset
    let kmer_hash_set_medium = generate_kmer_hash_set(10_000); // Medium dataset

    let hv_d = 4096; // Hypervector dimension

    // Benchmark small dataset
    c.bench_function("encode_hash_hd_simd_small", |b| {
        b.iter(|| encode_hash_hd_simd(black_box(&kmer_hash_set_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd2_small", |b| {
        b.iter(|| encode_hash_hd_simd2(black_box(&kmer_hash_set_small), hv_d))
    });
    c.bench_function("encode_hash_hd_rapid_small", |b| {
        b.iter(|| encode_hash_hd_rapid(black_box(&kmer_hash_set_small), hv_d))
    });
    c.bench_function("encode_hash_hd_wy_small", |b| {
        b.iter(|| encode_hash_hd_wy(black_box(&kmer_hash_set_small), hv_d))
    });
    c.bench_function("encode_hash_hd_std_small", |b| {
        b.iter(|| encode_hash_hd_std(black_box(&kmer_hash_set_small), hv_d))
    });
    c.bench_function("encode_hash_hd_small_small", |b| {
        b.iter(|| encode_hash_hd_small(black_box(&kmer_hash_set_small), hv_d))
    });

    // // Benchmark medium dataset
    // c.bench_function("encode_hash_hd_rapid_medium", |b| {
    //     b.iter(|| encode_hash_hd_rapid(black_box(&kmer_hash_set_medium), hv_d))
    // });
    // c.bench_function("encode_hash_hd_wy_medium", |b| {
    //     b.iter(|| encode_hash_hd_wy(black_box(&kmer_hash_set_medium), hv_d))
    // });
    // c.bench_function("encode_hash_hd_std_medium", |b| {
    //     b.iter(|| encode_hash_hd_std(black_box(&kmer_hash_set_medium), hv_d))
    // });
    // c.bench_function("encode_hash_hd_small_medium", |b| {
    //     b.iter(|| encode_hash_hd_small(black_box(&kmer_hash_set_medium), hv_d))
    // });
}

// Define benchmark group
criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10); // Set sample size
    targets = bench_encode_hash_hd
);

// Run benchmarks
criterion_main!(benches);
