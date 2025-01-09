#![feature(array_chunks)]
#![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]
use std::simd::prelude::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};

// Generate a random vector of f32 values
fn rand_vec(len: usize) -> Vec<f32> {
    let mut rng = thread_rng();
    let side = Uniform::new(0.0, 10.0);

    (0..len).map(|_| rng.sample(side) as f32).collect()
}

// Calculate the L2 norm using map and sum
fn norm_map(a: &[f32]) -> f32 {
    a.iter().map(|x| x.powi(2)).sum::<f32>().sqrt()
}

// Calculate the L2 norm using fold
fn norm_fold(a: &[f32]) -> f32 {
    a.iter().fold(0.0f32, |s, x| s + x.powi(2)).sqrt()
}

const LANES: usize = 8;
// Calculate the L2 norm using SIMD
fn norm_simd(a: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();

    // Initialize sums with the extra elements
    let mut sums = [0.0; LANES];
    for (x, d) in std::iter::zip(a_extra, &mut sums) {
        *d = x * x;
    }

    // Accumulate sums using SIMD
    let mut sums = f32x8::from_array(sums);
    a_chunks.into_iter().for_each(|x| {
        sums += f32x8::from_array(*x) * f32x8::from_array(*x);
    });

    sums.reduce_sum().sqrt()
}

pub fn bench_rand(c: &mut Criterion) {
    c.bench_function("rand_vec", |b| b.iter(|| rand_vec(black_box(10005))));
}

pub fn bench_norm(c: &mut Criterion) {
    let v1 = rand_vec(10005);

    // Verify that all implementations produce the same result
    assert_eq!(norm_map(&v1), norm_fold(&v1));
    approx::assert_relative_eq!(norm_map(&v1), norm_simd(&v1), epsilon = 0.01);

    // Benchmark each implementation
    c.bench_function("norm_map", |b| b.iter(|| norm_map(black_box(&v1))));
    c.bench_function("norm_fold", |b| b.iter(|| norm_fold(black_box(&v1))));
    c.bench_function("norm_simd", |b| b.iter(|| norm_simd(black_box(&v1))));
}

criterion_group!(benches, bench_rand, bench_norm);
criterion_main!(benches);

// Ryzen 7 8745HS
// msvc
// norm_map                time:   [6.2935 µs 6.3123 µs 6.3328 µs]
// norm_fold               time:   [6.2957 µs 6.3161 µs 6.3383 µs]
// norm_simd               time:   [808.49 ns 810.71 ns 813.05 ns]
// wsl2
// norm_map                time:   [6.3657 µs 6.3785 µs 6.3924 µs]
// norm_fold               time:   [6.3869 µs 6.4031 µs 6.4202 µs]
// norm_simd               time:   [820.59 ns 824.79 ns 831.52 ns]
