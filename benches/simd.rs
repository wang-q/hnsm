#![feature(array_chunks)]
#![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]
use std::simd::prelude::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};

fn rand_vec(len: usize) -> Vec<f32> {
    let mut rng = thread_rng();
    let side = Uniform::new(0.0, 10.0);

    let mut vec: Vec<f32> = Vec::new();
    for _ in 0..len {
        let f = rng.sample(side) as f32;
        vec.push(f)
    }

    vec
}

fn norm_map(a: &[f32]) -> f32 {
    a.iter().map(|x| x.powi(2)).sum::<f32>().sqrt()
}

fn norm_fold(a: &[f32]) -> f32 {
    a.iter().fold(0.0f32, |s, x| s + x.powi(2)).sqrt()
}

const LANES: usize = 8;
fn norm_simd(a: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();

    let mut sums = [0.0; LANES];
    for (x, d) in std::iter::zip(a_extra, &mut sums) {
        *d = x * x;
    }

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
    assert_eq!(norm_map(&v1), norm_fold(&v1));
    approx::assert_relative_eq!(norm_map(&v1), norm_simd(&v1), epsilon = 0.01);

    c.bench_function("norm_map", |b| b.iter(|| norm_map(black_box(&v1))));

    c.bench_function("norm_fold", |b| b.iter(|| norm_fold(black_box(&v1))));

    c.bench_function("norm_simd", |b| b.iter(|| norm_simd(black_box(&v1))));
}

criterion_group!(benches, bench_rand, bench_norm);
criterion_main!(benches);

// msvc
// norm_map                time:   [8.9184 µs 8.9294 µs 8.9423 µs]
// norm_fold               time:   [8.9407 µs 8.9542 µs 8.9664 µs]
// norm_simd               time:   [1.1196 µs 1.1214 µs 1.1230 µs]

// wsl2
// norm_map                time:   [8.9086 µs 8.9226 µs 8.9381 µs]
// norm_fold               time:   [8.9320 µs 8.9456 µs 8.9621 µs]
// norm_simd               time:   [1.1420 µs 1.1442 µs 1.1470 µs]
