// #![feature(array_chunks)]
// #![feature(slice_as_chunks)]
// // Add these imports to use the stdsimd library
// #![feature(portable_simd)]
// use std::simd::prelude::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::Rng;
use std::collections::{BTreeSet, HashSet};

fn generate_random_set(len: usize, size: usize) -> Vec<Vec<u64>> {
    let mut rng = rand::thread_rng();
    let side = rand::distributions::Uniform::new(0, u64::MAX);

    let mut vec = Vec::new();
    for _ in 0..size {
        let set: Vec<u64> = (0..len).map(|_| rng.sample(side)).collect();
        vec.push(set)
    }

    vec
}

const LEN: usize = 1005;
const SIZE: usize = 105;
fn btree_inter(c: &mut Criterion) {
    let vec : Vec<BTreeSet<u64>>= generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_inter", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let inter: BTreeSet<_> = set1.intersection(set2).cloned().collect();
            black_box(inter);
        })
    });
}

fn btree_union(c: &mut Criterion) {
    let vec : Vec<BTreeSet<u64>>= generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_union", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let union: BTreeSet<_> = set1.union(set2).cloned().collect();
            black_box(union);
        })
    });
}

fn hashset_inter(c: &mut Criterion) {
    let vec : Vec<HashSet<u64>>= generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_inter", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let inter: HashSet<_> = set1.intersection(set2).cloned().collect();
            black_box(inter);
        })
    });
}

fn hashset_union(c: &mut Criterion) {
    let vec : Vec<HashSet<u64>>= generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_union", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let union: HashSet<_> = set1.union(set2).cloned().collect();
            black_box(union);
        })
    });
}

criterion_group!(benches, btree_inter, btree_union, hashset_inter, hashset_union);
criterion_main!(benches);

// msvc
// norm_map                time:   [8.9184 µs 8.9294 µs 8.9423 µs]
// norm_fold               time:   [8.9407 µs 8.9542 µs 8.9664 µs]
// norm_simd               time:   [1.1196 µs 1.1214 µs 1.1230 µs]

// wsl2
// norm_map                time:   [8.9086 µs 8.9226 µs 8.9381 µs]
// norm_fold               time:   [8.9320 µs 8.9456 µs 8.9621 µs]
// norm_simd               time:   [1.1420 µs 1.1442 µs 1.1470 µs]
