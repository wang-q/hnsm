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
    let side = rand::distributions::Uniform::new(0, u16::MAX as u64);

    let mut vec = Vec::new();
    for _ in 0..size {
        let set: Vec<u64> = (0..len).map(|_| rng.sample(side)).collect();
        vec.push(set)
    }

    vec
}

const LEN: usize = 5005;
const SIZE: usize = 105;

fn btree_inter(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_inter", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let inter = set1.intersection(set2).cloned().count();
            black_box(inter);
        })
    });
}

fn btree_union(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_union", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let union = set1.union(set2).cloned().count();
            black_box(union);
        })
    });
}

fn hashset_inter(c: &mut Criterion) {
    let vec: Vec<HashSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_inter", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let inter = set1.intersection(set2).count();
            black_box(inter);
        })
    });
}

fn hashset_union(c: &mut Criterion) {
    let vec: Vec<HashSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_union", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let union = set1.union(set2).count();
            black_box(union);
        })
    });
}

fn btree_jaccard(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).cloned().count();
            let union = set1.union(set2).cloned().count();

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn btree_jaccard2(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_jaccard2", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).cloned().count();
            let union = set1.len() + set2.len() + inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn hashset_jaccard2(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_jaccard2", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).count();
            let union = set1.len() + set2.len() + inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

criterion_group!(
    benches,
    btree_inter,
    btree_union,
    hashset_inter,
    hashset_union,
    btree_jaccard,
    btree_jaccard2,
    hashset_jaccard2
);
criterion_main!(benches);

// msvc
// btree_jaccard2          time:   [44.467 µs 45.254 µs 46.136 µs]
// hashset_jaccard2        time:   [48.958 µs 50.457 µs 52.138 µs]
