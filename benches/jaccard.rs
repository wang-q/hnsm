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

fn btree_access(c: &mut Criterion) {
    let vec: Vec<BTreeSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("btree_access", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();
            let sum = set1.len() + set2.len();

            black_box(sum);
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
            let union = set1.len() + set2.len() - inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn hashset_jaccard(c: &mut Criterion) {
    let vec: Vec<HashSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("hashset_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).count();
            let union = set1.len() + set2.len() - inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn tinyset_jaccard(c: &mut Criterion) {
    let vec: Vec<tinyset::SetU64> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("tinyset_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.iter().filter(|&x| set2.contains(x)).count();
            let union = set1.len() + set2.len() - inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn nohash_jaccard(c: &mut Criterion) {
    let vec: Vec<nohash_hasher::IntSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("nohash_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).count();
            let union = set1.len() + set2.len() - inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn rapidhash_jaccard(c: &mut Criterion) {
    let vec: Vec<rapidhash::RapidHashSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("rapidhash_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).count();
            let union = set1.len() + set2.len() - inter;

            let jaccard = (inter as f64) / (union as f64);
            black_box(jaccard);
        })
    });
}

fn rapidinlinehash_jaccard(c: &mut Criterion) {
    let vec: Vec<rapidhash::RapidInlineHashSet<u64>> = generate_random_set(LEN, SIZE)
        .iter()
        .map(|set| set.iter().cloned().collect())
        .collect();

    c.bench_function("rapidinlinehash_jaccard", |b| {
        b.iter(|| {
            let mut rng = rand::thread_rng();
            let side = rand::distributions::Uniform::new(0, SIZE);
            let set1 = vec.get(rng.sample(side)).unwrap();
            let set2 = vec.get(rng.sample(side)).unwrap();

            let inter = set1.intersection(set2).count();
            let union = set1.len() + set2.len() - inter;

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
    btree_access,
    btree_jaccard,
    hashset_jaccard,
    tinyset_jaccard,
    nohash_jaccard,
    rapidhash_jaccard,
    rapidinlinehash_jaccard,
);
criterion_main!(benches);

// msvc
// btree_access            time:   [8.5894 ns 8.6377 ns 8.6902 ns]
// btree_jaccard           time:   [56.012 µs 57.249 µs 58.617 µs]
// hashset_jaccard         time:   [74.317 µs 75.589 µs 76.957 µs]
// tinyset_jaccard         time:   [53.552 µs 53.821 µs 54.103 µs]
// nohash_jaccard          time:   [60.036 µs 60.417 µs 60.790 µs]
// rapidhash_jaccard       time:   [22.560 µs 22.733 µs 22.913 µs]
// rapidinlinehash_jaccard time:   [22.773 µs 23.040 µs 23.351 µs]
