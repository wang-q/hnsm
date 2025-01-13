use std::simd::prelude::*;

const LANES: usize = 8; // 32 * 8 = 256, AVX2

// https://www.maartengrootendorst.com/blog/distances/
// https://crates.io/crates/semanticsimilarity_rs
pub fn euclidean_distance(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        let diff = x - y;
        *d = diff * diff;
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        let diff = f32x8::from_array(*x) - f32x8::from_array(*y);
        sums += diff * diff;
    });

    sums.reduce_sum().sqrt()
}

pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = x * y;
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::from_array(*x) * f32x8::from_array(*y);
    });

    sums.reduce_sum()
}

#[inline]
pub fn norm_l2(a: &[f32]) -> f32 {
    norm_l2_sq(a).sqrt()
}

pub fn norm_l2_sq(a: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();

    let mut sums = [0.0; LANES];
    for (x, d) in std::iter::zip(a_extra, &mut sums) {
        *d = x * x;
    }

    let mut sums = f32x8::from_array(sums);
    a_chunks.into_iter().for_each(|x| {
        sums += f32x8::from_array(*x) * f32x8::from_array(*x);
    });

    sums.reduce_sum()
}

pub fn jaccard_intersection(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = f32::min(*x, *y);
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::simd_min(f32x8::from_array(*x), f32x8::from_array(*y));
    });

    sums.reduce_sum()
}

pub fn jaccard_union(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = f32::max(*x, *y);
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::simd_max(f32x8::from_array(*x), f32x8::from_array(*y));
    });

    sums.reduce_sum()
}
