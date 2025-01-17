use std::simd::prelude::*;

/// Number of lanes in the SIMD vector.
/// Each SIMD vector can process 8 `f32` elements at once.
/// 32 * 8 = 256, AVX2
const LANES: usize = 8;

// https://www.maartengrootendorst.com/blog/distances/
// https://crates.io/crates/semanticsimilarity_rs

/// Computes the Euclidean distance between two vectors `a` and `b`.
///
/// # Arguments
/// * `a` - The first vector.
/// * `b` - The second vector.
///
/// # Returns
/// The Euclidean distance between `a` and `b`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let b = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
/// let distance = hnsm::euclidean_distance(&a, &b);
/// assert_eq!(format!("{:.4}", distance), "18.1659".to_string());
/// ```
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

/// Computes the dot product of two vectors `a` and `b`.
///
/// # Arguments
/// * `a` - The first vector.
/// * `b` - The second vector.
///
/// # Returns
/// The dot product of `a` and `b`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let b = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
/// let dot = hnsm::dot_product(&a, &b);
/// assert_eq!(dot, 220.0);
/// ```
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

/// Computes the L2 norm (Euclidean norm) of a vector `a`.
///
/// # Arguments
/// * `a` - The vector.
///
/// # Returns
/// The L2 norm of `a`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let norm = hnsm::norm_l2(&a);
/// assert_eq!(format!("{:.4}", norm), "19.6214".to_string());
/// ```
#[inline]
pub fn norm_l2(a: &[f32]) -> f32 {
    norm_l2_sq(a).sqrt()
}

/// Computes the squared L2 norm of a vector `a`.
///
/// # Arguments
/// * `a` - The vector.
///
/// # Returns
/// The squared L2 norm of `a`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let norm_sq = hnsm::norm_l2_sq(&a);
/// assert_eq!(norm_sq, 385.0);
/// ```
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

/// Computes the Jaccard intersection of two vectors `a` and `b`.
/// The Jaccard intersection is the sum of the minimum values of corresponding elements.
///
/// # Arguments
/// * `a` - The first vector.
/// * `b` - The second vector.
///
/// # Returns
/// The Jaccard intersection of `a` and `b`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let b = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
/// let intersection = hnsm::jaccard_intersection(&a, &b);
/// assert_eq!(intersection, 30.0);
/// ```
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

/// Computes the Jaccard union of two vectors `a` and `b`.
/// The Jaccard union is the sum of the maximum values of corresponding elements.
///
/// # Arguments
/// * `a` - The first vector.
/// * `b` - The second vector.
///
/// # Returns
/// The Jaccard union of `a` and `b`.
///
/// # Examples
/// ```
/// let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
/// let b = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
/// let union = hnsm::jaccard_union(&a, &b);
/// assert_eq!(union, 80.0);
/// ```
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
