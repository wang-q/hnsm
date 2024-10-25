//! # What is Principal Coordinates Analysis (PCoA)
//!
//! The objective of PCoA is to accurately represent the **distances** between objects using the
//! smallest possible number of dimensions.
//!
//! PCoA with Euclidean distances == classical Multidimensional Scaling (MDS) == PCA
//!
//! Key steps:
//! 1. Distance Matrix. PCoA starts with a distance or dissimilarity matrix that quantifies how
//! similar or different the data points are from each other. Dissimilarity can take various forms
//! as long as it meets the criteria of a metric space. $D = [d_{ij}]$
//!
//! 2. Double-Centering. $B=-{\frac{1}{2}}CD^{(2)}C$ using the centering matrix
//! $C=I-{\frac {1}{n}}J_{n}$ where $n$ is the number of objects, $I$ is the $n \times n$
//! identity matrix, and $J_{n}$ is an $n \times n$ matrix of all ones.
//!
//! 3. Eigenvalue Decomposition. The eigenvalues indicate the amount of variance captured by each
//! dimension, while the eigenvectors provide the coordinates of the data points in the new space.
//!
//! ## Dimensionality Reduction
//!
//! PCoA transforms high-dimensional data into a lower-dimensional space (typically2D or3D)
//! while preserving the relationships between data points, facilitating easier visualization and
//! interpretation.
//!
//!

/*
# R code
matrix <- matrix(c(
    0., 7., 5., 5.,
    7., 0., 4., 9.,
    5., 4., 0., 3.,
    5., 9., 3., 0.
    ), nrow = 4, ncol = 4)

fit <- cmdscale(matrix, k = 2, eig = TRUE, x.ret = TRUE)

fit$x

*/
//
// use ndarray::{Array2, Axis, Zip};
// use ndarray_linalg::{Eigh, UPLO};
//
// pub struct Mds {
//     dim: usize,
// }
//
// impl Mds {
//     pub fn new(dim: usize) -> Self {
//         Mds { dim }
//     }
//
//     pub fn double_centering(&self, matrix: &mut Array2<f32>) {
//         matrix.mapv_inplace(|v| -(v * v) / 2.);
//
//         let row_means = matrix.mean_axis(Axis(0)).unwrap();
//         let col_means = matrix.mean_axis(Axis(1)).unwrap();
//         let matrix_mean = matrix.mean().unwrap();
//
//         Zip::indexed(matrix).for_each(|(row, col), val| {
//             *val = *val - row_means[col] - col_means[row] + matrix_mean
//         });
//     }
//
//     pub fn eigen(&self, matrix: &mut Array2<f32>) {
//         let (eigvals, eigvecs) = matrix.eigh(UPLO::Lower)?;
//
//         eprintln!("eigvals = {:#?}", eigvals);
//
//     }
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn centered() {
//         let mut matrix: Array2<f32> = ndarray::arr2(&[
//             [0., 7., 5., 5.],
//             [7., 0., 4., 9.],
//             [5., 4., 0., 3.],
//             [5., 9., 3., 0.],
//         ]);
//         let exp: Vec<f32> = vec![
//             11.9375, -6.6875, -6.6875, 1.4375, //
//             -6.6875, 23.6875, 3.6875, -20.6875, //
//             -6.6875, 3.6875, -0.3125, 3.3125, //
//             1.4375, -20.6875, 3.3125, 15.9375, //
//         ];
//
//         eprintln!("matrix = {:#?}", matrix);
//         let mut mds = Mds::new(2);
//         mds.double_centering(&mut matrix);
//         eprintln!("matrix = {:#?}", matrix);
//
//         let res = matrix.iter().cloned().collect::<Vec<f32>>();
//
//         assert_eq!(res, exp);
//     }
// }
