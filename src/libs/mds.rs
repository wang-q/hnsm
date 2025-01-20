//! # What is Principal Coordinates Analysis (PCoA)
//!
//! The objective of PCoA is to accurately represent the **distances** between objects using the
//! smallest possible number of dimensions.
//!
//! PCoA with Euclidean distances == classical Multidimensional Scaling (MDS) == PCA
//!
//! Key steps:
//!
//! 1. Distance Matrix. PCoA starts with a distance or dissimilarity matrix that quantifies how
//!    similar or different the data points are from each other. Dissimilarity can take various forms
//!    as long as it meets the criteria of a metric space. `$D=[d_{ij}]$`
//!
//! 2. Double-Centering. `$B=-{\frac{1}{2}}CD^{(2)}C$` using the centering matrix
//!    `$C=I-{\frac{1}{n}}J_{n}$` where `$n$` is the number of objects, `$I$` is the `$n \times n$`
//!    identity matrix, and `$J_{n}$` is an `$n \times n$` matrix of all ones.
//!
//! 3. Eigenvalue Decomposition. The eigenvalues indicate the amount of variance captured by each
//!    dimension, while the eigenvectors provide the coordinates of the data points in the new space.
//!
//! ## Dimensionality Reduction
//!
//! PCoA transforms high-dimensional data into a lower-dimensional space (typically 2D or 3D)
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

// pub struct Mds {
//     dim: usize,
// }
//
// impl Mds {
//     pub fn new(dim: usize) -> Self {
//         Mds { dim }
//     }
//
//     pub fn double_centering(&self, matrix: &faer::Mat<f64>) -> faer::Mat<f64> {
//         let mut centered = matrix.clone();
//
//         let ncol = centered.ncols();
//         let mut col_mean = faer::Row::zeros(ncol);
//         faer::stats::row_mean(
//             col_mean.as_mut(),
//             centered.as_ref(),
//             faer::stats::NanHandling::Ignore,
//         );
//
//         let nrow = centered.nrows();
//         let mut row_mean = faer::Row::zeros(nrow);
//         faer::stats::row_mean(
//             row_mean.as_mut(),
//             centered.as_ref(),
//             faer::stats::NanHandling::Ignore,
//         );
//
//         let grand_mean = centered.mean();
//
//         for j in 0..ncol {
//             for i in 0..nrow {
//                 centered[(i, j)] -= row_mean[i] + col_mean[j] - grand_mean;
//             }
//         }
//
//         centered
//     }
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn test_centered() {
//         let mut matrix = faer::mat![
//             [0., 7., 5., 5.],
//             [7., 0., 4., 9.],
//             [5., 4., 0., 3.],
//             [5., 9., 3., 0.],
//         ];
//         let exp: faer::mat![
//             [11.9375, -6.6875, -6.6875, 1.4375],
//             [-6.6875, 23.6875, 3.6875, -20.6875],
//             [-6.6875, 3.6875, -0.3125, 3.3125],
//             [1.4375, -20.6875, 3.3125, 15.9375],
//         ];
//
//         eprintln!("matrix = {:#?}", matrix);
//         let mut mds = Mds::new(2);
//         let res = mds.double_centering(&mut matrix);
//         eprintln!("matrix = {:#?}", matrix);
//         assert_eq!(res, exp);
//     }
// }
