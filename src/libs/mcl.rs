//! Implementation of the [MCL](https://micans.org/mcl/) (Markov Clustering) algorithm.
//!
//! MCL is a fast and scalable unsupervised cluster algorithm for graphs based on
//! simulation of (stochastic) flow in graphs.
//!
//! Key concepts:
//! * **Expansion**: Squaring the matrix to simulate flow spreading (random walks).
//! * **Inflation**: Raising elements to a power to boost strong connections and weaken weak ones.
//!
//! # Example
//!
//! ```
//! use hnsm::Mcl;
//! use intspan::ScoringMatrix;
//!
//! // Create a matrix for 4 nodes (A, B, C, D)
//! // Defaults: self-loop=1.0, missing=0.0
//! let mut sm = ScoringMatrix::<f32>::with_size_and_defaults(4, 1.0, 0.0);
//!
//! // Connect A-B
//! sm.set(0, 1, 1.0);
//! sm.set(1, 0, 1.0);
//!
//! // Connect C-D
//! sm.set(2, 3, 1.0);
//! sm.set(3, 2, 1.0);
//!
//! let mcl = Mcl::new(2.0);
//! let clusters = mcl.perform_clustering(&sm);
//!
//! // Should find 2 clusters: {0,1} and {2,3}
//! assert_eq!(clusters.len(), 2);
//! ```

use std::collections::HashMap;

/// Markov Clustering Algorithm configuration and execution
pub struct Mcl {
    inflation: f64,
    prune_limit: f64,
    max_iter: usize,
}

impl Mcl {
    /// Create a new MCL instance with specified inflation parameter.
    ///
    /// # Arguments
    ///
    /// * `inflation` - Controls the granularity of clusters.
    ///   * Low values (e.g. 1.4) lead to coarser clusters.
    ///   * High values (e.g. 4.0) lead to tighter/more clusters.
    ///   * A common default is 2.0.
    pub fn new(inflation: f64) -> Self {
        Self {
            inflation,
            prune_limit: 1e-5,
            max_iter: 100,
        }
    }

    /// Set the threshold for pruning small values in the matrix.
    ///
    /// Matrix entries smaller than this value will be set to zero to maintain sparsity
    /// and improve performance. Default is 1e-5.
    pub fn set_prune_limit(&mut self, limit: f64) {
        self.prune_limit = limit;
    }

    /// Set the maximum number of iterations for the algorithm.
    ///
    /// If convergence is not reached within this limit, the algorithm stops.
    /// Default is 100.
    pub fn set_max_iter(&mut self, max_iter: usize) {
        self.max_iter = max_iter;
    }

    /// Perform MCL clustering on the given ScoringMatrix.
    ///
    /// # Returns
    ///
    /// A vector of clusters, where each cluster is a vector of node indices
    /// corresponding to the input `ScoringMatrix`.
    pub fn perform_clustering(&self, sm: &intspan::ScoringMatrix<f32>) -> Vec<Vec<usize>> {
        let mut matrix = SparseMat::from_scoring_matrix(sm);
        matrix.normalize();

        let mut changed = true;
        let mut iter = 0;

        while changed && iter < self.max_iter {
            let prev_matrix = matrix.clone();

            // Expansion (Power 2)
            matrix = matrix.expand();

            // Inflation (Element-wise power + Normalize)
            matrix.inflate(self.inflation);

            // Pruning
            matrix.prune(self.prune_limit);

            if matrix.is_converged(&prev_matrix) {
                changed = false;
            }
            iter += 1;
        }

        // Interpret Clusters
        let mut graph = petgraph::graphmap::UnGraphMap::<usize, ()>::new();
        // Add edges for all non-zero entries in the result matrix
        // The attractors (nodes with self-loops) will gather their attracted nodes
        for j in 0..matrix.size {
            for &(i, _) in &matrix.cols[j] {
                graph.add_edge(i, j, ());
            }
        }

        petgraph::algo::tarjan_scc(&graph)
    }
}

// Simple sparse matrix (Column-Major)
#[derive(Clone)]
struct SparseMat {
    size: usize,
    cols: Vec<Vec<(usize, f64)>>,
}

impl SparseMat {
    fn from_scoring_matrix(sm: &intspan::ScoringMatrix<f32>) -> Self {
        let size = sm.size();
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); size];

        // Iterating N^2 is required as ScoringMatrix is generic
        // Assuming reasonably sparse or small N
        for i in 0..size {
            for j in 0..size {
                let val = sm.get(i, j) as f64;
                if val.abs() > 1e-5 {
                    cols[j].push((i, val));
                }
            }
        }
        // Ensure sorted rows
        for col in &mut cols {
            col.sort_by_key(|(r, _)| *r);
        }
        Self { size, cols }
    }

    fn normalize(&mut self) {
        for col in &mut self.cols {
            let sum: f64 = col.iter().map(|(_, v)| *v).sum();
            if sum > 1e-9 {
                for (_, v) in col.iter_mut() {
                    *v /= sum;
                }
            }
        }
    }

    fn expand(&self) -> Self {
        let mut new_cols = vec![Vec::new(); self.size];

        // M_new = M * M
        // Col j of M_new = M * col_j(M)
        for j in 0..self.size {
            let mut accumulator: HashMap<usize, f64> = HashMap::new();

            // For each non-zero entry (k, val_k) in col j of M
            for &(k, val_k) in &self.cols[j] {
                // Add val_k * col_k(M) to accumulator
                if let Some(col_k) = self.cols.get(k) {
                    for &(row_i, val_i) in col_k {
                        *accumulator.entry(row_i).or_insert(0.0) += val_i * val_k;
                    }
                }
            }

            // Convert accumulator to vec
            let mut col: Vec<(usize, f64)> = accumulator.into_iter().collect();
            // Sort by row index
            col.sort_by_key(|(r, _)| *r);
            new_cols[j] = col;
        }
        Self {
            size: self.size,
            cols: new_cols,
        }
    }

    fn inflate(&mut self, power: f64) {
        for col in &mut self.cols {
            for (_, v) in col.iter_mut() {
                *v = v.powf(power);
            }
        }
        self.normalize();
    }

    fn prune(&mut self, threshold: f64) {
        for col in &mut self.cols {
            col.retain(|(_, v)| *v > threshold);
        }
    }

    fn is_converged(&self, other: &Self) -> bool {
        if self.size != other.size {
            return false;
        }

        // Compare structure and values
        for j in 0..self.size {
            let col1 = &self.cols[j];
            let col2 = &other.cols[j];

            if col1.len() != col2.len() {
                return false;
            }

            for ((r1, v1), (r2, v2)) in col1.iter().zip(col2.iter()) {
                if r1 != r2 || (*v1 - *v2).abs() > 1e-5 {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use intspan::ScoringMatrix;

    #[test]
    fn test_mcl_simple_clusters() {
        // Create 2 disconnected cliques: A-B-C and D-E
        let mut sm = ScoringMatrix::<f32>::with_size_and_defaults(5, 1.0, 0.0);

        // Clique 1: 0, 1, 2
        sm.set(0, 1, 1.0);
        sm.set(1, 0, 1.0);
        sm.set(0, 2, 1.0);
        sm.set(2, 0, 1.0);
        sm.set(1, 2, 1.0);
        sm.set(2, 1, 1.0);

        // Clique 2: 3, 4
        sm.set(3, 4, 1.0);
        sm.set(4, 3, 1.0);

        let mcl = Mcl::new(2.0);
        let clusters = mcl.perform_clustering(&sm);

        assert_eq!(clusters.len(), 2);

        // Verify cluster contents
        let mut c1 = clusters.iter().find(|c| c.contains(&0)).unwrap().clone();
        c1.sort();
        assert_eq!(c1, vec![0, 1, 2]);

        let mut c2 = clusters.iter().find(|c| c.contains(&3)).unwrap().clone();
        c2.sort();
        assert_eq!(c2, vec![3, 4]);
    }

    #[test]
    fn test_mcl_parameter_config() {
        let mut mcl = Mcl::new(2.0);
        mcl.set_prune_limit(1e-4);
        mcl.set_max_iter(50);

        assert_eq!(mcl.prune_limit, 1e-4);
        assert_eq!(mcl.max_iter, 50);
    }
}
