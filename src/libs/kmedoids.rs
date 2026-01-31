//! Implementation of the K-Medoids clustering algorithm.
//!
//! K-Medoids is a clustering algorithm related to K-Means and the Medoidshift algorithm.
//! Both the K-Means and K-Medoids algorithms are partitional (breaking the dataset up into groups)
//! and attempt to minimize the distance between points labeled to be in a cluster and a point
//! designated as the center of that cluster. In contrast to the K-Means algorithm, K-Medoids
//! chooses actual data points as centers (medoids), and is thereby more robust to noise and outliers.
//!
//! # Example
//!
//! ```
//! use hnsm::KMedoids;
//! use intspan::ScoringMatrix;
//!
//! // Create a distance matrix for 5 points
//! let mut sm = ScoringMatrix::<f32>::with_size_and_defaults(5, 0.0, 100.0);
//!
//! // Cluster 1: {0, 1}
//! sm.set(0, 1, 1.0);
//!
//! // Cluster 2: {2, 3, 4}
//! sm.set(2, 3, 1.0);
//! sm.set(2, 4, 1.0);
//! sm.set(3, 4, 1.0);
//!
//! // Distances between clusters are large (default 100.0)
//!
//! let kmedoids = KMedoids::new(2, 100, 5);
//! let clusters = kmedoids.perform_clustering(&sm);
//!
//! assert_eq!(clusters.len(), 2);
//! ```

use rand::prelude::*;

/// K-Medoids Clustering (Lloyd-like algorithm)
pub struct KMedoids {
    k: usize,
    max_iter: usize,
    runs: usize,
}

impl KMedoids {
    /// Create a new KMedoids instance
    ///
    /// # Arguments
    ///
    /// * `k` - Number of clusters
    /// * `max_iter` - Maximum number of iterations per run
    /// * `runs` - Number of random initializations to perform
    pub fn new(k: usize, max_iter: usize, runs: usize) -> Self {
        Self { k, max_iter, runs }
    }

    /// Perform clustering on the given distance matrix
    pub fn perform_clustering(&self, matrix: &intspan::ScoringMatrix<f32>) -> Vec<Vec<usize>> {
        let n = matrix.size();
        if n == 0 {
            return vec![];
        }
        if self.k >= n {
            return (0..n).map(|i| vec![i]).collect();
        }

        let mut best_cost = f32::MAX;
        let mut best_assignment = vec![0; n];

        let mut rng = rand::rng();
        let indices: Vec<usize> = (0..n).collect();

        for _ in 0..self.runs {
            // 1. Initialize medoids
            let mut medoids: Vec<usize> =
                indices.choose_multiple(&mut rng, self.k).cloned().collect();

            let mut assignment = vec![0; n];
            let mut iter = 0;

            // Loop until convergence or max_iter
            loop {
                let mut changed = false;

                // 2. Assignment step
                for i in 0..n {
                    let mut min_dist = f32::MAX;
                    let mut closest_c_idx = 0;

                    for (c_idx, &medoid) in medoids.iter().enumerate() {
                        let d = matrix.get(i, medoid);
                        if d < min_dist {
                            min_dist = d;
                            closest_c_idx = c_idx;
                        }
                    }
                    if assignment[i] != closest_c_idx {
                        assignment[i] = closest_c_idx;
                        changed = true;
                    }
                }

                if !changed || iter >= self.max_iter {
                    break;
                }

                // 3. Update step
                let mut clusters = vec![Vec::new(); self.k];
                for (i, &c_idx) in assignment.iter().enumerate() {
                    clusters[c_idx].push(i);
                }

                for (c_idx, points) in clusters.iter().enumerate() {
                    if points.is_empty() {
                        continue;
                    }

                    // Find new medoid (min sum of distances)
                    let mut min_sum_dist = f32::MAX;
                    let mut new_medoid = medoids[c_idx];

                    for &candidate in points {
                        let mut sum_dist = 0.0;
                        for &peer in points {
                            sum_dist += matrix.get(candidate, peer);
                        }
                        if sum_dist < min_sum_dist {
                            min_sum_dist = sum_dist;
                            new_medoid = candidate;
                        }
                    }
                    medoids[c_idx] = new_medoid;
                }

                iter += 1;
            }

            // Calculate total cost
            let mut total_cost = 0.0;
            for i in 0..n {
                let medoid = medoids[assignment[i]];
                total_cost += matrix.get(i, medoid);
            }

            if total_cost < best_cost {
                best_cost = total_cost;
                best_assignment = assignment;
            }
        }

        // Convert to result format
        let mut res_clusters = vec![Vec::new(); self.k];
        for (i, &c_idx) in best_assignment.iter().enumerate() {
            res_clusters[c_idx].push(i);
        }

        res_clusters.into_iter().filter(|c| !c.is_empty()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use intspan::ScoringMatrix;

    #[test]
    fn test_kmedoids_simple() {
        // Create a distance matrix for 4 points
        // 0-1 are close, 2-3 are close
        let mut sm = ScoringMatrix::<f32>::with_size_and_defaults(4, 0.0, 10.0);
        sm.set(0, 1, 1.0);
        sm.set(2, 3, 1.0);

        let kmedoids = KMedoids::new(2, 100, 10);
        let clusters = kmedoids.perform_clustering(&sm);

        assert_eq!(clusters.len(), 2);

        // Check content of clusters
        let c1 = &clusters[0];
        let c2 = &clusters[1];

        // One cluster should contain 0,1 and other 2,3
        let has_0 = c1.contains(&0) || c2.contains(&0);
        let has_2 = c1.contains(&2) || c2.contains(&2);
        assert!(has_0 && has_2);
    }

    #[test]
    fn test_kmedoids_trivial() {
        // k=1
        let sm = ScoringMatrix::<f32>::with_size_and_defaults(3, 0.0, 1.0);
        let kmedoids = KMedoids::new(1, 10, 1);
        let clusters = kmedoids.perform_clustering(&sm);

        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 3);
    }

    #[test]
    fn test_kmedoids_k_equals_n() {
        // k=n
        let sm = ScoringMatrix::<f32>::with_size_and_defaults(3, 0.0, 1.0);
        let kmedoids = KMedoids::new(3, 10, 1);
        let clusters = kmedoids.perform_clustering(&sm);

        assert_eq!(clusters.len(), 3);
    }
}
