//! Implementation of the [DBSCAN](https://en.wikipedia.org/wiki/DBSCAN) clustering algorithm.
//!
//! Key features:
//! * Density-based clustering
//! * Automatic noise detection
//! * No predefined cluster count
//! * Handles arbitrary cluster shapes
//!
//! Parameters:
//! * eps: Neighborhood radius
//! * min_points: Core point threshold
//!
//! Output formats:
//! * Cluster labels: Some(id) or None (noise)
//! * Cluster groups: Vec<Vec<point_indices>>
//! * Representative pairs: Vec<(center, member)>
// Adopt from https://blog.petrzemek.net/2017/01/01/implementing-dbscan-from-distance-matrix-in-rust/
use crate::ScoringMatrix;
use std::collections::{HashMap, VecDeque};

#[derive(Debug)]
pub struct Dbscan<T> {
    eps: T,
    min_points: usize,
    clusters: Vec<Option<usize>>,
    visited: Vec<bool>,
    current_cluster: usize,
}

impl<T> Dbscan<T>
where
    T: Default + Copy + PartialOrd + std::ops::AddAssign + num_traits::ToPrimitive,
{
    /// Creates a new DBSCAN instance.
    ///
    /// # Parameters
    ///
    /// * `eps` - The maximum distance between two points for them to be in the
    ///   same neighborhood.
    /// * `min_points` - The minimal number of points in a neighborhood for a
    ///   point to be considered as a core point.
    pub fn new(eps: T, min_points: usize) -> Self {
        Dbscan {
            eps,
            min_points,
            clusters: Vec::new(),
            visited: Vec::new(),
            current_cluster: 0,
        }
    }

    /// Performs DBSCAN clustering from the given distance matrix.
    ///
    /// # Returns
    ///
    /// Returns cluster labels for each point in the dataset. Noisy samples are
    /// set to `None`.
    ///
    /// ```
    /// # use hnsm::Dbscan;
    /// # use hnsm::ScoringMatrix;
    ///
    /// let mut dbscan = Dbscan::new(1, 2);
    /// let mut m = ScoringMatrix::<i8>::new(5, 0, 100);
    /// m.set(0, 1, 1);
    /// m.set(0, 2, 9);
    /// m.set(0, 3, 9);
    /// m.set(0, 4, 9);
    /// m.set(1, 2, 9);
    /// m.set(1, 3, 9);
    /// m.set(1, 4, 9);
    /// m.set(2, 3, 1);
    /// m.set(2, 4, 9);
    /// m.set(3, 4, 9);
    ///
    /// let clustering = dbscan.perform_clustering(&m);
    ///
    /// assert_eq!(clustering[0], Some(0));
    /// assert_eq!(clustering[1], Some(0));
    /// assert_eq!(clustering[2], Some(1));
    /// assert_eq!(clustering[3], Some(1));
    /// assert_eq!(clustering[4], None);
    /// ```
    ///
    /// In the above example, points `0` and `1` form a single cluster, points
    /// `2` and `3` form a different cluster, and point `4` does not belong any
    /// cluster (it is a noise point).
    pub fn perform_clustering(&mut self, matrix: &ScoringMatrix<T>) -> &Vec<Option<usize>> {
        self.clusters = vec![None; matrix.size()];
        self.visited = vec![false; matrix.size()];
        self.current_cluster = 0;

        for point in 0..matrix.size() {
            if self.visited[point] {
                continue;
            }

            self.visited[point] = true;
            let neighbors = self.region_query(matrix, point);
            if neighbors.len() >= self.min_points {
                self.expand_cluster(matrix, point, neighbors);
                self.current_cluster += 1;
            }
        }

        self.clusters.as_ref()
    }

    fn expand_cluster(
        &mut self,
        matrix: &ScoringMatrix<T>,
        point: usize,
        mut neighbors: VecDeque<usize>,
    ) {
        self.clusters[point] = Some(self.current_cluster);

        while let Some(other_point) = neighbors.pop_front() {
            if !self.visited[other_point] {
                self.visited[other_point] = true;
                let mut other_neighbors = self.region_query(matrix, other_point);
                if other_neighbors.len() >= self.min_points {
                    neighbors.append(&mut other_neighbors);
                }
            }
            if self.clusters[other_point].is_none() {
                self.clusters[other_point] = Some(self.current_cluster);
            }
        }
    }

    fn region_query(&self, matrix: &ScoringMatrix<T>, point: usize) -> VecDeque<usize> {
        let mut neighbors = VecDeque::new();
        for other_point in 0..matrix.size() {
            let dist = matrix.get(point, other_point);
            if dist <= self.eps {
                neighbors.push_back(other_point);
            }
        }
        neighbors
    }

    fn all_clusters(&self) -> (HashMap<usize, Vec<usize>>, Vec<usize>) {
        let mut cluster_map: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut noise_points: Vec<usize> = Vec::new();

        for (point, cluster) in self.clusters.iter().enumerate() {
            match cluster {
                Some(cluster_id) => {
                    cluster_map.entry(*cluster_id).or_default().push(point);
                }
                None => {
                    noise_points.push(point);
                }
            }
        }
        (cluster_map, noise_points)
    }

    pub fn results_cluster(&self) -> Vec<Vec<usize>> {
        let (cluster_map, noise_points) = self.all_clusters();
        let mut res: Vec<Vec<usize>> = vec![];

        for (_, points) in cluster_map.iter() {
            res.push(points.clone());
        }
        for p in noise_points {
            res.push(vec![p]);
        }

        res
    }

    /// Finds and prints the representative point of each cluster.
    pub fn results_pair(&self, matrix: &ScoringMatrix<T>) -> Vec<(usize, usize)> {
        let (cluster_map, noise_points) = self.all_clusters();

        // representative point, point
        let mut res: Vec<(usize, usize)> = vec![];

        for (_, points) in cluster_map.iter() {
            let mut sum_distance_of: HashMap<usize, f64> = HashMap::new();
            for &point in points {
                // Calculate the sum of distances from this point to all others in the cluster
                let mut sum_distance = T::default();
                for &other_point in points {
                    sum_distance += matrix.get(point, other_point);
                }
                sum_distance_of.insert(point, sum_distance.to_f64().unwrap());
            }
            let representative = sum_distance_of
                .iter()
                .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(&key, _)| key)
                .unwrap();

            for &point in points {
                res.push((representative, point));
            }
        }

        for p in noise_points {
            res.push((p, p));
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_points_are_in_single_cluster_when_their_distance_is_zero() {
        let mut dbscan = Dbscan::new(1, 2);
        let m = ScoringMatrix::<i8>::new(2, 0, 1);

        let clustering = dbscan.perform_clustering(&m);

        assert_eq!(clustering[0], Some(0));
        assert_eq!(clustering[1], Some(0));
    }

    #[test]
    fn test_points_are_correctly_clustered_based_on_their_distance() {
        let mut dbscan = Dbscan::new(1, 2);
        let mut m = ScoringMatrix::<i8>::new(5, 0, 100);
        m.set(0, 1, 1);
        m.set(0, 2, 9);
        m.set(0, 3, 9);
        m.set(0, 4, 9);
        m.set(1, 2, 9);
        m.set(1, 3, 9);
        m.set(1, 4, 9);
        m.set(2, 3, 1);
        m.set(2, 4, 9);
        m.set(3, 4, 9);

        let clustering = dbscan.perform_clustering(&m);

        assert_eq!(clustering[0], Some(0));
        assert_eq!(clustering[1], Some(0));
        assert_eq!(clustering[2], Some(1));
        assert_eq!(clustering[3], Some(1));
        assert_eq!(clustering[4], None);
    }

    #[test]
    fn test_neighboring_points_are_put_into_cluster_even_if_they_have_been_visited() {
        // In 2D, the points in this test are placed as follows:
        //
        //    0
        //      1
        //        2
        //
        // Epsilon is set to 1 and min_points to 3. When the first point is
        // checked (0), it is marked as visited. Since it has only a single
        // neighbor, the two points (0 and 1) cannot form a cluster because
        // min_points is 3. Then, the algorithm continues to point 1. It has
        // two neighbors (0 and
        // 2), so the three points (0, 1, 2) can form a cluster. In this test,
        // we ensure that even when the first point (0) has already been
        // marked as visited, it is put into the cluster because it is not
        // yet a member of any other cluster.
        let mut dbscan = Dbscan::new(1, 3);
        let mut m = ScoringMatrix::<i8>::new(3, 0, 100);
        m.set(0, 1, 1);
        m.set(0, 2, 2);
        m.set(1, 2, 1);

        let clustering = dbscan.perform_clustering(&m);

        assert_eq!(clustering[0], Some(0));
        assert_eq!(clustering[1], Some(0));
        assert_eq!(clustering[2], Some(0));
    }

    #[test]
    fn test_points_that_do_not_belong_to_any_cluster_are_none() {
        let mut dbscan = Dbscan::new(1, 2);
        let m = ScoringMatrix::<i8>::new(1, 0, 100);

        let clustering = dbscan.perform_clustering(&m);

        assert_eq!(clustering[0], None);
    }
}
