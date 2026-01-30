// use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub struct ChainOpt {
    pub gap_open_penalty: f32,
    pub gap_extension_penalty: f32,
    pub bp_gap_size: i32,
    pub max_match_score: f32,
    pub max_dist_between_matches: i32,
    pub min_alignment_score: f32,
}

impl Default for ChainOpt {
    fn default() -> Self {
        Self {
            gap_open_penalty: -1.0,
            gap_extension_penalty: -5.0,
            bp_gap_size: 10000,
            max_match_score: 50.0,
            max_dist_between_matches: 100000,
            min_alignment_score: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Anchor {
    pub id: usize, // identifier to map back to original data
    pub x: i32,
    pub y: i32,
    pub score: f32,
}

#[derive(Debug, Clone)]
pub struct Chain {
    pub indices: Vec<usize>,
    pub score: f32,
    pub path_scores: Vec<f32>, // Accumulated score at each step
}

pub struct DagChainer {
    pub options: ChainOpt,
}

impl DagChainer {
    pub fn new(options: ChainOpt) -> Self {
        Self { options }
    }

    /// Find chains in a set of anchors.
    /// Anchors must be sorted by X then Y before calling this function.
    pub fn find_chains(&self, anchors: &[Anchor]) -> Vec<Chain> {
        let n = anchors.len();
        if n == 0 {
            return Vec::new();
        }

        let mut is_used = vec![false; n];
        let mut chains = Vec::new();

        loop {
            let mut path_scores = vec![0.0; n];
            let mut from_indices = vec![-1; n];

            // DP: Compute scores for all unused nodes
            for j in 0..n {
                if is_used[j] {
                    continue;
                }

                path_scores[j] = anchors[j].score; // Base score

                // Look back at previous nodes
                for i in (0..j).rev() {
                    if is_used[i] {
                        continue;
                    }

                    let del_x = anchors[j].x - anchors[i].x - 1;
                    let del_y = anchors[j].y - anchors[i].y - 1;

                    // Must be increasing in both X and Y
                    if del_x < 0 || del_y < 0 {
                        continue;
                    }

                    // Check maximum distances
                    if del_x > self.options.max_dist_between_matches
                        || del_y > self.options.max_dist_between_matches
                    {
                        if del_x > self.options.max_dist_between_matches {
                            break;
                        }
                        continue;
                    }

                    let num_gaps = ((del_x + del_y + (del_x - del_y).abs()) as f32
                        / (2 * self.options.bp_gap_size) as f32
                        + 0.5) as i32;

                    let mut new_score = path_scores[i] + anchors[j].score;

                    if num_gaps > 0 {
                        new_score += self.options.gap_open_penalty
                            + (num_gaps as f32 * self.options.gap_extension_penalty);
                    }

                    if new_score > path_scores[j] {
                        path_scores[j] = new_score;
                        from_indices[j] = i as i32;
                    }
                }
            }

            // Find the best ending node among unused nodes
            let mut best_score = -1.0;
            let mut best_idx = usize::MAX;

            for (i, &score) in path_scores.iter().enumerate() {
                if !is_used[i] && score > best_score {
                    best_score = score;
                    best_idx = i;
                }
            }

            // Check if the best chain meets criteria
            if best_score < self.options.min_alignment_score {
                break;
            }

            // Reconstruct the best chain
            let mut indices = Vec::new();
            let mut current = best_idx;

            while from_indices[current] >= 0 {
                indices.push(current);
                current = from_indices[current] as usize;
            }
            indices.push(current);
            indices.reverse(); // Start -> End

            // Collect path scores
            let chain_path_scores: Vec<f32> = indices.iter().map(|&idx| path_scores[idx]).collect();

            // Mark nodes as used
            for &idx in &indices {
                is_used[idx] = true;
            }

            chains.push(Chain {
                indices,
                score: best_score,
                path_scores: chain_path_scores,
            });
        }

        chains
    }
}
