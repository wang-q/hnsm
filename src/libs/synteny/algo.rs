use crate::libs::bloom::BloomFilter;
use crate::libs::hash::{seq_sketch, MinimizerInfo};
use crate::libs::synteny::block::SyntenyBlock;
use crate::libs::synteny::graph::SyntenyGraph;
use intspan::IntSpan;
use log::info;
use std::collections::HashMap;

pub struct SyntenyFinder {
    pub k: usize,
    pub rounds: Vec<usize>,
    pub min_weight: usize,
    pub max_freq: u32,
    pub block_size: usize,
    pub chain_gap: u32,
}

impl SyntenyFinder {
    pub fn new(
        k: usize,
        rounds: Vec<usize>,
        min_weight: usize,
        max_freq: u32,
        block_size: usize,
        chain_gap: u32,
    ) -> Self {
        Self {
            k,
            rounds,
            min_weight,
            max_freq,
            block_size,
            chain_gap,
        }
    }

    /// Run synteny detection with iterative refinement.
    ///
    /// # Arguments
    /// * `provider` - A closure that iterates over sequences.
    ///                It takes a callback `emit(name, seq)` and calls it for each sequence.
    /// * `callback` - A closure called for each identified synteny block.
    pub fn run<P, F>(&self, mut provider: P, mut callback: F) -> anyhow::Result<()>
    where
        P: FnMut(&mut dyn FnMut(&str, &[u8])) -> anyhow::Result<()>,
        F: FnMut(usize, &SyntenyBlock),
    {
        let mut coverage: HashMap<u32, IntSpan> = HashMap::new();

        for &raw_w in &self.rounds {
            let mut w = raw_w;
            // Ensure window size is odd for symmetric minimizers if required by implementation
            if w % 2 == 0 {
                w += 1;
            }

            info!("Starting round with window size w={}", w);

            // 1. Count minimizer frequencies (First Pass)
            // Use Bloom Filter to filter out singletons
            let mut counts: HashMap<u64, u32> = HashMap::new();
            let mut bloom = BloomFilter::new(100_000_000, 0.01);
            let mut global_seq_id = 0;
            let mut total_minimizers = 0;

            info!("Pass 1: Counting minimizers...");

            provider(&mut |_, seq| {
                global_seq_id += 1;

                let is_covered = |pos: u32| {
                    coverage
                        .get(&global_seq_id)
                        .map(|s| s.contains(pos as i32))
                        .unwrap_or(false)
                };

                // Use a permissive filter for counting
                for m in seq_sketch(seq, global_seq_id, self.k, w, |_| true) {
                    if !is_covered(m.pos) {
                        total_minimizers += 1;
                        if bloom.contains(m.hash) {
                            counts.entry(m.hash).and_modify(|c| *c += 1).or_insert(2);
                        } else {
                            bloom.insert(m.hash);
                        }
                    }
                }
            })?;

            info!("Total minimizers: {}", total_minimizers);
            info!("Repetitive minimizers (frequency >= 2): {}", counts.len());
            
            info!("Pass 2: Building graph...");

            // 2. Build Graph (Second Pass)
            let mut graph = SyntenyGraph::new();
            global_seq_id = 0;

            provider(&mut |_, seq| {
                global_seq_id += 1;

                let is_covered = |pos: u32| {
                    coverage
                        .get(&global_seq_id)
                        .map(|s| s.contains(pos as i32))
                        .unwrap_or(false)
                };

                // Filter by frequency
                let mins = seq_sketch(seq, global_seq_id, self.k, w, |h| {
                    if let Some(&c) = counts.get(&h) {
                        c <= self.max_freq
                    } else {
                        false
                    }
                });

                // Filter by mask and add to graph
                let mins_filtered: Vec<MinimizerInfo> =
                    mins.into_iter().filter(|m| !is_covered(m.pos)).collect();

                graph.add_minimizers(&mins_filtered, self.chain_gap);
            })?;

            info!("Graph built. Edges: {}", graph.graph.edge_count());
            info!(
                "Pruning low weight edges (min_weight={})...",
                self.min_weight
            );

            // 3. Prune low weight edges
            graph.prune_low_weight_edges(self.min_weight);

            info!("Edges after pruning: {}", graph.graph.edge_count());
            
            // 3.5 Transitive Reduction
            info!("Performing transitive reduction...");
            graph.transitive_reduction();
            info!("Edges after transitive reduction: {}", graph.graph.edge_count());

            info!("Finding linear paths...");

            // 4. Find linear paths and convert to blocks
            let paths = graph.get_linear_paths();
            info!("Found {} linear paths. Converting to blocks...", paths.len());
            let mut blocks_found = 0;

            for path in paths {
                if path.is_empty() {
                    continue;
                }
                let block = SyntenyBlock::from_path(&graph, &path);

                // Filter by block size
                let max_len = block
                    .ranges
                    .values()
                    .map(|r| r.end.saturating_sub(r.start))
                    .max()
                    .unwrap_or(0);

                if (max_len as usize) < self.block_size {
                    continue;
                }

                // Update mask
                for (seq_id, range) in &block.ranges {
                    coverage
                        .entry(*seq_id)
                        .or_default()
                        .add_pair(range.start as i32, range.end as i32);
                }

                callback(w, &block);
                blocks_found += 1;
            }

            info!("Round complete. Found {} blocks.", blocks_found);
        }
        Ok(())
    }
}
