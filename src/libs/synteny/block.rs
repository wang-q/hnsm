use crate::libs::synteny::graph::SyntenyGraph;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct BlockRange {
    pub seq_id: u32,
    pub start: u32,
    pub end: u32,
    pub strand: bool,
    pub count: usize,
}

#[derive(Debug, Clone)]
pub struct SyntenyBlock {
    pub ranges: HashMap<u32, BlockRange>,
}

impl SyntenyBlock {
    pub fn new() -> Self {
        Self {
            ranges: HashMap::new(),
        }
    }

    /// Reconstruct synteny blocks from a linear path of minimizer hashes.
    pub fn from_path(graph: &SyntenyGraph, path: &[u64]) -> Self {
        let mut block = SyntenyBlock::new();
        if path.is_empty() {
            return block;
        }

        let node_indices: Vec<NodeIndex> = path
            .iter()
            .filter_map(|h| graph.node_map.get(h).cloned())
            .collect();

        if node_indices.len() != path.len() {
            return block;
        }

        // Special case for single node
        if node_indices.len() == 1 {
            let u = node_indices[0];
            let node = &graph.graph[u];
            for occ in &node.occurrences {
                block.update_range(occ.seq_id, occ.pos, occ.pos, occ.strand);
            }
            return block;
        }

        for i in 0..node_indices.len() - 1 {
            let u = node_indices[i];
            let v = node_indices[i + 1];
            let u_node = &graph.graph[u];
            let v_node = &graph.graph[v];

            for edge in graph.graph.edges_directed(u, petgraph::Direction::Outgoing) {
                if edge.target() == v {
                    let weight = edge.weight();
                    let seq_id = weight.seq_id;
                    let dist = weight.distance;

                    // Find occurrence pair (occ_u, occ_v) such that occ_v.pos - occ_u.pos == dist
                    // Optimization: Use binary search since occurrences are sorted by seq_id then pos
                    
                    let find_range = |occs: &[crate::libs::synteny::graph::Occurrence], sid: u32| -> std::ops::Range<usize> {
                        let start = occs.partition_point(|occ| occ.seq_id < sid);
                        let len = occs[start..].partition_point(|occ| occ.seq_id <= sid);
                        start..start + len
                    };

                    let u_range = find_range(&u_node.occurrences, seq_id);
                    let v_range = find_range(&v_node.occurrences, seq_id);

                    if u_range.is_empty() || v_range.is_empty() {
                        continue;
                    }

                    let u_slice = &u_node.occurrences[u_range];
                    let v_slice = &v_node.occurrences[v_range];

                    let mut start_pos = None;
                    let mut end_pos = None;
                    let mut strand = true;

                    for occ_u in u_slice {
                        // We need occ_v.pos == occ_u.pos + dist
                        let target_pos = occ_u.pos + dist;
                        
                        if let Ok(idx) = v_slice.binary_search_by_key(&target_pos, |occ| occ.pos) {
                            start_pos = Some(occ_u.pos);
                            end_pos = Some(v_slice[idx].pos);
                            strand = occ_u.strand;
                            break;
                        }
                    }

                    if let (Some(s), Some(e)) = (start_pos, end_pos) {
                        block.update_range(seq_id, s, e, strand);
                    }
                }
            }
        }

        block
    }

    fn update_range(&mut self, seq_id: u32, pos1: u32, pos2: u32, strand: bool) {
        let start = std::cmp::min(pos1, pos2);
        let end = std::cmp::max(pos1, pos2);

        self.ranges
            .entry(seq_id)
            .and_modify(|r| {
                if start < r.start {
                    r.start = start;
                }
                if end > r.end {
                    r.end = end;
                }
                r.count += 1;
            })
            .or_insert(BlockRange {
                seq_id,
                start,
                end,
                strand,
                count: 1,
            });
    }
}
