use crate::libs::hash::MinimizerInfo;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::{EdgeRef, NodeIndexable};
use petgraph::Direction;
use std::collections::{HashMap, HashSet, VecDeque};

#[derive(Debug, Clone)]
pub struct Node {
    pub hash: u64,
    pub occurrences: Vec<Occurrence>,
}

#[derive(Debug, Clone, Copy)]
pub struct Occurrence {
    pub seq_id: u32,
    pub pos: u32,
    pub strand: bool,
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub seq_id: u32,
    pub distance: u32,
}

pub struct SyntenyGraph {
    pub graph: DiGraph<Node, Edge>,
    pub node_map: HashMap<u64, NodeIndex>,
}

impl SyntenyGraph {
    pub fn new() -> Self {
        Self {
            graph: DiGraph::new(),
            node_map: HashMap::new(),
        }
    }

    /// Add minimizers from a single sequence to the graph.
    /// Assumes minimizers are sorted by position.
    pub fn add_minimizers(&mut self, minimizers: &[MinimizerInfo], max_indel: u32) {
        if minimizers.is_empty() {
            return;
        }

        let mut prev_idx: Option<NodeIndex> = None;
        let mut prev_info: Option<&MinimizerInfo> = None;

        for info in minimizers {
            let idx = *self.node_map.entry(info.hash).or_insert_with(|| {
                self.graph.add_node(Node {
                    hash: info.hash,
                    occurrences: Vec::new(),
                })
            });

            // Add occurrence
            self.graph[idx].occurrences.push(Occurrence {
                seq_id: info.seq_id,
                pos: info.pos,
                strand: info.strand,
            });

            // Add edge from previous minimizer
            if let (Some(u), Some(prev)) = (prev_idx, prev_info) {
                if prev.seq_id == info.seq_id {
                    let dist = info.pos.saturating_sub(prev.pos);
                    if dist <= max_indel {
                        self.graph.add_edge(
                            u,
                            idx,
                            Edge {
                                seq_id: info.seq_id,
                                distance: dist,
                            },
                        );
                    }
                }
            }

            prev_idx = Some(idx);
            prev_info = Some(info);
        }
    }

    /// Prune edges that have low support (fewer than `min_weight` occurrences).
    /// Since we store one edge per occurrence, we count parallel edges.
    pub fn prune_low_weight_edges(&mut self, min_weight: usize) {
        let mut good_connections = HashSet::new();

        // Iterate over all nodes to find edges
        for node_idx in self.graph.node_indices() {
            // Group edges by target
            // Use a map to count: Target -> Count
            let mut target_counts: HashMap<NodeIndex, usize> = HashMap::new();

            for edge in self.graph.edges_directed(node_idx, Direction::Outgoing) {
                *target_counts.entry(edge.target()).or_insert(0) += 1;
            }

            for (target, count) in target_counts {
                if count >= min_weight {
                    good_connections.insert((node_idx, target));
                }
            }
        }

        // Use retain_edges to safely remove edges
        // We keep an edge if (source, target) is in good_connections
        self.graph.retain_edges(|g, e| {
            let (source, target) = g.edge_endpoints(e).unwrap();
            good_connections.contains(&(source, target))
        });
    }

    /// Perform transitive reduction on the graph.
    /// Removes edge u -> v if there is another path from u to v.
    pub fn transitive_reduction(&mut self) {
        // Iterate over all nodes
        for node in self.graph.node_indices() {
            let children_set: HashSet<NodeIndex> = self.graph.neighbors(node).collect();
            if children_set.len() < 2 {
                continue;
            }
            let mut children: Vec<NodeIndex> = children_set.into_iter().collect();
            children.sort(); // Deterministic order

            // Track removed children to handle mutual dependencies (cycles) correctly
            let mut removed_children = HashSet::new();

            // Check each child
            for (i, &child) in children.iter().enumerate() {
                // Try to find a path from any OTHER child to this child
                // If u -> other -> ... -> child exists, then u -> child is redundant
                for (j, &other) in children.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    // If we already decided to remove u->other, we can't use it as a bridge
                    if removed_children.contains(&other) {
                        continue;
                    }

                    // Limit BFS depth to avoid performance hit on large graphs
                    // Synteny skips shouldn't be excessively long
                    if self.has_path(other, child, 50) {
                        // Found path u -> other -> ... -> child
                        // And u -> other exists (not removed)
                        // So u -> child is redundant

                        // Remove all edges between node and child
                        while let Some(e) = self.graph.find_edge(node, child) {
                            self.graph.remove_edge(e);
                        }
                        removed_children.insert(child);
                        break; // Edge removed, move to next child
                    }
                }
            }
        }
    }

    // Simple BFS with depth limit
    fn has_path(&self, start: NodeIndex, end: NodeIndex, max_depth: usize) -> bool {
        if start == end {
            return true;
        }

        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        queue.push_back((start, 0));
        visited.insert(start);

        while let Some((curr, depth)) = queue.pop_front() {
            if curr == end {
                return true;
            }
            if depth >= max_depth {
                continue;
            }

            for neighbor in self.graph.neighbors(curr) {
                if visited.insert(neighbor) {
                    queue.push_back((neighbor, depth + 1));
                }
            }
        }
        false
    }

    /// Find linear paths (synteny blocks) in the graph.
    /// A linear path is a sequence of nodes v1 -> v2 -> ... -> vn where:
    /// - v_i has exactly one outgoing neighbor v_{i+1} (after pruning)
    /// - v_{i+1} has exactly one incoming neighbor v_i
    /// Returns a list of paths, where each path is a list of minimizer hashes.
    pub fn get_linear_paths(&self) -> Vec<Vec<u64>> {
        let node_bound = self.graph.node_bound();
        let invalid = NodeIndex::end(); 
        
        // Adjacency tables: None=0, Some(x)=1 unique, Some(invalid)=many/conflict
        let mut adj_out = vec![None; node_bound];
        let mut adj_in = vec![None; node_bound];
        
        // 1. Global Edge Scan (O(E))
        // Collapses parallel edges and detects branching
        for edge in self.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let u_idx = u.index();
            let v_idx = v.index();
            
            // Update Out u
            match adj_out[u_idx] {
                None => adj_out[u_idx] = Some(v),
                Some(curr) => {
                    if curr != v && curr != invalid {
                        adj_out[u_idx] = Some(invalid); // Branching out
                    }
                }
            }
            
            // Update In v
            match adj_in[v_idx] {
                None => adj_in[v_idx] = Some(u),
                Some(curr) => {
                    if curr != u && curr != invalid {
                        adj_in[v_idx] = Some(invalid); // Branching in
                    }
                }
            }
        }
        
        let mut paths = Vec::new();
        let mut visited = vec![false; node_bound];
        let all_nodes: Vec<NodeIndex> = self.graph.node_indices().collect();
        
        // 2. Find paths starting from valid heads
        for &node in &all_nodes {
            if visited[node.index()] { continue; }
            
            let u_in = adj_in[node.index()];
            // Skip isolated nodes (no in, no out)
            if u_in.is_none() && adj_out[node.index()].is_none() {
                continue;
            }
            
            // Is Start Node?
            let is_start = if let Some(parent) = u_in {
                if parent == invalid {
                    true // Multiple parents -> Start (Merge point)
                } else {
                    // Single parent. Check if parent branches.
                    let p_out = adj_out[parent.index()];
                    if p_out == Some(invalid) {
                        true // Parent branches -> Start
                    } else if p_out != Some(node) {
                        true // Parent points elsewhere
                    } else {
                        false // Internal node (1-to-1)
                    }
                }
            } else {
                true // 0 parents -> Start
            };
            
            if is_start {
                let mut path = Vec::new();
                let mut curr = node;
                
                loop {
                    if visited[curr.index()] { break; }
                    visited[curr.index()] = true;
                    path.push(self.graph[curr].hash);
                    
                    // Move next
                    if let Some(next) = adj_out[curr.index()] {
                        if next == invalid { break; } // Branching out
                        
                        let next_in = adj_in[next.index()];
                        if next_in == Some(invalid) {
                            break; // Next has multiple parents
                        }
                        
                        // Check if next's unique parent is us
                        if next_in != Some(curr) {
                             break;
                        }
                        
                        curr = next;
                    } else {
                        break; // End of path
                    }
                }
                
                if !path.is_empty() {
                    paths.push(path);
                }
            }
        }
        
        // 3. Handle Pure Cycles (Rings) or Remnants
        for &node in &all_nodes {
            if visited[node.index()] { continue; }
            
            // Skip isolated nodes here too!
            if adj_in[node.index()].is_none() && adj_out[node.index()].is_none() {
                continue;
            }
            
            let mut path = Vec::new();
            let mut curr = node;
            
            loop {
                if visited[curr.index()] { break; }
                visited[curr.index()] = true;
                path.push(self.graph[curr].hash);
                
                if let Some(next) = adj_out[curr.index()] {
                     if next == invalid { break; } 
                     curr = next;
                } else {
                    break;
                }
            }
            if !path.is_empty() {
                paths.push(path);
            }
        }
        
        paths
    }
}
