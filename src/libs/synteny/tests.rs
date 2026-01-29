use crate::libs::hash::MinimizerInfo;
use crate::libs::synteny::graph::SyntenyGraph;

#[test]
fn test_synteny_graph_linear_path() {
    let mut graph = SyntenyGraph::new();

    // Seq 1: 10 -> 20 -> 30 -> 40
    let seq1 = vec![
        MinimizerInfo {
            hash: 10,
            seq_id: 1,
            pos: 100,
            strand: true,
        },
        MinimizerInfo {
            hash: 20,
            seq_id: 1,
            pos: 200,
            strand: true,
        },
        MinimizerInfo {
            hash: 30,
            seq_id: 1,
            pos: 300,
            strand: true,
        },
        MinimizerInfo {
            hash: 40,
            seq_id: 1,
            pos: 400,
            strand: true,
        },
    ];

    // Seq 2: 10 -> 20 -> 30 -> 50
    let seq2 = vec![
        MinimizerInfo {
            hash: 10,
            seq_id: 2,
            pos: 100,
            strand: true,
        },
        MinimizerInfo {
            hash: 20,
            seq_id: 2,
            pos: 200,
            strand: true,
        },
        MinimizerInfo {
            hash: 30,
            seq_id: 2,
            pos: 300,
            strand: true,
        },
        MinimizerInfo {
            hash: 50,
            seq_id: 2,
            pos: 500,
            strand: true,
        },
    ];

    graph.add_minimizers(&seq1, 100000);
    graph.add_minimizers(&seq2, 100000);

    // Initial check: edges count
    // 10->20: 2 edges
    // 20->30: 2 edges
    // 30->40: 1 edge
    // 30->50: 1 edge
    // Total edges: 6
    assert_eq!(graph.graph.edge_count(), 6);

    // Prune edges with weight < 2
    graph.prune_low_weight_edges(2);

    // Remaining edges: 10->20 (2), 20->30 (2). Total 4.
    assert_eq!(graph.graph.edge_count(), 4);

    // Find paths
    let paths = graph.get_linear_paths();

    // Should find one path: 10 -> 20 -> 30
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0], vec![10, 20, 30]);
}

#[test]
fn test_synteny_graph_cycle() {
    let mut graph = SyntenyGraph::new();

    // Seq 1: 10 -> 20 -> 10 (cycle)
    let seq1 = vec![
        MinimizerInfo {
            hash: 10,
            seq_id: 1,
            pos: 100,
            strand: true,
        },
        MinimizerInfo {
            hash: 20,
            seq_id: 1,
            pos: 200,
            strand: true,
        },
        MinimizerInfo {
            hash: 10,
            seq_id: 1,
            pos: 300,
            strand: true,
        },
    ];

    graph.add_minimizers(&seq1, 100000);

    // Prune with weight 1 (keep everything)
    graph.prune_low_weight_edges(1);

    let paths = graph.get_linear_paths();

    // Cycle logic might pick any start node if pure cycle
    // 10 -> 20 -> 10
    // Nodes: 10, 20.
    // 10: Out(20), In(20)
    // 20: Out(10), In(10)
    // All degree 1.
    // Logic should handle it.

    assert_eq!(paths.len(), 1);
    // Path could be 10->20->10 or 20->10->20 depending on start
    // My implementation breaks at visited, so it should be length 3 (node count) if we include closing node?
    // Wait, traverse_path loop:
    // push(curr)
    // next = ...
    // curr = next
    // if visited.contains(curr) break

    // Trace:
    // Start 10.
    // Path: [10]
    // Next: 20. Visited? No.
    // Curr = 20.
    // Loop.
    // Path: [10, 20]
    // Next: 10. Visited? Yes.
    // Break.

    // Result: [10, 20].
    // Note: The edge 20->10 exists, but we stop when we see 10 again.
    // So we get [10, 20].

    assert_eq!(paths[0].len(), 2);
    assert!(paths[0].contains(&10));
    assert!(paths[0].contains(&20));
}

#[test]
fn test_synteny_finder_run() -> anyhow::Result<()> {
    use crate::libs::synteny::algo::SyntenyFinder;
    let finder = SyntenyFinder::new(5, vec![5], 2, 100, 0, 100000, false);
    let seq1 = b"ACGTACGTACGTACGTACGT";
    let seq2 = b"ACGTACGTACGTACGTACGT";
    let mut blocks = Vec::new();

    finder.run(
        |emit| {
            emit("seq1", seq1);
            emit("seq2", seq2);
            Ok(())
        },
        |w, block| {
            blocks.push((w, block.clone()));
        },
    )?;

    // Should find synteny
    // With k=5, w=5:
    // Seq len 20.
    // Should produce minimizers.
    // Since identical, should form a perfect path.
    // min_weight=2 means we need at least 2 genomes (we have 2).

    assert!(!blocks.is_empty());

    // Check block content
    // We expect at least one block
    let (_w, block) = &blocks[0];
    assert!(block.ranges.contains_key(&1)); // seq1
    assert!(block.ranges.contains_key(&2)); // seq2

    Ok(())
}
