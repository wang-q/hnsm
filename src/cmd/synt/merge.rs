use clap::*;
use hnsm::libs::synteny::io::{read_blocks, write_blocks, Block, Segment};

pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merge fragmented synteny blocks")
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input synteny file (Format: hnsm Block TSV)"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .default_value("merged.tsv")
                .help("Output filename"),
        )
        .arg(
            Arg::new("max_gap")
                .long("max-gap")
                .default_value("100000")
                .value_parser(value_parser!(u64))
                .help("Maximum gap size allowed for merging"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();
    let max_gap = *matches.get_one::<u64>("max_gap").unwrap();

    eprintln!("Loading blocks from {}...", infile);
    let blocks = read_blocks(infile)?;
    eprintln!("Loaded {} blocks.", blocks.len());

    let merged = merge_blocks(blocks, max_gap);
    eprintln!("Merged into {} blocks.", merged.len());

    write_blocks(&merged, outfile)?;
    eprintln!("Written to {}", outfile);

    Ok(())
}

fn merge_blocks(mut blocks: Vec<Block>, max_gap: u64) -> Vec<Block> {
    // Sort primarily by first range (Reference)
    // We assume the first range in each block corresponds to the same genome/role
    // or at least that blocks are "comparable" via their first range.
    blocks.sort_by(|a, b| {
        if a.ranges.is_empty() || b.ranges.is_empty() {
            return std::cmp::Ordering::Equal;
        }
        let r1 = &a.ranges[0];
        let r2 = &b.ranges[0];
        
        r1.seq_name.cmp(&r2.seq_name)
            .then(r1.start.cmp(&r2.start))
    });

    let mut merged = Vec::new();
    if blocks.is_empty() {
        return merged;
    }

    let mut current = blocks[0].clone();

    for next in blocks.into_iter().skip(1) {
        if can_merge(&current, &next, max_gap) {
            current = merge_two(&current, &next);
        } else {
            // Assign new ID to finalized block
            current.id = merged.len() + 1;
            merged.push(current);
            current = next;
        }
    }
    current.id = merged.len() + 1;
    merged.push(current);

    merged
}

fn can_merge(b1: &Block, b2: &Block, max_gap: u64) -> bool {
    // 1. Must have same number of ranges
    if b1.ranges.len() != b2.ranges.len() {
        return false;
    }

    // 2. Iterate through all ranges and check consistency
    // We assume ranges in a block are ordered (e.g., Query, Target1, Target2...)
    // If not, we should try to match by seq_name.
    // For safety, let's match by index first, but verify seq_name.
    
    for (i, r1) in b1.ranges.iter().enumerate() {
        let r2 = &b2.ranges[i];
        
        // Genome mismatch
        if r1.seq_name != r2.seq_name {
            return false;
        }

        // Strand mismatch
        if r1.strand != r2.strand {
            return false;
        }

        // Distance check
        let gap = if r1.strand == '+' {
            // Forward: b2 should be after b1
            if r2.start < r1.end {
                // Overlap or wrong order
                // Allow small overlap? Usually chaining requires strict order.
                // But simplified: check signed distance.
                // If r2.start < r1.end, it's an overlap.
                // Overlap is usually fine for merging if it's small?
                // Or maybe strictly >?
                // Let's allow overlap if it's not "contained".
                // Actually, if we are merging "blocks", they should be sequential.
                // Let's use absolute distance of (r2.start as i64 - r1.end as i64).
                
                if r2.start < r1.start { return false; } // Definitely wrong order
                0 // Overlap treated as 0 gap
            } else {
                r2.start - r1.end
            }
        } else {
            // Reverse: b2 should be "before" b1 in coordinates (logical after)
            // Wait, if strand is -, it means the segment is inverted relative to the "positive" reference.
            // But the coordinates are still genomic (start < end).
            // If we are traversing the query in + direction, and we see a synteny block on - strand.
            // Block 1: Q:100-200, T:500-400 (stored as 400-500, -)
            // Block 2: Q:300-400, T:300-200 (stored as 200-300, -)
            
            // T1: 400-500 (-). T2: 200-300 (-).
            // Logical flow on T: 500 -> 400 -> 300 -> 200.
            // So T1.start (400) > T2.end (300).
            // Gap = T1.start - T2.end.
            
            if r2.end > r1.start {
                 // Wrong order for negative strand traversal
                 // T1 (400-500), T2 (600-700) -> T1..T2 is + direction.
                 // But we want - direction.
                 return false;
            }
            r1.start - r2.end
        };

        if gap > max_gap {
            return false;
        }
    }

    true
}

fn merge_two(b1: &Block, b2: &Block) -> Block {
    let mut new_ranges = Vec::new();
    
    for (i, r1) in b1.ranges.iter().enumerate() {
        let r2 = &b2.ranges[i];
        
        // Merge range: min start, max end
        // This is safe because we checked collinearity and overlap/gap.
        new_ranges.push(Segment {
            seq_name: r1.seq_name.clone(),
            start: r1.start.min(r2.start),
            end: r1.end.max(r2.end),
            strand: r1.strand,
            score: r1.score + r2.score, // Sum scores? Or max? Usually sum of counts.
        });
    }

    Block {
        id: b1.id, // ID will be reassigned later
        ranges: new_ranges,
    }
}
