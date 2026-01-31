use clap::*;
use hnsm::libs::synteny::io::{read_blocks, write_blocks, Block, Segment};

pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merge fragmented synteny blocks")
        .after_help(
            "Default parameters based on divergence (-d):\n  \
             < 1.0%:   --chain-gap 10000\n  \
             1.0-10.0%: --chain-gap 100000\n  \
             > 10.0%:  --chain-gap 1000000\n\
             \n\
             If --divergence is not specified, --chain-gap defaults to 100000.",
        )
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
                .help("Output filename. [stdout] for screen")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("divergence")
                .short('d')
                .long("divergence")
                .value_parser(value_parser!(f64))
                .help("Approx. maximum percent sequence divergence (e.g., 1.0 for 1%)"),
        )
        .arg(
            Arg::new("chain_gap")
                .long("chain-gap")
                .value_parser(value_parser!(u64))
                .help("Maximum gap size allowed for merging [default: 100000]"),
        )
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .short('v')
                .action(clap::ArgAction::SetTrue)
                .help("Verbose output"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();

    // Determine default chain_gap based on divergence
    // ntSynt:
    // < 1%: 10,000
    // 1-10%: 100,000
    // > 10%: 1,000,000
    // Default if no divergence: 100,000 (middle ground)
    let divergence = matches.get_one::<f64>("divergence");
    let default_gap = if let Some(d) = divergence {
        if *d < 1.0 {
            10_000
        } else if *d <= 10.0 {
            100_000
        } else {
            1_000_000
        }
    } else {
        100_000
    };

    let chain_gap = *matches.get_one::<u64>("chain_gap").unwrap_or(&default_gap);

    if matches.get_flag("verbose") {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Info)
            .init();
    } else {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Warn)
            .init();
    }

    log::info!("Loading blocks from {}...", infile);
    let mut blocks = read_blocks(infile)?;
    log::info!("Loaded {} blocks.", blocks.len());

    // Normalize block ranges: sort by seq_name to ensure consistent order for N-way comparison
    for block in &mut blocks {
        block.ranges.sort_by(|a, b| a.seq_name.cmp(&b.seq_name));
    }

    let merged = merge_blocks(blocks, chain_gap);
    log::info!("Merged into {} blocks.", merged.len());

    write_blocks(&merged, outfile)?;
    // log::info!("Written to {}", outfile);

    Ok(())
}

fn merge_blocks(mut blocks: Vec<Block>, chain_gap: u64) -> Vec<Block> {
    // Sort primarily by first range (Reference)
    // We assume the first range in each block corresponds to the same genome/role
    // or at least that blocks are "comparable" via their first range.
    blocks.sort_by(|a, b| {
        if a.ranges.is_empty() || b.ranges.is_empty() {
            return std::cmp::Ordering::Equal;
        }
        let r1 = &a.ranges[0];
        let r2 = &b.ranges[0];

        r1.seq_name.cmp(&r2.seq_name).then(r1.start.cmp(&r2.start))
    });

    let mut merged = Vec::new();
    if blocks.is_empty() {
        return merged;
    }

    let mut current = blocks[0].clone();

    for next in blocks.into_iter().skip(1) {
        if can_merge(&current, &next, chain_gap) {
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

fn can_merge(b1: &Block, b2: &Block, chain_gap: u64) -> bool {
    // 1. Must have same number of ranges
    if b1.ranges.len() != b2.ranges.len() {
        return false;
    }

    // 2. Iterate through all ranges and check consistency
    // Ranges are now sorted by seq_name, so we can safely compare by index.

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
            // Forward strand case:
            // The "next" block (r2) should physically follow the "current" block (r1).
            // Genomic coordinates: r1.start < r1.end <= r2.start < r2.end
            if r2.start < r1.end {
                // Overlap or wrong order
                if r2.start < r1.start {
                    return false;
                } // Definitely wrong order (r2 starts before r1)
                0 // Overlap treated as 0 gap (allow merging overlapping blocks)
            } else {
                r2.start - r1.end
            }
        } else {
            // Reverse strand case:
            // The "next" block (r2) should logically follow r1 in the query's traversal order.
            // But since it's on the negative strand, "logical next" means "physical previous" (upstream).
            // Genomic coordinates: r2.start < r2.end <= r1.start < r1.end
            // Example:
            // Block 1 (r1): 400-500 (-)  <- Physically downstream
            // Block 2 (r2): 200-300 (-)  <- Physically upstream (logical next)
            if r2.end > r1.start {
                // Wrong order: r2 is physically after r1, which violates the negative strand collinearity
                return false;
            }
            r1.start - r2.end
        };

        if gap > chain_gap {
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
