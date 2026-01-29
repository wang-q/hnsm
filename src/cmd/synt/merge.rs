use clap::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use anyhow::Context;

#[derive(Debug, Clone)]
struct Anchor {
    q_name: String,
    q_start: u64,
    q_end: u64,
    t_name: String,
    t_start: u64,
    t_end: u64,
    strand: char, // '+' or '-'
    score: f64,
}

impl Anchor {
    fn from_line(line: &str) -> Option<Self> {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() < 7 {
            return None;
        }

        let q_name = fields[0].to_string();
        let q_start = fields[1].parse().ok()?;
        let q_end = fields[2].parse().ok()?;
        let t_name = fields[3].to_string();
        let t_start = fields[4].parse().ok()?;
        let t_end = fields[5].parse().ok()?;
        let strand = fields[6].chars().next()?;
        let score = if fields.len() > 7 {
            fields[7].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        Some(Anchor {
            q_name,
            q_start,
            q_end,
            t_name,
            t_start,
            t_end,
            strand,
            score,
        })
    }

}

pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merge fragmented synteny blocks")
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input synteny file (TSV: q_name, q_start, q_end, t_name, t_start, t_end, strand, score)"),
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

    eprintln!("Loading anchors from {}...", infile);
    let anchors = load_anchors(infile)?;
    eprintln!("Loaded {} anchors.", anchors.len());

    let merged = merge_anchors(anchors, max_gap);
    eprintln!("Merged into {} blocks.", merged.len());

    write_anchors(&merged, outfile)?;
    eprintln!("Written to {}", outfile);

    Ok(())
}

fn load_anchors<P: AsRef<Path>>(path: P) -> anyhow::Result<Vec<Anchor>> {
    let file = File::open(path).context("Failed to open input file")?;
    let reader = BufReader::new(file);
    let mut anchors = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some(anchor) = Anchor::from_line(&line) {
            anchors.push(anchor);
        }
    }
    Ok(anchors)
}

fn write_anchors<P: AsRef<Path>>(anchors: &[Anchor], path: P) -> anyhow::Result<()> {
    let file = File::create(path).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "#q_name\tq_start\tq_end\tt_name\tt_start\tt_end\tstrand\tscore")?;
    for a in anchors {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            a.q_name, a.q_start, a.q_end, a.t_name, a.t_start, a.t_end, a.strand, a.score
        )?;
    }
    Ok(())
}

fn merge_anchors(mut anchors: Vec<Anchor>, max_gap: u64) -> Vec<Anchor> {
    // Sort primarily by query position
    anchors.sort_by(|a, b| {
        a.q_name
            .cmp(&b.q_name)
            .then(a.t_name.cmp(&b.t_name))
            .then(a.strand.cmp(&b.strand))
            .then(a.q_start.cmp(&b.q_start))
    });

    let mut merged = Vec::new();
    if anchors.is_empty() {
        return merged;
    }

    let mut current = anchors[0].clone();

    for next in anchors.into_iter().skip(1) {
        if can_merge(&current, &next, max_gap) {
            current = merge_two(&current, &next);
        } else {
            merged.push(current);
            current = next;
        }
    }
    merged.push(current);

    merged
}

fn can_merge(prev: &Anchor, next: &Anchor, max_gap: u64) -> bool {
    if prev.q_name != next.q_name || prev.t_name != next.t_name || prev.strand != next.strand {
        return false;
    }

    // Check query gap (next is guaranteed >= prev.q_start by sort)
    // Overlap is allowed (gap < 0), effectively 0 distance
    let q_gap = next.q_start.saturating_sub(prev.q_end);
    if q_gap > max_gap {
        return false;
    }

    // Check target gap based on strand
    if prev.strand == '+' {
        // Monotonic increasing on target
        if next.t_start < prev.t_start {
             // Order violation on target side for + strand
             return false; 
        }
        let t_gap = next.t_start.saturating_sub(prev.t_end);
        if t_gap > max_gap {
            return false;
        }
    } else {
        // Monotonic decreasing on target (for increasing query)
        // prev: Q[100-200], T[500-600] (-)
        // next: Q[300-400], T[300-400] (-) -> Valid
        // next.t_end should be < prev.t_start approximately
        
        if next.t_end > prev.t_end {
            // Order violation: next block on target is "after" previous block,
            // but for '-' strand it should be "before" (smaller coords)
            return false;
        }

        let t_gap = prev.t_start.saturating_sub(next.t_end);
        if t_gap > max_gap {
            return false;
        }
    }

    true
}

fn merge_two(prev: &Anchor, next: &Anchor) -> Anchor {
    let mut new_anchor = prev.clone();
    
    // Query always extends min to max
    new_anchor.q_start = prev.q_start.min(next.q_start);
    new_anchor.q_end = prev.q_end.max(next.q_end);

    // Target depends on strand, but since we store start/end as absolute coords:
    // For +, start is min, end is max
    // For -, start is min, end is max (just the span covers both)
    // Wait, the block definition is usually start < end.
    // The "span" of the merged block is the bounding box.
    new_anchor.t_start = prev.t_start.min(next.t_start);
    new_anchor.t_end = prev.t_end.max(next.t_end);

    new_anchor.score += next.score;
    
    new_anchor
}
