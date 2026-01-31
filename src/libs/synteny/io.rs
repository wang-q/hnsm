use anyhow::Context;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Segment {
    pub seq_name: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub score: f64,
}

impl Segment {
    pub fn from_line(line: &str) -> Option<(usize, Self)> {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() < 3 {
            return None;
        }
        // # Block_ID Range Score ...
        let block_id = fields[0].parse().ok()?;
        let range_str = fields[1];
        let score = fields[2].parse().unwrap_or(0.0);

        // Parse range_str: SeqName(Strand):Start-End
        let (left, right) = range_str.rsplit_once(':')?;
        let (start_str, end_str) = right.split_once('-')?;
        let start = start_str.parse().ok()?;
        let end = end_str.parse().ok()?;

        if !left.ends_with(')') {
            return None;
        }
        let left_len = left.len();
        if left_len < 3 {
            return None;
        }
        if left.chars().nth(left_len - 3)? != '(' {
            return None;
        }
        let strand = left.chars().nth(left_len - 2)?;
        let seq_name = &left[..left_len - 3];

        Some((
            block_id,
            Segment {
                seq_name: seq_name.to_string(),
                start,
                end,
                strand,
                score,
            },
        ))
    }
}

#[derive(Debug, Clone)]
pub struct Block {
    pub id: usize,
    pub ranges: Vec<Segment>,
}

pub fn read_blocks<P: AsRef<Path>>(path: P) -> anyhow::Result<Vec<Block>> {
    let file = File::open(&path).context("Failed to open input file")?;
    let reader = BufReader::new(file);
    read_blocks_from_reader(reader)
}

pub fn read_blocks_from_reader<R: BufRead>(reader: R) -> anyhow::Result<Vec<Block>> {
    let mut blocks = Vec::new();
    let mut current_block_id = None;
    let mut current_ranges = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        if let Some((block_id, range)) = Segment::from_line(&line) {
            match current_block_id {
                Some(id) if id == block_id => {
                    current_ranges.push(range);
                }
                Some(_) => {
                    // Block ID changed, save previous block
                    if !current_ranges.is_empty() {
                        blocks.push(Block {
                            id: blocks.len() + 1, // Re-index locally? Or keep original ID?
                            // merge.rs re-indexes locally: id: blocks.len() + 1
                            // But usually preserving IDs is better if we want to trace back.
                            // However, if input IDs are non-contiguous or unordered, re-indexing ensures consistency in Vec.
                            // Let's keep re-indexing behavior from merge.rs for now.
                            ranges: current_ranges,
                        });
                    }
                    current_ranges = Vec::new();
                    current_block_id = Some(block_id);
                    current_ranges.push(range);
                }
                None => {
                    current_block_id = Some(block_id);
                    current_ranges.push(range);
                }
            }
        }
    }

    // Process last block
    if !current_ranges.is_empty() {
        blocks.push(Block {
            id: blocks.len() + 1,
            ranges: current_ranges,
        });
    }

    Ok(blocks)
}

pub fn write_blocks(blocks: &[Block], path: &str) -> anyhow::Result<()> {
    let mut writer = intspan::writer(path);

    writeln!(writer, "# Block_ID\tRange\tScore")?;
    for block in blocks {
        let flip = if let Some(first) = block.ranges.first() {
            first.strand == '-'
        } else {
            false
        };

        for range in &block.ranges {
            let strand = if flip {
                match range.strand {
                    '+' => '-',
                    '-' => '+',
                    _ => range.strand,
                }
            } else {
                range.strand
            };

            writeln!(
                writer,
                "{}\t{}({}):{}-{}\t{:.1}",
                block.id, range.seq_name, strand, range.start, range.end, range.score
            )?;
        }
    }
    Ok(())
}
