use clap::{Arg, ArgAction, Command};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use hnsm::libs::synteny::io::{read_blocks, Block};

pub fn make_subcommand() -> Command {
    Command::new("view")
        .about("Visualize synteny blocks as SVG")
        .after_help(
            r#"
EXAMPLES:
    # Basic usage
    hnsm synt view blocks.tsv -o plot.svg

    # With size files for accurate chromosome lengths
    hnsm synt view blocks.tsv genome1.size.tsv genome2.size.tsv -o plot.svg

    # Custom size
    hnsm synt view blocks.tsv -o plot.svg --width 1200 --height 200
"#,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input synteny blocks file (.tsv)"),
        )
        .arg(
            Arg::new("size_files")
                .action(ArgAction::Append)
                .index(2)
                .help("Optional size files for chromosome lengths"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .help("Output filename. [stdout] for screen")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("width")
                .long("width")
                .default_value("1000")
                .value_parser(clap::value_parser!(u32))
                .help("Width of the SVG canvas"),
        )
        .arg(
            Arg::new("height")
                .long("height")
                .default_value("300")
                .value_parser(clap::value_parser!(u32))
                .help("Height of each track (distance between genomes)"),
        )
        .arg(
            Arg::new("no_label")
                .long("no-label")
                .action(ArgAction::SetTrue)
                .help("Do not draw labels"),
        )
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();
    let size_files: Vec<&String> = matches
        .get_many::<String>("size_files")
        .unwrap_or_default()
        .collect();
    let width = *matches.get_one::<u32>("width").unwrap() as f64;
    let track_height = *matches.get_one::<u32>("height").unwrap() as f64;
    let no_label = matches.get_flag("no_label");

    // 1. Read blocks
    let blocks = read_blocks(infile)?;

    // 2. Parse size files (if any) or infer from blocks
    let mut chrom_lengths: HashMap<String, u64> = HashMap::new();
    let mut genome_chroms: HashMap<String, Vec<String>> = HashMap::new(); // genome -> [chrom_names]

    // Load size files
    let mut seq_to_genome: HashMap<String, String> = HashMap::new();
    let mut genome_order: Vec<String> = Vec::new();

    if !size_files.is_empty() {
        for (_, size_path) in size_files.iter().enumerate() {
            let genome_name = std::path::Path::new(size_path)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .split('.')
                .next()
                .unwrap_or("unknown")
                .to_string();
            
            genome_order.push(genome_name.clone());
            let entry = genome_chroms.entry(genome_name.clone()).or_default();

            let file = File::open(size_path)?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line?;
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 2 {
                    let seq_name = fields[0].to_string();
                    let len: u64 = fields[1].parse()?;
                    chrom_lengths.insert(seq_name.clone(), len);
                    seq_to_genome.insert(seq_name.clone(), genome_name.clone());
                    entry.push(seq_name);
                }
            }
        }
    } else {
        // Infer mode
        // 1. Collect all seq names and max lengths from blocks
        let mut seq_max_len: HashMap<String, u64> = HashMap::new();

        for block in &blocks {
            for range in &block.ranges {
                let current_max = seq_max_len.entry(range.seq_name.clone()).or_insert(0);
                if range.end > *current_max {
                    *current_max = range.end;
                }
            }
        }

        chrom_lengths = seq_max_len.clone();

        // 2. Group sequences into "Genomes" by name prefix
        // e.g. "mg1655.NC_000913" -> Genome "mg1655"
        let all_seqs: Vec<String> = seq_max_len.keys().cloned().collect();
        
        for seq in all_seqs {
            let genome_name = seq.split('.').next().unwrap_or(&seq).to_string();
            seq_to_genome.insert(seq.clone(), genome_name.clone());
            genome_chroms.entry(genome_name).or_default().push(seq);
        }

        // Determine genome order (sorted alphabetically)
        let mut genomes: Vec<String> = genome_chroms.keys().cloned().collect();
        genomes.sort(); 
        
        genome_order = genomes;
    }

    // 3. Layout Tracks
    let margin_x = 100.0;
    let margin_y = 50.0;
    // let track_height = (height - 2.0 * margin_y) / (genome_order.len().max(1) as f64);
    let height = 2.0 * margin_y + track_height * (genome_order.len().max(1) as f64);
    
    // Find global max length to scale x
    let mut global_max_len = 0;
    for chroms in genome_chroms.values() {
        let total_len: u64 = chroms.iter().map(|c| *chrom_lengths.get(c).unwrap_or(&0)).sum();
        // Or if we stack them horizontally?
        // Usually ribbon plots stack chromosomes of one genome horizontally with gaps.
        if total_len > global_max_len {
            global_max_len = total_len;
        }
    }
    
    let scale_x = (width - 2.0 * margin_x) / (global_max_len as f64).max(1.0);
    
    // Calculate offsets for chromosomes on each track
    let mut chrom_offsets: HashMap<String, f64> = HashMap::new();
    let chrom_gap = global_max_len as f64 * 0.01; // 1% gap

    for g_name in &genome_order {
        let mut current_x = margin_x;
        if let Some(chroms) = genome_chroms.get(g_name) {
            // Sort chromosomes? Alphabetical for now if not specified
            let mut sorted_chroms = chroms.clone();
            sorted_chroms.sort();
            
            for chrom in sorted_chroms {
                chrom_offsets.insert(chrom.clone(), current_x);
                let len = *chrom_lengths.get(&chrom).unwrap_or(&0) as f64;
                current_x += len * scale_x + chrom_gap * scale_x;
            }
        }
    }

    let layout = Layout {
        width,
        height,
        margin_x,
        margin_y,
        track_height,
        scale_x,
        chrom_offsets,
        chrom_lengths,
        genome_order,
        genome_chroms,
        seq_to_genome,
        global_max_len,
    };

    // 4. Generate SVG
    let mut writer = intspan::writer(outfile);
    writeln!(writer, r#"<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">"#, layout.width, layout.height)?;
    writeln!(writer, r#"<style>text {{ font-family: sans-serif; font-size: 12px; }}</style>"#)?;
    writeln!(writer, r#"<rect width="100%" height="100%" fill="white" />"#)?;

    draw_tracks(&mut writer, &layout, no_label)?;
    draw_ribbons(&mut writer, &layout, &blocks)?;
    draw_scale_bar(&mut writer, &layout)?;

    writeln!(writer, "</svg>")?;

    Ok(())
}

struct Layout {
    width: f64,
    height: f64,
    margin_x: f64,
    margin_y: f64,
    track_height: f64,
    scale_x: f64,
    chrom_offsets: HashMap<String, f64>,
    chrom_lengths: HashMap<String, u64>,
    genome_order: Vec<String>,
    genome_chroms: HashMap<String, Vec<String>>,
    seq_to_genome: HashMap<String, String>,
    global_max_len: u64,
}

fn draw_tracks(
    writer: &mut impl Write,
    layout: &Layout,
    no_label: bool,
) -> std::io::Result<()> {
    for (i, g_name) in layout.genome_order.iter().enumerate() {
        let y = layout.margin_y + i as f64 * layout.track_height + layout.track_height / 2.0;
        
        // Draw Label
        if !no_label {
            writeln!(writer, r#"<text x="{}" y="{}" text-anchor="end" dominant-baseline="middle" font-weight="bold">{}</text>"#, layout.margin_x - 10.0, y + 5.0, g_name)?;
        }

        if let Some(chroms) = layout.genome_chroms.get(g_name) {
            let mut sorted_chroms = chroms.clone();
            sorted_chroms.sort();
            
            for chrom in sorted_chroms {
                let x = *layout.chrom_offsets.get(&chrom).unwrap();
                let len = *layout.chrom_lengths.get(&chrom).unwrap_or(&0) as f64;
                let w = len * layout.scale_x;
                
                // Draw Chromosome Bar
                writeln!(writer, r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="10" fill="#ddd" stroke="#999" rx="5" />"##, x, y, w)?;
                
                // Draw Chromosome Name
                if !no_label {
                    let display_name = chrom.strip_prefix(&format!("{}.", g_name)).unwrap_or(&chrom);
                    writeln!(writer, r##"<text x="{:.1}" y="{:.1}" dy="25" text-anchor="middle" font-size="10" fill="#666">{}</text>"##, x + w/2.0, y, display_name)?;
                }
            }
        }
    }
    Ok(())
}

fn draw_ribbons(
    writer: &mut impl Write,
    layout: &Layout,
    blocks: &[Block],
) -> std::io::Result<()> {
    let colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"];
    
    for (block_idx, block) in blocks.iter().enumerate() {
        if block.ranges.len() < 2 { continue; }
        
        // Assume pairwise for simplicity: Range 0 -> Range 1
        let r1 = &block.ranges[0];
        let r2 = &block.ranges[1];
        
        let g1 = layout.seq_to_genome.get(&r1.seq_name);
        let g2 = layout.seq_to_genome.get(&r2.seq_name);
        
        if g1.is_none() || g2.is_none() { continue; }
        let g1 = g1.unwrap();
        let g2 = g2.unwrap();
        
        // Get Y coordinates
        let y1_idx = layout.genome_order.iter().position(|g| g == g1).unwrap();
        let y2_idx = layout.genome_order.iter().position(|g| g == g2).unwrap();
        
        if y1_idx == y2_idx { continue; } // Intra-genome not supported well yet
        
        let y1 = layout.margin_y + y1_idx as f64 * layout.track_height + layout.track_height / 2.0 + 10.0; // Bottom of bar
        let y2 = layout.margin_y + y2_idx as f64 * layout.track_height + layout.track_height / 2.0; // Top of bar
        
        // Get X coordinates
        let x1_off = *layout.chrom_offsets.get(&r1.seq_name).unwrap();
        let x2_off = *layout.chrom_offsets.get(&r2.seq_name).unwrap();
        
        let x1_start = x1_off + r1.start as f64 * layout.scale_x;
        let x1_end = x1_off + r1.end as f64 * layout.scale_x;
        
        // Handle strand for r2
        let (x2_start, x2_end) = if r1.strand == r2.strand {
            (x2_off + r2.start as f64 * layout.scale_x, x2_off + r2.end as f64 * layout.scale_x)
        } else {
            // Invert visualization for inverted alignment?
            // Usually we draw a "twist".
            (x2_off + r2.end as f64 * layout.scale_x, x2_off + r2.start as f64 * layout.scale_x)
        };
        
        let color = colors[block_idx % colors.len()];
        let opacity = 0.5;
        
        // Bezier Path
        let h = (y2 - y1) / 2.0;
        
        writeln!(writer, r#"<path d="M {:.1} {:.1} C {:.1} {:.1}, {:.1} {:.1}, {:.1} {:.1} L {:.1} {:.1} C {:.1} {:.1}, {:.1} {:.1}, {:.1} {:.1} Z" fill="{}" fill-opacity="{}" stroke="none" />"#,
            x1_start, y1,
            x1_start, y1 + h, x2_start, y2 - h, x2_start, y2,
            x2_end, y2,
            x2_end, y2 - h, x1_end, y1 + h, x1_end, y1,
            color, opacity
        )?;
    }
    Ok(())
}

fn draw_scale_bar(
    writer: &mut impl Write,
    layout: &Layout,
) -> std::io::Result<()> {
    if layout.global_max_len > 0 {
        let target_width = (layout.width - 2.0 * layout.margin_x) / 5.0;
        let target_bp = target_width / layout.scale_x;
        let magnitude = 10f64.powf(target_bp.log10().floor());
        let normalized = target_bp / magnitude;
        let bar_bp = if normalized >= 5.0 {
            5.0 * magnitude
        } else if normalized >= 2.0 {
            2.0 * magnitude
        } else {
            1.0 * magnitude
        };
        let bar_width = bar_bp * layout.scale_x;
        let bar_label = if bar_bp >= 1_000_000.0 {
            format!("{:.0} Mb", bar_bp / 1_000_000.0)
        } else if bar_bp >= 1_000.0 {
            format!("{:.0} kb", bar_bp / 1_000.0)
        } else {
            format!("{:.0} bp", bar_bp)
        };

        let bar_x = layout.width - layout.margin_x - bar_width;
        let bar_y = layout.height - 20.0;
        
        writeln!(writer, r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="black" stroke-width="2" />"##, 
            bar_x, bar_y, bar_x + bar_width, bar_y)?;
        writeln!(writer, r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="black" stroke-width="2" />"##, 
            bar_x, bar_y - 5.0, bar_x, bar_y + 5.0)?; // Left tick
        writeln!(writer, r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="black" stroke-width="2" />"##, 
            bar_x + bar_width, bar_y - 5.0, bar_x + bar_width, bar_y + 5.0)?; // Right tick
        writeln!(writer, r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="12">{}</text>"##, 
            bar_x + bar_width / 2.0, bar_y - 10.0, bar_label)?;
    }
    Ok(())
}
