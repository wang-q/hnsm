use clap::{Arg, ArgAction, Command};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use hnsm::libs::synteny::io::{read_blocks, Block};

pub fn make_subcommand() -> Command {
    Command::new("circle")
        .about("Visualize synteny blocks as circular Circos-style SVG")
        .after_help(
            r#"
EXAMPLES:
    # Basic usage
    hnsm synt circle blocks.tsv -o plot.svg

    # With size files for accurate chromosome lengths
    hnsm synt circle blocks.tsv genome1.size.tsv genome2.size.tsv -o plot.svg

    # Custom size
    hnsm synt circle blocks.tsv -o plot.svg --width 800
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
                .default_value("800")
                .value_parser(clap::value_parser!(u32))
                .help("Width/Height of the SVG canvas (square)"),
        )
        .arg(
            Arg::new("track_width")
                .long("track-width")
                .default_value("20")
                .value_parser(clap::value_parser!(f64))
                .help("Width of the chromosome track ring"),
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
    let track_width = *matches.get_one::<f64>("track_width").unwrap();
    let no_label = matches.get_flag("no_label");

    // 1. Read blocks
    let blocks = read_blocks(infile)?;

    // 2. Parse size files (if any) or infer from blocks
    let mut chrom_lengths: HashMap<String, u64> = HashMap::new();
    let mut genome_chroms: HashMap<String, Vec<String>> = HashMap::new(); // genome -> [chrom_names]
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

            let file = File::open(size_path).expect("Could not open size file");
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line?;
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    let name = parts[0].to_string();
                    let len: u64 = parts[1].parse().unwrap();
                    chrom_lengths.insert(name.clone(), len);
                    entry.push(name.clone());
                    seq_to_genome.insert(name.clone(), genome_name.clone());
                }
            }
        }
    } else {
        // Infer from blocks
        for block in &blocks {
            for range in &block.ranges {
                let current = chrom_lengths.get(&range.seq_name).cloned().unwrap_or(0);
                if range.end > current {
                    chrom_lengths.insert(range.seq_name.clone(), range.end);
                }
                
                // Try to guess genome name (e.g. "Genome.Chrom")
                if !seq_to_genome.contains_key(&range.seq_name) {
                    let g_name = if let Some((g, _)) = range.seq_name.split_once('.') {
                        g.to_string()
                    } else {
                        "Genome".to_string()
                    };
                    seq_to_genome.insert(range.seq_name.clone(), g_name.clone());
                    
                    if !genome_chroms.contains_key(&g_name) {
                        genome_order.push(g_name.clone());
                    }
                    genome_chroms.entry(g_name).or_default().push(range.seq_name.clone());
                }
            }
        }
        // Deduplicate chroms in genome_chroms
        for chroms in genome_chroms.values_mut() {
            chroms.sort();
            chroms.dedup();
        }
    }
    
    // Auto-prefix logic for size files
    if !size_files.is_empty() {
        // Collect all seq_names from blocks
        let mut block_seqs = std::collections::HashSet::new();
        for block in &blocks {
            for range in &block.ranges {
                block_seqs.insert(range.seq_name.clone());
            }
        }
        
        // Check if we need to prepend genome name
        let mut new_chrom_lengths = HashMap::new();
        let mut new_seq_to_genome = HashMap::new();
        let mut new_genome_chroms = HashMap::new();

        for g_name in &genome_order {
            let chroms = genome_chroms.get(g_name).unwrap();
            let mut new_chroms_list = Vec::new();
            
            for chrom in chroms {
                let len = *chrom_lengths.get(chrom).unwrap();
                let full_name = format!("{}.{}", g_name, chrom);
                
                if block_seqs.contains(&full_name) && !block_seqs.contains(chrom) {
                     // Need prefix
                     new_chrom_lengths.insert(full_name.clone(), len);
                     new_seq_to_genome.insert(full_name.clone(), g_name.clone());
                     new_chroms_list.push(full_name);
                } else {
                     // No prefix or already has it
                     new_chrom_lengths.insert(chrom.clone(), len);
                     new_seq_to_genome.insert(chrom.clone(), g_name.clone());
                     new_chroms_list.push(chrom.clone());
                }
            }
            new_genome_chroms.insert(g_name.clone(), new_chroms_list);
        }
        
        chrom_lengths = new_chrom_lengths;
        seq_to_genome = new_seq_to_genome;
        genome_chroms = new_genome_chroms;
    }

    // 3. Layout
    // Calculate total length to map to 360 degrees (minus gaps)
    let total_len: u64 = chrom_lengths.values().sum();
    
    let gap_degrees = 2.0; // Gap between chromosomes in degrees
    let total_chroms: usize = genome_chroms.values().map(|v| v.len()).sum();
    let total_gap_degrees = gap_degrees * total_chroms as f64;
    
    let available_degrees = 360.0 - total_gap_degrees;
    let scale_angle = available_degrees.to_radians() / total_len as f64; // radians per bp

    let mut chrom_angles: HashMap<String, (f64, f64)> = HashMap::new();
    let mut current_angle = -PI / 2.0; // Start at 12 o'clock (-90 degrees)
    
    let gap_rad = gap_degrees.to_radians();

    for g_name in &genome_order {
        let chroms = genome_chroms.get(g_name).unwrap();
        for chrom in chroms {
            let len = *chrom_lengths.get(chrom).unwrap_or(&1000);
            let start_a = current_angle;
            let sweep = len as f64 * scale_angle;
            let end_a = start_a + sweep;
            
            chrom_angles.insert(chrom.clone(), (start_a, end_a));
            
            current_angle = end_a + gap_rad;
        }
    }

    let layout = Layout {
        width,
        height: width,
        cx: width / 2.0,
        cy: width / 2.0,
        radius: width / 2.0 - 50.0, // Margin
        track_width,
        chrom_angles,
        seq_to_genome,
        scale_angle,
    };

    // 4. Output
    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(std::io::stdout())
    } else {
        Box::new(std::io::BufWriter::new(File::create(outfile)?))
    };

    writeln!(writer, r#"<svg width="{}" height="{}" viewBox="0 0 {} {}" xmlns="http://www.w3.org/2000/svg">"#, width, width, width, width)?;
    writeln!(writer, r#"<style>text {{ font-family: sans-serif; }}</style>"#)?;
    
    // Draw tracks
    for g_name in &genome_order {
        let chroms = genome_chroms.get(g_name).unwrap();
        for chrom in chroms {
            draw_chrom_track(&mut writer, &layout, chrom, no_label)?;
        }
    }
    
    // Draw ribbons
    draw_ribbons(&mut writer, &layout, &blocks)?;

    writeln!(writer, "</svg>")?;

    Ok(())
}

struct Layout {
    #[allow(dead_code)]
    width: f64,
    #[allow(dead_code)]
    height: f64,
    cx: f64,
    cy: f64,
    radius: f64,
    track_width: f64,
    chrom_angles: HashMap<String, (f64, f64)>,
    seq_to_genome: HashMap<String, String>,
    scale_angle: f64,
}

fn polar_to_cartesian(cx: f64, cy: f64, radius: f64, angle: f64) -> (f64, f64) {
    (
        cx + radius * angle.cos(),
        cy + radius * angle.sin(),
    )
}

fn draw_chrom_track(writer: &mut impl Write, layout: &Layout, chrom: &str, no_label: bool) -> std::io::Result<()> {
    if let Some((start_a, end_a)) = layout.chrom_angles.get(chrom) {
        // SVG Arc: A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        
        let r_outer = layout.radius;
        let r_inner = layout.radius - layout.track_width;

        let (x1_out, y1_out) = polar_to_cartesian(layout.cx, layout.cy, r_outer, *start_a);
        let (x2_out, y2_out) = polar_to_cartesian(layout.cx, layout.cy, r_outer, *end_a);
        let (x1_in, y1_in) = polar_to_cartesian(layout.cx, layout.cy, r_inner, *start_a);
        let (x2_in, y2_in) = polar_to_cartesian(layout.cx, layout.cy, r_inner, *end_a);

        let large_arc = if end_a - start_a > PI { 1 } else { 0 };

        // Path: Move to start_outer -> Arc to end_outer -> Line to end_inner -> Arc to start_inner -> Close
        writeln!(
            writer,
            r##"<path d="M {:.1} {:.1} A {:.1} {:.1} 0 {} 1 {:.1} {:.1} L {:.1} {:.1} A {:.1} {:.1} 0 {} 0 {:.1} {:.1} Z" fill="#ddd" stroke="#999" />"##,
            x1_out, y1_out,
            r_outer, r_outer,
            large_arc,
            x2_out, y2_out,
            x2_in, y2_in,
            r_inner, r_inner,
            large_arc,
            x1_in, y1_in
        )?;

        // Label
        if !no_label {
            let mid_angle = (start_a + end_a) / 2.0;
            let label_r = layout.radius + 15.0;
            let (lx, ly) = polar_to_cartesian(layout.cx, layout.cy, label_r, mid_angle);
            
            // Rotation for readability
            let mut rotation = mid_angle * 180.0 / PI;
            if rotation > 90.0 && rotation < 270.0 {
                rotation += 180.0;
            }

            // Genome Name? Or Chrom name?
            let g_name = layout.seq_to_genome.get(chrom).unwrap();
            let display_name = chrom.strip_prefix(&format!("{}.", g_name)).unwrap_or(chrom);

            writeln!(
                writer,
                r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" dominant-baseline="middle" transform="rotate({:.1}, {:.1}, {:.1})" font-size="10" fill="#333">{}</text>"##,
                lx, ly, rotation, lx, ly, display_name
            )?;
        }
    }
    Ok(())
}

fn draw_ribbons(writer: &mut impl Write, layout: &Layout, blocks: &[Block]) -> std::io::Result<()> {
    let colors = [
        "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
    ];

    let r = layout.radius - layout.track_width - 2.0; // Slightly inside track

    for (block_idx, block) in blocks.iter().enumerate() {
        if block.ranges.len() < 2 {
            continue;
        }

        let r1 = &block.ranges[0];
        let r2 = &block.ranges[1];

        if !layout.chrom_angles.contains_key(&r1.seq_name) || !layout.chrom_angles.contains_key(&r2.seq_name) {
            continue;
        }

        let (start1_base, _) = layout.chrom_angles.get(&r1.seq_name).unwrap();
        let (start2_base, _) = layout.chrom_angles.get(&r2.seq_name).unwrap();

        let a1_start = start1_base + r1.start as f64 * layout.scale_angle;
        let a1_end = start1_base + r1.end as f64 * layout.scale_angle;

        // Handle strand
        let (a2_start, a2_end) = if r1.strand == r2.strand {
             // Direct: Connect End1->End2, Start2->Start1
             // Arc2 goes End2 -> Start2 (Reverse)
             (
                start2_base + r2.end as f64 * layout.scale_angle,
                start2_base + r2.start as f64 * layout.scale_angle,
             )
        } else {
             // Inverted: Connect End1->Start2, End2->Start1
             // Arc2 goes Start2 -> End2 (Forward)
             (
                start2_base + r2.start as f64 * layout.scale_angle,
                start2_base + r2.end as f64 * layout.scale_angle,
             )
        };
        
        // Points on the inner rim
        let (x1s, y1s) = polar_to_cartesian(layout.cx, layout.cy, r, a1_start);
        let (x1e, y1e) = polar_to_cartesian(layout.cx, layout.cy, r, a1_end);
        let (x2s, y2s) = polar_to_cartesian(layout.cx, layout.cy, r, a2_start);
        let (x2e, y2e) = polar_to_cartesian(layout.cx, layout.cy, r, a2_end);

        let color = colors[block_idx % colors.len()];
        let opacity = 0.5;

        // Arc flags
        let large_arc1 = if (a1_end - a1_start).abs() > PI { 1 } else { 0 };
        let sweep1 = if a1_end > a1_start { 1 } else { 0 }; // Usually 1 if we go clockwise
        
        let large_arc2 = if (a2_end - a2_start).abs() > PI { 1 } else { 0 };
        let sweep2 = if a2_end > a2_start { 1 } else { 0 }; 

        // For curves, using cubic bezier with control points towards center
        let cp_r = 0.0; // Pull to center
        let (cpx, cpy) = polar_to_cartesian(layout.cx, layout.cy, cp_r, 0.0); // Center

        writeln!(
            writer,
            r##"<path d="M {:.1} {:.1} A {:.1} {:.1} 0 {} {} {:.1} {:.1} Q {:.1} {:.1}, {:.1} {:.1} A {:.1} {:.1} 0 {} {} {:.1} {:.1} Q {:.1} {:.1}, {:.1} {:.1} Z" fill="{}" fill-opacity="{}" stroke="none" />"##,
            x1s, y1s,
            r, r, large_arc1, sweep1, x1e, y1e,
            cpx, cpy, x2s, y2s,
            r, r, large_arc2, sweep2, x2e, y2e,
            cpx, cpy, x1s, y1s,
            color, opacity
        )?;
    }
    Ok(())
}
