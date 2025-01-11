use clap::*;
use rust_xlsxwriter::*;
use std::cmp::max;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("xlsx")
        .about("List variations (substitutions/indels)")
        .after_help(
            r###"
* <infiles> are paths to block fasta files, .fas.gz is supported
    * infile == stdin means reading from STDIN

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("wrap")
                .long("wrap")
                .value_parser(value_parser!(u16))
                .num_args(1)
                .default_value("50")
                .help("Wrap length"),
        )
        .arg(
            Arg::new("indel")
                .long("indel")
                .action(ArgAction::SetTrue)
                .help("List indels"),
        )
        .arg(
            Arg::new("outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("There are outgroups at the end of each block"),
        )
        .arg(
            Arg::new("nosingle")
                .long("nosingle")
                .action(ArgAction::SetTrue)
                .help("Omit singleton SNPs and indels"),
        )
        .arg(
            Arg::new("nocomplex")
                .long("nocomplex")
                .action(ArgAction::SetTrue)
                .help("Omit complex SNPs and indels"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Minimal frequency"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Maximal frequency"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("variations.xlsx")
                .help("Output filename"),
        )
}

/// Enum to represent variations (substitutions or indels)
#[derive(Debug)]
enum Variation {
    Substitution(hnsm::Substitution),
    Indel(hnsm::Indel),
}

#[derive(Debug)]
struct Opt {
    sec_cursor: u32,     // Current section's starting row
    col_cursor: u16,     // Current column
    sec_height: u32,     // Height of each section
    max_name_len: usize, // Maximum name length
    wrap: u16,           // Wrap length
    color_loop: u32,     // Number of background colors for variations
    seq_count: u32,      // Number of sequences (excluding outgroup if applicable)
    is_outgroup: bool,   // Whether outgroups are present
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outfile = args.get_one::<String>("outfile").unwrap();

    let opt_wrap = *args.get_one::<u16>("wrap").unwrap();
    let is_outgroup = args.get_flag("outgroup");

    let is_indel = args.get_flag("indel");
    let is_nosingle = args.get_flag("nosingle");
    let is_nocomplex = args.get_flag("nocomplex");
    let opt_min = args.get_one::<f64>("min").cloned();
    let opt_max = args.get_one::<f64>("max").cloned();

    //----------------------------
    // Ops
    //----------------------------

    // Create workbook and worksheet objects
    let mut workbook = Workbook::new();
    let mut worksheet = workbook.add_worksheet();

    let format_of: BTreeMap<String, Format> = create_formats();
    // eprintln!("format_of = {:#?}", format_of.keys());

    let mut opt = Opt {
        sec_cursor: 1,
        col_cursor: 1,
        sec_height: 0,
        max_name_len: 1,
        wrap: opt_wrap,
        color_loop: 15,
        seq_count: 0,
        is_outgroup,
    };

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            // Get variations (substitutions and indels)
            let vars = get_vars(
                &seqs,
                is_outgroup,
                is_indel,
                is_nosingle,
                is_nocomplex,
                opt_min,
                opt_max,
            )?;

            opt.seq_count = seqs.len() as u32;
            opt.sec_height = opt.seq_count + 2; // 1 for pos, 1 for spacing
            opt.col_cursor = 1;

            // each section
            // Write names
            paint_name(&mut worksheet, &format_of.clone(), &mut opt, &block)?;

            if opt.is_outgroup {
                opt.seq_count -= 1;
            }

            // Write variations
            // BTreeMap has sorted keys
            for (_, var) in vars {
                match var {
                    Variation::Substitution(sub) => {
                        paint_sub(&mut worksheet, &format_of.clone(), &mut opt, &sub).unwrap()
                    }
                    Variation::Indel(indel) => {
                        paint_indel(&mut worksheet, format_of.clone(), &mut opt, &indel)?
                    } // Indel
                } // Match

                // Increase column cursor
                opt.col_cursor += 1;
                // Wrap
                if opt.col_cursor > opt.wrap {
                    opt.col_cursor = 1;
                    opt.sec_cursor += 1;
                }
            }

            opt.sec_cursor += 1;
        } // block
    }

    worksheet.set_column_width(0, opt.max_name_len as f64)?;
    for i in 1..=(opt.wrap + 3) {
        worksheet.set_column_width(i, 1.6)?;
    }

    // Save the file to disk.
    workbook.save(outfile)?;

    Ok(())
}

fn paint_name(
    worksheet: &mut Worksheet,
    format_of: &BTreeMap<String, Format>,
    opt: &mut Opt,
    block: &hnsm::FasBlock,
) -> anyhow::Result<()> {
    for i in 1..=block.entries.len() {
        let pos_row = opt.sec_height * (opt.sec_cursor - 1);

        let rg = block.entries[i - 1].range().to_string();
        worksheet.write_with_format(
            pos_row + i as u32,
            0,
            rg.clone(),
            format_of.clone().get("name").unwrap(),
        )?;

        // record max length
        opt.max_name_len = max(rg.len(), opt.max_name_len);
    }
    Ok(())
}

fn paint_indel(
    worksheet: &mut Worksheet,
    format_of: BTreeMap<String, Format>,
    opt: &mut Opt,
    indel: &hnsm::Indel,
) -> anyhow::Result<()> {
    let mut pos_row = opt.sec_height * (opt.sec_cursor - 1);

    // how many column does this indel take up
    let col_taken = indel.length.min(3) as u16;

    // if exceed the wrap limit, start a new section
    if opt.col_cursor + col_taken > opt.wrap {
        opt.col_cursor = 1;
        opt.sec_cursor += 1;
        pos_row = opt.sec_height * (opt.sec_cursor - 1);
    }

    // Write indel type and length
    let indel_string = format!("{}{}", indel.itype, indel.length);
    let format = {
        let bg_idx = if indel.occurred == "unknown" {
            "unknown".to_string()
        } else {
            let idx = u32::from_str_radix(&indel.occurred, 2)? % opt.color_loop;
            idx.to_string()
        };
        let format_key = format!("indel_{}", bg_idx);
        format_of.get(&format_key).unwrap()
    };

    for i in 1..=opt.seq_count {
        let mut flag_draw = false;
        if indel.occurred == "unknown" {
            flag_draw = true;
        } else {
            let occ = indel.occurred.chars().nth(i as usize - 1).unwrap();
            if occ == '1' {
                flag_draw = true;
            }
        }

        if !flag_draw {
            continue;
        }

        if col_taken == 1 {
            // Write position
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;

            // indel occurred lineages
            worksheet.write_with_format(pos_row + i, opt.col_cursor, &indel_string, format)?;
        } else if col_taken == 2 {
            // Write indel start position
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;
            // Write indel end position
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 1,
                indel.end,
                format_of.get("pos").unwrap(),
            )?;

            // Merge indel positions
            worksheet.merge_range(
                pos_row + i,
                opt.col_cursor,
                pos_row + i,
                opt.col_cursor + 1,
                &indel_string,
                format,
            )?;
        } else {
            // Write indel start position
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;
            // Write middle hyphen
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 1,
                "|",
                format_of.get("pos").unwrap(),
            )?;
            // Write indel end position
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 2,
                indel.end,
                format_of.get("pos").unwrap(),
            )?;

            // Merge indel positions
            worksheet.merge_range(
                pos_row + i,
                opt.col_cursor,
                pos_row + i,
                opt.col_cursor + 2,
                &indel_string,
                format,
            )?;
        }
    }
    Ok(())
}

fn paint_sub(
    worksheet: &mut Worksheet,
    format_of: &BTreeMap<String, Format>,
    opt: &mut Opt,
    sub: &hnsm::Substitution,
) -> anyhow::Result<()> {
    let pos_row = opt.sec_height * (opt.sec_cursor - 1);

    // Write position
    worksheet.write_with_format(
        pos_row,
        opt.col_cursor,
        sub.pos,
        format_of.get("pos").unwrap(),
    )?;

    for i in 1..=opt.seq_count {
        let base = sub.bases.chars().nth(i as usize - 1).unwrap();
        let occurred = if sub.pattern == "unknown" {
            '0'
        } else {
            sub.pattern.chars().nth(i as usize - 1).unwrap()
        };

        let base_color = if occurred == '1' {
            let bg_idx = u32::from_str_radix(&sub.pattern, 2)? % opt.color_loop;
            format!("sub_{}_{}", base, bg_idx)
        } else {
            format!("sub_{}_unknown", base)
        };
        let format = format_of.get(&base_color).unwrap();
        worksheet.write_with_format(pos_row + i, opt.col_cursor, base.to_string(), format)?;
    }

    // Outgroup bases with no bg colors
    if opt.is_outgroup {
        let base_color = format!("sub_{}_unknown", sub.obase);
        let format = format_of.get(&base_color).unwrap();
        worksheet.write_with_format(
            pos_row + opt.seq_count + 1,
            opt.col_cursor,
            sub.obase.clone(),
            format,
        )?;
    }
    Ok(())
}

/// Get variations (substitutions and indels)
fn get_vars(
    seqs: &[&[u8]],
    is_outgroup: bool,
    is_indel: bool,
    no_single: bool,
    no_complex: bool,
    min_freq: Option<f64>,
    max_freq: Option<f64>,
) -> anyhow::Result<BTreeMap<i32, Variation>> {
    let mut vars = BTreeMap::new();

    let mut seq_count = seqs.len();
    let out_seq = if is_outgroup {
        seq_count -= 1;
        Some(seqs[seq_count])
    } else {
        None
    };

    // Get substitutions
    let subs = if is_outgroup {
        let mut unpolarized = hnsm::get_subs(&seqs[..seq_count])?;
        hnsm::polarize_subs(&mut unpolarized, out_seq.unwrap());
        unpolarized
    } else {
        hnsm::get_subs(&seqs)?
    };

    for sub in subs {
        // Filter substitutions
        if no_single && sub.freq <= 1 {
            continue;
        }
        if no_complex && sub.freq == -1 {
            continue;
        }
        if let Some(min) = min_freq {
            if (sub.freq as f64) / (seq_count as f64) < min {
                continue;
            }
        }
        if let Some(max) = max_freq {
            if (sub.freq as f64) / (seq_count as f64) > max {
                continue;
            }
        }

        vars.insert(sub.pos, Variation::Substitution(sub));
    }

    // Get indels
    if is_indel {
        let indels = if is_outgroup {
            let mut unpolarized = hnsm::get_indels(&seqs[..seq_count])?;
            hnsm::polarize_indels(&mut unpolarized, out_seq.unwrap());
            unpolarized
        } else {
            hnsm::get_indels(&seqs)?
        };

        for indel in indels {
            // Filter indels
            if no_single && indel.freq <= 1 {
                continue;
            }
            if no_complex && indel.freq == -1 {
                continue;
            }
            if let Some(min) = min_freq {
                if (indel.freq as f64) / (seq_count as f64) < min {
                    continue;
                }
            }
            if let Some(max) = max_freq {
                if (indel.freq as f64) / (seq_count as f64) > max {
                    continue;
                }
            }

            vars.insert(indel.start, Variation::Indel(indel));
        }
    }

    Ok(vars)
}

fn create_formats() -> BTreeMap<String, Format> {
    let mut format_of: BTreeMap<String, Format> = BTreeMap::new();

    // species names
    format_of.insert(
        "name".to_string(),
        Format::new().set_font_name("Courier New").set_font_size(10),
    );

    // align positions of variations
    format_of.insert(
        "pos".to_string(),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(8)
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center)
            .set_rotation(90),
    );

    // the standard Excel colors in the range 8..63

    // 15 colors
    let bg_colors: Vec<u32> = vec![
        0xC0C0C0, // Gray-25%, silver, 22
        0xFFFF99, // Light Yellow, 43
        0xCCFFCC, // Light Green, 42
        0xCCFFFF, // Lite Turquoise, 27
        0x99CCFF, // Pale Blue, 44
        0xCC99FF, // Lavender, 46
        0xFFCC99, // Tan, 47
        0x9999FF, // Periwinkle, 24
        0x33CCCC, // Aqua, 49
        0xFFCC00, // Gold, 51
        0xFF99CC, // Rose, 45
        0xFF9900, // Light Orange, 52
        0xFFFFCC, // Ivory, 26
        0xFF8080, // Coral, 29
        0xCCCCFF, // Ice Blue, 31
                  // 0x0066CC,       // Ocean Blue, 30
                  // 0xCCFFFF,       // Light Turquoise, again, 41
                  // 0x3366FF,       // Light Blue, 48
                  // 0x99CC00,       // Lime, 50
                  // 0x666699,       // Blue-Gray, 54
                  // 0x333399,       // Indigo, 62
    ];

    // font colors
    let sub_fc_of: BTreeMap<String, u32> = BTreeMap::from([
        ("A".to_string(), 0x003300), // Dark Green, 58
        ("C".to_string(), 0x000080), // Dark Blue, Navy, 18
        ("G".to_string(), 0x660066), // Dark Purple, 28
        ("T".to_string(), 0x800000), // Dark Red, Brown, 16
        ("N".to_string(), 0x000000), // Black, 8
        ("-".to_string(), 0x000000), // Black, 8
    ]);

    // sub _ base       _ color rotation
    // sub - font color - background color
    for fc in sub_fc_of.keys() {
        format_of.insert(
            format!("sub_{}_{}", fc, "unknown"),
            Format::new()
                .set_font_name("Courier New")
                .set_font_size(10)
                .set_align(FormatAlign::VerticalCenter)
                .set_align(FormatAlign::Center)
                .set_font_color(*sub_fc_of.get(fc).unwrap())
                .set_background_color(Color::White),
        );

        for i in 0..bg_colors.len() {
            let key = format!("sub_{}_{}", fc, i);
            format_of.insert(
                key,
                Format::new()
                    .set_font_name("Courier New")
                    .set_font_size(10)
                    .set_align(FormatAlign::VerticalCenter)
                    .set_align(FormatAlign::Center)
                    .set_font_color(*sub_fc_of.get(fc).unwrap())
                    .set_background_color(*bg_colors.get(i).unwrap()),
            );
        }
    }

    format_of.insert(
        "sub_-".to_string(),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(10)
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center),
    );

    for i in 0..bg_colors.len() {
        let key = format!("indel_{}", i);
        format_of.insert(
            key,
            Format::new()
                .set_font_name("Courier New")
                .set_font_size(10)
                .set_bold()
                .set_align(FormatAlign::VerticalCenter)
                .set_align(FormatAlign::Center)
                .set_background_color(*bg_colors.get(i).unwrap()),
        );
    }
    format_of.insert(
        format!("indel_{}", "unknown"),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(10)
            .set_bold()
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center)
            .set_background_color(Color::White),
    );

    format_of
}
