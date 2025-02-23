use clap::*;
use std::io::{BufRead, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("convert")
        .about("Conversion between pairwise distances and a distance matrix")
        .after_help(
            r###"
Conversion modes:
    * matrix: Convert pairwise distances to a full distance matrix.
    * lower: Convert pairwise distances to a lower-triangular matrix.
    * pair: Convert a (lower-triangular) relaxed PHYLIP distance matrix to pairwise distances.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise distances or a PHYLIP distance matrix"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("matrix"),
                    builder::PossibleValue::new("lower"),
                    builder::PossibleValue::new("pair"),
                ])
                .default_value("matrix")
                .help("Conversion mode"),
        )
        .arg(
            Arg::new("same")
                .long("same")
                .num_args(1)
                .default_value("0.0")
                .value_parser(value_parser!(f32))
                .help("Default score of identical element pairs"),
        )
        .arg(
            Arg::new("missing")
                .long("missing")
                .num_args(1)
                .default_value("1.0")
                .value_parser(value_parser!(f32))
                .help("Default score of missing pairs"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_mode = args.get_one::<String>("mode").unwrap();

    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    if opt_mode.as_str() == "pair" {
        // Process PHYLIP matrix to pairwise distances
        let reader = intspan::reader(infile);
        let mut lines = reader.lines();
        let mut names = Vec::new();

        // Attempt to read the first line as the number of sequences
        if let Some(Ok(line)) = lines.next() {
            if line.trim().parse::<usize>().is_err() {
                // If the first line is not a number, treat it as a data line
                process_phylip_line(&line, &mut names, &mut writer)?;
            }
        }

        for line in lines.map_while(Result::ok) {
            process_phylip_line(&line, &mut names, &mut writer)?;
        }

        return Ok(());
    }

    // Convert pairwise distances to a matrix
    let (pair_scores, index_name) = hnsm::load_pair_scores(infile);
    let matrix = hnsm::populate_matrix(&pair_scores, &index_name, opt_same, opt_missing);
    let size = matrix.size();

    for i in 0..size {
        writer.write_fmt(format_args!("{}", index_name.get(i).unwrap()))?;
        match opt_mode.as_str() {
            "matrix" => {
                for j in 0..size {
                    writer.write_fmt(format_args!("\t{}", matrix.get(i, j)))?;
                }
            }
            "lower" => {
                for j in 0..i {
                    writer.write_fmt(format_args!("\t{}", matrix.get(i, j)))?;
                }
            }
            _ => unreachable!(),
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}

// Process a single line of the PHYLIP matrix and output pairwise distances
fn process_phylip_line(
    line: &str,
    names: &mut Vec<String>,
    writer: &mut Box<dyn Write>,
) -> anyhow::Result<()> {
    let parts: Vec<&str> = line.trim().split_whitespace().collect();
    if !parts.is_empty() {
        let name = parts[0].to_string();
        names.push(name.clone());

        // Read lower-triangle distances
        let distances: Vec<f32> = parts[1..=names.len()]
            .iter()
            .map(|&s| s.parse().unwrap())
            .collect();

        for (i, &distance) in distances.iter().enumerate() {
            writer.write_fmt(format_args!("{}\t{}\t{}\n", names[i], name, distance))?;
        }
    }

    Ok(())
}
