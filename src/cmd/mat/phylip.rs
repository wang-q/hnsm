use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("phylip")
        .about("Convert pairwise distances to a phylip distance matrix")
        .after_help(
            r###"
Conversion modes:
    * full:  a full distance matrix
    * lower: a lower-triangular matrix
    * strict: a strict phylip distance matrix

Input format:
    * Tab-separated values (TSV)
    * Three columns: name1, name2, distance

Output format:
    * PHYLIP distance matrix format
    * First line contains the number of sequences
    * Each subsequent line contains a sequence name and distances
    * In strict mode:
        - Names are limited to 10 characters
        - Names are left-aligned and padded with spaces
        - Distances are space-separated with 6 decimal places

Examples:
    1. Create a full matrix:
       hnsm mat phylip input.tsv -o output.phy

    2. Create a lower-triangular matrix:
       hnsm mat phylip input.tsv --mode lower -o output.phy

    3. Create a strict PHYLIP matrix:
       hnsm mat phylip input.tsv --mode strict -o output.phy
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise distances"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("full"),
                    builder::PossibleValue::new("lower"),
                    builder::PossibleValue::new("strict"),
                ])
                .default_value("full")
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
    // Load matrix from pairwise distances
    let matrix = hnsm::NamedMatrix::from_pair_scores(infile, opt_same, opt_missing);
    let names = matrix.get_names();
    let size = matrix.size();

    // Write sequence count
    writer.write_fmt(format_args!("{:>4}\n", size))?;

    for i in 0..size {
        match opt_mode.as_str() {
            "full" => {
                writer.write_fmt(format_args!("{}", names[i]))?;
                for j in 0..size {
                    writer.write_fmt(format_args!("\t{}", matrix.get(i, j)))?;
                }
            }
            "lower" => {
                writer.write_fmt(format_args!("{}", names[i]))?;
                for j in 0..i {
                    writer.write_fmt(format_args!("\t{}", matrix.get(i, j)))?;
                }
            }
            "strict" => {
                // Strict mode: names limited to 10 chars, left-aligned with space padding
                writer.write_fmt(format_args!(
                    "{:<10}",
                    names[i].chars().take(10).collect::<String>()
                ))?;
                for j in 0..size {
                    writer.write_fmt(format_args!(" {:.6}", matrix.get(i, j)))?;
                }
            }
            _ => unreachable!(),
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}
