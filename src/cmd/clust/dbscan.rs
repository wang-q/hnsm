use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dbscan")
        .about("DBSCAN clustering based on pairwise distances")
        .after_help(
            r###"
Density-based spatial clustering of applications with noise (DBSCAN).

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise distances in .tsv format"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("cluster"),
                    builder::PossibleValue::new("pair"),
                ])
                .default_value("cluster")
                .help("Output format for clustering results"),
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
            Arg::new("eps")
                .long("eps")
                .num_args(1)
                .default_value("0.05")
                .value_parser(value_parser!(f32))
                .help("The maximum distance between two points for DBSCAN clustering"),
        )
        .arg(
            Arg::new("min_points")
                .long("min_points")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Minimum number of points to form a dense region in DBSCAN"),
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

    let opt_format = args.get_one::<String>("format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let opt_eps = *args.get_one::<f32>("eps").unwrap();
    let opt_min_points = *args.get_one::<usize>("min_points").unwrap();

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------

    // Load matrix from pairwise distances
    let (matrix, names) = hnsm::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing);

    let mut dbscan = hnsm::Dbscan::new(opt_eps, opt_min_points);
    let _ = dbscan.perform_clustering(&matrix);
    match opt_format.as_str() {
        "cluster" => {
            let clusters = dbscan.results_cluster();
            for c in clusters {
                writer.write_fmt(format_args!(
                    "{}\n",
                    c.iter()
                        .map(|&num| names[num].clone())
                        .collect::<Vec<_>>()
                        .join("\t")
                ))?;
            }
        }
        "pair" => {
            let rep_points = dbscan.results_pair(&matrix);
            for (rep, point) in rep_points {
                writer.write_fmt(format_args!("{}\t{}\n", names[rep], names[point]))?;
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
