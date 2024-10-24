use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cluster")
        .about("Clustering based on pairwise distances")
        .after_help(
            r###"
modes:
    * matrix: write out a distance matrix
    * dbscan

format:
    * cluster: a line contains points of one cluster
    * pair: lines of multiple (representative point, cluster member) pairs

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("matrix"),
                    builder::PossibleValue::new("dbscan"),
                ])
                .default_value("matrix")
                .help("Clustering method"),
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
                .help("Output formats"),
        )
        .arg(
            Arg::new("eps")
                .long("eps")
                .num_args(1)
                .default_value("0.05")
                .value_parser(value_parser!(f32))
                .help("The maximum distance between two points"),
        )
        .arg(
            Arg::new("min_points")
                .long("min_points")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("core point"),
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
    let opt_format = args.get_one::<String>("format").unwrap();

    let opt_eps = *args.get_one::<f32>("eps").unwrap();
    let opt_min_points = *args.get_one::<usize>("min_points").unwrap();

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    // Read pair scores from a TSV file
    let pair_scores = hnsm::load_file(infile);
    let (matrix, index_name) = hnsm::populate(&pair_scores);
    let size = matrix.size();

    match opt_mode.as_str() {
        "matrix" => {
            // Print the scoring matrix
            for i in 0..size {
                writer.write_fmt(format_args!("{}\t", index_name.get(i).unwrap()))?;
                for j in 0..size {
                    writer.write_fmt(format_args!("{}\t", matrix.get(i, j)))?;
                }
                writer.write_fmt(format_args!("\n"))?;
            }
        }
        "dbscan" => {
            let mut dbscan = hnsm::Dbscan::new(opt_eps, opt_min_points);
            let _ = dbscan.perform_clustering(&matrix);
            match opt_format.as_str() {
                "cluster" => {
                    let clusters = dbscan.results_cluster();
                    for c in clusters {
                        writer.write_fmt(format_args!(
                            "{}\n",
                            c.iter()
                                .map(|&num| index_name.get(num).unwrap().to_string())
                                .collect::<Vec<_>>()
                                .join("\t")
                        ))?;
                    }
                }
                "pair" => {
                    let rep_points = dbscan.results_pair(&matrix);
                    for (rep, point) in rep_points {
                        writer.write_fmt(format_args!(
                            "{}\t{}\n",
                            index_name.get(rep).unwrap(),
                            index_name.get(point).unwrap()
                        ))?;
                    }
                }
                _ => unreachable!(),
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
