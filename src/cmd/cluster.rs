use clap::*;
use std::collections::HashMap;
use std::io::BufRead;

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
    let pair_scores = load_file(infile);
    let (matrix, index_name) = populate(&pair_scores);
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
            let mut dbscan = hnsm::DBSCAN::new(opt_eps, opt_min_points);
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

fn populate(
    pair_scores: &Vec<((String, String), f32)>,
) -> (hnsm::SymmetricMatrix<f32>, Vec<String>) {
    // Create a mapping from string identifiers to indices
    let mut index_map = HashMap::new();
    let mut index_name = vec![];
    let mut current_index = 0usize;

    for ((n1, n2), _) in pair_scores {
        if !index_map.contains_key(n1) {
            index_map.insert(n1.clone(), current_index);
            current_index += 1;
            index_name.push(n1.clone());
        }
        if !index_map.contains_key(n2) {
            index_map.insert(n2.clone(), current_index);
            current_index += 1;
            index_name.push(n2.clone());
        }
    }

    // Determine the size of the matrix
    let size = index_map.len();

    // Create a new scoring matrix
    let mut matrix: hnsm::SymmetricMatrix<f32> = hnsm::SymmetricMatrix::new(size);

    // Populate the scoring matrix with the pair scores
    for ((n1, n2), score) in pair_scores {
        let i = index_map[n1];
        let j = index_map[n2];
        matrix.set(i, j, *score);
    }

    (matrix, index_name)
}

fn load_file(infile: &str) -> Vec<((String, String), f32)> {
    let mut pair_scores = Vec::new();
    let reader = intspan::reader(infile);
    for line in reader.lines().map_while(Result::ok) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let n1 = fields[0].to_string();
            let n2 = fields[1].to_string();
            let score: f32 = fields[2].parse::<f32>().unwrap();
            pair_scores.push(((n1, n2), score));
        }
    }
    pair_scores
}
