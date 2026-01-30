use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dbscan")
        .about("DBSCAN clustering based on pairwise distances")
        .after_help(
            r###"
Density-based spatial clustering of applications with noise (DBSCAN).

Note: The input file should contain pairwise distances (lower is better), NOT similarities.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Note:
For the 'pair' format, the representative point is the medoid (point with minimum sum of distances to other cluster members).
If there are ties, the alphabetically first member is chosen.

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
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    let opt_format = args.get_one::<String>("format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let opt_eps = *args.get_one::<f32>("eps").unwrap();
    let opt_min_points = *args.get_one::<usize>("min_points").unwrap();

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // 2. Load Matrix
    //----------------------------

    // Load matrix from pairwise distances
    let (matrix, names) = intspan::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing);

    //----------------------------
    // 3. Clustering
    //----------------------------
    let mut dbscan = hnsm::Dbscan::new(opt_eps, opt_min_points);
    dbscan.perform_clustering(&matrix);
    let mut clusters = dbscan.results_cluster();

    // Sort members within each cluster
    for c in &mut clusters {
        c.sort_by_key(|&idx| &names[idx]);
    }

    // Sort clusters: first by first member name, then by size (descending)
    clusters.sort_by(|a, b| {
        let name_a = &names[a[0]];
        let name_b = &names[b[0]];
        match b.len().cmp(&a.len()) {
            std::cmp::Ordering::Equal => name_a.cmp(name_b),
            other => other,
        }
    });

    //----------------------------
    // 4. Output
    //----------------------------
    match opt_format.as_str() {
        "cluster" => {
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
            for component in clusters {
                // Find medoid
                let mut best_rep = *component.first().unwrap();
                let mut min_sum = f32::MAX;

                for &candidate in &component {
                    let mut current_sum = 0.0;
                    for &member in &component {
                        current_sum += matrix.get(candidate, member);
                    }
                    if current_sum < min_sum {
                        min_sum = current_sum;
                        best_rep = candidate;
                    }
                    // Implicit tie-break (first one wins, which is alphabetical)
                }

                let rep_name = &names[best_rep];
                for &member_idx in &component {
                    writer.write_fmt(format_args!("{}\t{}\n", rep_name, names[member_idx]))?;
                }
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
