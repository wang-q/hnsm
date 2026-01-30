use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("k-medoids")
        .about("K-Medoids clustering")
        .visible_alias("km")
        .after_help(
            r###"
K-Medoids clustering algorithm (PAM/Lloyd-like).

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
            Arg::new("k")
                .long("k")
                .short('k')
                .num_args(1)
                .required(true)
                .value_parser(value_parser!(usize))
                .help("Number of clusters"),
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
            Arg::new("runs")
                .long("runs")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Number of random initializations"),
        )
        .arg(
            Arg::new("max_iter")
                .long("max_iter")
                .num_args(1)
                .default_value("100")
                .value_parser(value_parser!(usize))
                .help("Maximum number of iterations"),
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

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // 1. Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_k = *args.get_one::<usize>("k").unwrap();
    let opt_format = args.get_one::<String>("format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let runs = *args.get_one::<usize>("runs").unwrap();
    let max_iter = *args.get_one::<usize>("max_iter").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();

    let mut writer = intspan::writer(outfile);

    //----------------------------
    // 2. Load Matrix
    //----------------------------
    let (sm, names) = intspan::ScoringMatrix::from_pair_scores(infile, opt_same, opt_missing);

    //----------------------------
    // 3. Clustering
    //----------------------------
    let kmedoids = hnsm::KMedoids::new(opt_k, max_iter, runs);
    let mut clusters = kmedoids.perform_clustering(&sm);

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
            for component in clusters {
                let members: Vec<&str> = component
                    .iter()
                    .map(|&idx| names[idx].as_str())
                    .collect();
                // Already sorted

                writer.write_fmt(format_args!(
                    "{}\n",
                    members.join("\t")
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
                        current_sum += sm.get(candidate, member);
                    }
                    if current_sum < min_sum {
                        min_sum = current_sum;
                        best_rep = candidate;
                    }
                    // Since component is sorted, the first one with min_sum is alphabetically first.
                    // No need for explicit tie-break if we use strict <.
                    // But floating point equality...
                    else if (current_sum - min_sum).abs() < 1e-5 {
                        // Tie-break: Since we iterate in sorted order, current candidate is > best_rep.
                        // We want the smallest name. best_rep is already smaller.
                        // So do nothing.
                    }
                }
                
                let rep_name = &names[best_rep];
                let members: Vec<&str> = component
                    .iter()
                    .map(|&idx| names[idx].as_str())
                    .collect();
                // Already sorted

                for member in members {
                    writer.write_fmt(format_args!("{}\t{}\n", rep_name, member))?;
                }
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
