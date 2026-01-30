use clap::*;
use std::io::Write;
use intspan::ScoringMatrix;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mcl")
        .about("Markov Clustering Algorithm (MCL)")
        .after_help(
            r###"
MCL is a fast and scalable unsupervised cluster algorithm for graphs (also known as networks) based on simulation of (stochastic) flow in graphs.

It is particularly useful for clustering protein interaction networks or similarity networks.

Note: The input file should contain similarity scores (higher is better), NOT distances.

Output formats:
    * cluster: Each line contains points of one cluster.
    * pair: Each line contains a (representative point, cluster member) pair.

Reference:
Stijn van Dongen, Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, May 2000.
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise similarities (edge weights) in .tsv format"),
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
                .default_value("1.0")
                .value_parser(value_parser!(f32))
                .help("Default score of identical element pairs"),
        )
        .arg(
            Arg::new("missing")
                .long("missing")
                .num_args(1)
                .default_value("0.0")
                .value_parser(value_parser!(f32))
                .help("Default score of missing pairs"),
        )
        .arg(
            Arg::new("inflation")
                .long("inflation")
                .short('I')
                .num_args(1)
                .default_value("2.0")
                .value_parser(value_parser!(f64))
                .help("Inflation parameter. Controls the granularity of clusters. Higher values = tighter/more clusters."),
        )
        .arg(
            Arg::new("prune")
                .long("prune")
                .num_args(1)
                .default_value("1e-5")
                .value_parser(value_parser!(f64))
                .help("Pruning threshold. Matrix entries smaller than this will be set to zero."),
        )
        .arg(
            Arg::new("max_iter")
                .long("max_iter")
                .num_args(1)
                .default_value("100")
                .value_parser(value_parser!(usize))
                .help("Maximum number of iterations."),
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
    let infile = args.get_one::<String>("infile").unwrap();
    let opt_format = args.get_one::<String>("format").unwrap();
    let opt_same = *args.get_one::<f32>("same").unwrap();
    let opt_missing = *args.get_one::<f32>("missing").unwrap();
    let inflation = *args.get_one::<f64>("inflation").unwrap();
    let prune = *args.get_one::<f64>("prune").unwrap();
    let max_iter = *args.get_one::<usize>("max_iter").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();

    let mut writer = intspan::writer(outfile);

    // 1. Load Matrix
    // ScoringMatrix::from_pair_scores is only implemented for f32
    let (sm, names) = ScoringMatrix::<f32>::from_pair_scores(infile, opt_same, opt_missing);
    
    // 2. MCL Algorithm
    let mut mcl = hnsm::Mcl::new(inflation);
    mcl.set_prune_limit(prune);
    mcl.set_max_iter(max_iter);
    let clusters = mcl.perform_clustering(&sm);

    // 3. Output
    match opt_format.as_str() {
        "cluster" => {
            for component in clusters {
                let mut members: Vec<&str> = component
                    .iter()
                    .map(|&idx| names[idx].as_str())
                    .collect();
                members.sort(); // Sort for consistent output
        
                writer.write_fmt(format_args!(
                    "{}\n",
                    members.join("\t")
                ))?;
            }
        }
        "pair" => {
            for component in clusters {
                let mut members: Vec<&str> = component
                    .iter()
                    .map(|&idx| names[idx].as_str())
                    .collect();
                members.sort(); // Sort to pick a consistent representative
                
                if let Some(rep) = members.first().copied() {
                    for member in members {
                        writer.write_fmt(format_args!("{}\t{}\n", rep, member))?;
                    }
                }
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
