use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("manifold")
        .about("Manifold learning based on pairwise distances")
        .after_help(
            r###"
modes:
    * pcoa, a.k.a. classical multidimensional scaling (cMDS)
    * tsne
    * umap

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
                .value_parser([builder::PossibleValue::new("pcoa")])
                .default_value("pcoa")
                .help("Reduction method"),
        )
        .arg(
            Arg::new("dim")
                .long("dim")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("The number of dimensions"),
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

    let opt_dim = *args.get_one::<usize>("dim").unwrap();

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    // Read pair scores from a TSV file
    let pair_scores = hnsm::load_file(infile);
    let (matrix, index_name) = hnsm::populate(&pair_scores);
    let size = matrix.size();

    match opt_mode.as_str() {
        "pcoa" => {
            let mut dmatrix = pcoa::nalgebra::DMatrix::from_element(size, size, 1.);

            for row in 0..size {
                for col in 0..size {
                    dmatrix[(row, col)] = matrix.get(row, col);
                }
            }

            let coords_matrix = pcoa::apply_pcoa(dmatrix, opt_dim).expect("cannot apply PCoA");

            let coords_matrix = coords_matrix.transpose();
            let xs: Vec<_> = coords_matrix.column(0).iter().copied().collect();
            let ys: Vec<_> = coords_matrix.column(1).iter().copied().collect();

            for ((x, y), n) in std::iter::zip(xs, ys).zip(index_name) {
                writer.write_fmt(format_args!("{}\t{}\t{}\n", n, x, y))?;
            }
        }
        _ => unreachable!(),
    }

    Ok(())
}
