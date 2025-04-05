use clap::*;
use std::io::BufRead;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cc")
        .about("Connected components clustering")
        .after_help(
            r###"
Ignores scores and writes all connected components.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input file containing pairwise distances in .tsv format"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    let mut names = indexmap::IndexSet::new();

    let mut graph = petgraph::graphmap::UnGraphMap::<_, ()>::new();

    let reader = intspan::reader(infile);
    for line in reader.lines().map_while(Result::ok) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            names.insert(fields[0].to_string());
            names.insert(fields[1].to_string());
        }

        graph.add_edge(
            names.get_index_of(fields[0]).unwrap(),
            names.get_index_of(fields[1]).unwrap(),
            (),
        );
    }

    let mut scc = petgraph::algo::tarjan_scc(&graph);

    // First sort members within each component alphabetically
    for cc in &mut scc {
        cc.sort_by_key(|&idx| names.get_index(idx).unwrap().as_str());
    }

    // Then sort components by first member alphabetically
    scc.sort_by_key(|cc| names.get_index(cc[0]).unwrap().as_str());

    // Finally sort by size (descending) while maintaining alphabetical order for same size
    scc.sort_by_key(|cc| std::cmp::Reverse(cc.len()));

    // Output each component
    for cc in &scc {
        writer.write_fmt(format_args!(
            "{}\n",
            cc.iter()
                .map(|&idx| names.get_index(idx).unwrap().as_str())
                .collect::<Vec<_>>()
                .join("\t")
        ))?;
    }

    Ok(())
}
