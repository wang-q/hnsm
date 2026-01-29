use clap::*;

pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merge fragmented synteny blocks")
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input synteny file"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .default_value("merged.tsv")
                .help("Output filename"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let outfile = matches.get_one::<String>("outfile").unwrap();
    
    println!("Merging blocks from {} to {}", infile, outfile);
    Ok(())
}
