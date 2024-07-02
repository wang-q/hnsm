use clap::*;
use std::io::BufRead;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sixframe")
        .about("Six-Frame Translation")
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("header")
                .long("header")
                .short('H')
                .action(ArgAction::SetTrue)
                .help("Treat the first line of each file as a header"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Set the minimal length of AA sequence"),
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
    // Options
    //----------------------------
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        for line in reader.lines().map_while(Result::ok) {}
    }

    Ok(())
}
