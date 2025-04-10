use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("size")
        .about("Count total bases in FA file(s)")
        .after_help(
            r###"
This command counts the total number of bases in one or more FA files. It outputs the sequence name
and its length in a tab-separated format.

Examples:
1. Count bases in a single FASTA file:
   hnsm size input.fa

2. Count bases in multiple FASTA files:
   hnsm size input1.fa input2.fa

3. Save the output to a file:
   hnsm size input.fa -o output.tsv

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FA file(s) to process"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;

            writer.write_fmt(format_args!("{}\t{}\n", name, record.sequence().len()))?;
        }
    }

    Ok(())
}
