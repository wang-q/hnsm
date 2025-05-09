use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("one")
        .about("Extract one FA record by name")
        .after_help(
            r###"
This command extracts a single FA record from an input file based on the provided sequence name.

Notes:
* Case-sensitive name matching
* Stops after finding the first match
* Supports both plain text and gzipped (.gz) files

Examples:
1. Extract a record by name:
   hnsm one input.fa seq1

2. Save to a file:
   hnsm one input.fa seq1 -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FA file to process"),
        )
        .arg(
            Arg::new("name")
                .required(true)
                .index(2)
                .help("Name of the sequence to extract"),
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
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = noodles_fasta::io::Reader::new(reader);

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    let name = args.get_one::<String>("name").unwrap();

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let this_name = String::from_utf8(record.name().into())?;

        if this_name == *name {
            fa_out.write_record(&record)?;
            break;
        }
    }

    Ok(())
}
