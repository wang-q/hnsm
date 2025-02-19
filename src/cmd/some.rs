use clap::*;
use std::collections::HashSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("some")
        .about("Extract some FA records based on a list of names")
        .after_help(
            r###"
This command extracts FASTA records from an input file based on a list of sequence names.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* Supports both plain text and gzipped (.gz) files

Examples:
1. Extract sequences listed in list.txt:
   hnsm some input.fa list.txt

2. Extract sequences NOT in list.txt:
   hnsm some input.fa list.txt -i

3. Process gzipped files:
   hnsm some input.fa.gz list.txt -o output.fa.gz

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FA file to process"),
        )
        .arg(
            Arg::new("list.txt")
                .required(true)
                .index(2)
                .help("File containing one sequence name per line"),
        )
        .arg(
            Arg::new("invert")
                .long("invert")
                .short('i')
                .action(ArgAction::SetTrue)
                .help("Invert selection: output sequences NOT in the list"),
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
    let is_invert = args.get_flag("invert");

    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = noodles_fasta::io::Reader::new(reader);

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    //----------------------------
    // Load list
    //----------------------------
    let set_list: HashSet<String> =
        intspan::read_first_column(args.get_one::<String>("list.txt").unwrap())
            .into_iter()
            .collect();

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if set_list.contains(&name) != is_invert {
            fa_out.write_record(&record)?;
        }
    }

    Ok(())
}
