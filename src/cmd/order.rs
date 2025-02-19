use clap::*;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Extract some FA records in the order specified by a list")
        .after_help(
            r###"
This command extracts FA records from an input file in the order specified by a list of sequence names.

Notes:
* Case-sensitive name matching
* One sequence name per line in the list file
* Empty lines and lines starting with '#' are ignored
* All sequences are loaded into memory
* Supports both plain text and gzipped (.gz) files
* Missing sequences in the input file are silently skipped

Examples:
1. Extract sequences in order specified by list.txt:
   hnsm order input.fa list.txt

2. Process gzipped files:
   hnsm order input.fa.gz list.txt -o output.fa.gz

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

    let list: indexmap::IndexSet<_> =
        intspan::read_first_column(args.get_one::<String>("list.txt").unwrap())
            .into_iter()
            .collect();

    //----------------------------
    // Process
    //----------------------------
    // Load records into a BTreeMap for efficient lookup
    let mut record_of = BTreeMap::new();

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if list.contains(&name) {
            record_of.insert(name, record);
        }
    }

    for name in list.iter() {
        if let Some(record) = record_of.get(name) {
            fa_out.write_record(record)?;
        }
    }

    Ok(())
}
