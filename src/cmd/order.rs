use clap::*;
use noodles_fasta as fasta;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Extract some FA records in the order specified by a list")
        .after_help(
            r###"
This command extracts FA records from an input file in the order specified by a list of sequence names.
All sequences are loaded into memory, so this command may consume significant memory for large files.

Examples:
1. Extract sequences in the order specified by list.txt:
   hnsm order input.fa list.txt -o output.fa

2. Output to stdout:
   hnsm order input.fa list.txt

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
    let mut fa_in = fasta::io::Reader::new(reader);

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    let vec_list = intspan::read_first_column(args.get_one::<String>("list.txt").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    // Load records into a BTreeMap for efficient lookup
    let mut record_of = BTreeMap::new();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into()).unwrap();
        if vec_list.contains(&name) {
            record_of.insert(name, record);
        }
    }

    for el in vec_list.iter() {
        if record_of.contains_key(el) {
            let record = record_of.get(el).unwrap();
            fa_out.write_record(record)?;
        }
    }

    Ok(())
}
