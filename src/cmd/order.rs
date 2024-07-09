use clap::*;
use noodles_fasta as fasta;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("order")
        .about("Extract some FA records by the given order")
        .after_help(
            r###"
* Loads all sequences in memory, thus consuming more memory

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("list.txt")
                .required(true)
                .index(2)
                .help("One name per line"),
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
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = fasta::io::Reader::new(reader);

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_with_writer(writer);

    let vec_list = intspan::read_first_column(args.get_one::<String>("list.txt").unwrap());
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
