use clap::*;
use noodles_fasta as fasta;
use std::collections::HashSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse complement a FA file")
        .after_help(
            r###"
* [list] is optional, only RC sequences listed in this file

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
                .required(false)
                .index(2)
                .help("One name per line"),
        )
        .arg(
            Arg::new("consistent")
                .long("consistent")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Keep the name consistent (don't prepend RC_)"),
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

    let is_consistent = args.get_flag("consistent");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    let mut set_list: HashSet<String> = HashSet::new();
    if args.contains_id("list.txt") {
        set_list = intspan::read_first_column(args.get_one::<String>("list.txt").unwrap())
            .into_iter()
            .collect();
    }

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        let mut name = String::from_utf8(record.name().into()).unwrap();

        if args.contains_id("list.txt") && !set_list.contains(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        if !is_consistent {
            name = format!("RC_{}", name);
        }

        let definition = fasta::record::Definition::new(&*name, None);

        let seq_rc: fasta::record::Sequence = record
            .sequence()
            .complement()
            .rev()
            .collect::<Result<_, _>>()?;
        let record_rc = fasta::Record::new(definition, seq_rc);
        fa_out.write_record(&record_rc)?;
    }

    Ok(())
}
