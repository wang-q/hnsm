use clap::*;
use std::collections::HashSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rc")
        .about("Reverse complement sequences in FA file(s)")
        .after_help(
            r###"
This command reverse complements DNA sequences in FA files.

Features:
* Process all sequences or only selected ones
* Optionally prefix names with 'RC_'
* Handles IUPAC ambiguous codes correctly
* Preserves case (upper/lower) of bases

Notes:
* Case-sensitive name matching when using list
* Empty lines and lines starting with '#' are ignored in list
* Supports both plain text and gzipped (.gz) files
* Non-IUPAC characters are preserved as-is

Examples:
1. Reverse complement all sequences:
   hnsm rc input.fa -o output.fa

2. Only process listed sequences:
   hnsm rc input.fa list.txt -o output.fa

3. Keep original names (no 'RC_' prefix):
   hnsm rc input.fa -c -o output.fa

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
                .required(false)
                .index(2)
                .help("File containing one sequence name per line (optional)"),
        )
        .arg(
            Arg::new("consistent")
                .long("consistent")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Keep the name consistent (don't prepend 'RC_')"),
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

    let is_consistent = args.get_flag("consistent");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    let set_list: HashSet<String> = if args.contains_id("list.txt") {
        intspan::read_first_column(args.get_one::<String>("list.txt").unwrap())
            .into_iter()
            .collect()
    } else {
        HashSet::new()
    };

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if args.contains_id("list.txt") && !set_list.contains(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        let new_name = if is_consistent {
            name
        } else {
            format!("RC_{}", name)
        };

        let definition = noodles_fasta::record::Definition::new(&*new_name, None);
        let seq_rc: noodles_fasta::record::Sequence = record
            .sequence()
            .complement()
            .rev()
            .collect::<Result<_, _>>()?;
        let record_rc = noodles_fasta::Record::new(definition, seq_rc);
        fa_out.write_record(&record_rc)?;
    }

    Ok(())
}
