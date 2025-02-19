use clap::*;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replace headers of a FA file based on a TSV mapping")
        .after_help(
            r###"
This command replaces sequence headers in a FA file based on a TSV mapping file.
The TSV file should contain two or more columns: the original name and the replacement name.
If more than two columns are provided, the sequence will be duplicated for each replacement name.
Multiple lines of the same original_name will also duplicate the record.

The TSV file format:
    seq1    replace_name    more_replace_name
    seq2    replace_name
    seq2    another_replace_name

Examples:
1. Replace headers using a TSV file:
   hnsm replace input.fa replace.tsv -o output.fa

2. Only output sequences listed in the TSV file (like `hnsm some`):
   hnsm replace input.fa replace.tsv -s -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FA file to process"),
        )
        .arg(
            Arg::new("replace.tsv")
                .required(true)
                .index(2)
                .help("TSV file containing original and replacement names"),
        )
        .arg(
            Arg::new("some")
                .long("some")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Only output sequences listed in the TSV file, like `hnsm some`"),
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

    let replace_of = read_replaces(args.get_one::<String>("replace.tsv").unwrap());
    let is_some = args.get_flag("some");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;

        if let Some(new_names) = replace_of.get(&name) {
            for el in new_names {
                let definition = noodles_fasta::record::Definition::new(&**el, None);
                let record_replace =
                    noodles_fasta::Record::new(definition, record.sequence().clone());
                fa_out.write_record(&record_replace)?;
            }
        } else if !is_some {
            fa_out.write_record(&record)?;
        }
    }

    Ok(())
}

// Read the replacement mappings from a TSV file
fn read_replaces(input: &str) -> HashMap<String, Vec<String>> {
    let mut replaces: HashMap<String, Vec<String>> = HashMap::new();

    for line in intspan::read_lines(input) {
        let mut fields: Vec<&str> = line.split('\t').collect();

        let key = fields[0].to_string();
        let mut others = fields
            .split_off(1)
            .iter()
            .map(|s| (*s).to_string())
            .collect::<Vec<_>>();

        // Use the entry method to check if the key exists
        replaces.entry(key.clone()).or_default().append(&mut others);
    }

    replaces
}
