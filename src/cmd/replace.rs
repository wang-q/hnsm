use clap::*;
use noodles_fasta as fasta;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replace headers of a FA file")
        .after_help(
            r###"
* <replace.tsv> is a tab-separated file containing two or more fields

        original_name   replace_name    more_replace_name
        original_name   replace_name
        original_name   another_replace_name

* Three or more fields duplicates the record
* Multiple lines of the same original_name will also duplicate the record

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("replace.tsv")
                .required(true)
                .index(2)
                .help("Path to replace.tsv"),
        )
        .arg(
            Arg::new("some")
                .long("some")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Only output sequences in the list, like `hnsm some`"),
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

    let replace_of = read_replaces(args.get_one::<String>("replace.tsv").unwrap());
    let is_some = args.get_flag("some");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into()).unwrap();
        if replace_of.contains_key(&name) {
            for el in replace_of.get(&name).unwrap() {
                let definition = fasta::record::Definition::new(&**el, None);
                let record_replace = fasta::Record::new(definition, record.sequence().clone());
                // output the replaced record
                fa_out.write_record(&record_replace)?;
            }
        } else if is_some {
            continue;
        } else {
            // output the original record
            fa_out.write_record(&record)?;
        }
    }

    Ok(())
}

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
        replaces
            .entry(key.clone())
            .or_insert_with(Vec::new)
            .append(&mut others);
    }

    replaces
}
