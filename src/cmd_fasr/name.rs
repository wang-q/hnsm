use clap::*;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("name")
        .about("Output all species names from block FA files")
        .after_help(
            r###"
This subcommand extracts and outputs all species names from block FA files.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- By default, the subcommand outputs a list of unique species names.
- Use `--count` to also output the number of occurrences of each species name.

Examples:
1. Output all species names:
   fasr name tests/fasr/example.fas

2. Output species names with occurrence counts:
   fasr name tests/fasr/example.fas --count

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("count")
                .long("count")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Output species names with occurrence counts"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let is_count = args.get_flag("count");

    //----------------------------
    // Operating
    //----------------------------
    let mut names: Vec<String> = vec![];
    let mut count_of: BTreeMap<String, i32> = BTreeMap::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            for entry in &block.entries {
                let range = entry.range();

                // Collect unique species names
                if !names.contains(range.name()) {
                    names.push(range.name().to_string());
                }

                // Count occurrences of each species name
                count_of
                    .entry(range.name().to_string())
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
        }
    }

    //----------------------------
    // Output
    //----------------------------
    for name in &names {
        if is_count {
            let value = count_of.get(name).unwrap();
            writer.write_all(format!("{}\t{}\n", name, value).as_ref())?;
        } else {
            writer.write_all(format!("{}\n", name).as_ref())?;
        }
    }

    Ok(())
}
