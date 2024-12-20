use clap::*;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("name")
        .about("Output all species names")
        .after_help(
            r###"
* <infiles> are paths to block fasta files, .fas.gz is supported
* infile == stdin means reading from STDIN

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("count")
                .long("count")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Also count name occurrences"),
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

    let mut names: Vec<String> = vec![];
    let mut count_of: BTreeMap<String, i32> = BTreeMap::new();

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            for entry in &block.entries {
                let range = entry.range();

                if !names.contains(range.name()) {
                    names.push(range.name().to_string());
                }

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
