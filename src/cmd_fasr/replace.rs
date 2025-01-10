use clap::*;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("replace")
        .about("Replace headers in block FA files")
        .after_help(
            r###"
This subcommand replaces headers (name.range) in block FA files based on a TSV file.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- The TSV file should contain one or more fields:
  - `original_name  replace_name   more_replace_name`
  - One field: Deletes the entire alignment block for the specified species.
  - Three or more fields: Duplicates the entire alignment block for each replacement name.
- Does not support replacing multiple records in one block.

Examples:
1. Replace species names in a block FA file:
   fasr replace tests/fasr/replace.tsv tests/fasr/example.fas

"###,
        )
        .arg(
            Arg::new("replace.tsv")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to the TSV file containing replacement rules"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Input block FA file(s) to process"),
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

    let mut replace_of: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in intspan::read_lines(args.get_one::<String>("replace.tsv").unwrap()) {
        let rgs: Vec<_> = line.split('\t').collect();

        if rgs.is_empty() {
            continue;
        }

        let rg = rgs.first().unwrap().to_string();
        let replacements = rgs
            .iter()
            .skip(1)
            .map(|e| e.to_string())
            .collect::<Vec<String>>();
        replace_of.insert(rg.to_string(), replacements);
    }

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let originals = block.headers.clone();

            let matched: Vec<String> = replace_of
                .keys()
                .filter(|e| originals.contains(*e))
                .map(|e| e.to_string())
                .collect();

            // eprintln!("originals = {:#?}", originals);
            // eprintln!("matched = {:#?}", matched);

            if matched.is_empty() || matched.len() > 1 {
                // block untouched

                if matched.len() > 1 {
                    eprintln!("Doesn't support replacing multiple records in one block");
                }

                //----------------------------
                // Output
                //----------------------------
                for entry in &block.entries {
                    writer.write_all(entry.to_string().as_ref())?;
                }
                writer.write_all("\n".as_ref())?;
            } else {
                let original = matched.first().unwrap();
                let idx = block.headers.iter().position(|e| e == original).unwrap();
                for new in replace_of.get(original).unwrap() {
                    for (i, entry) in block.entries.iter().enumerate() {
                        //----------------------------
                        // Output
                        //----------------------------
                        if i == idx {
                            writer.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    new,
                                    String::from_utf8(entry.seq().to_vec()).unwrap()
                                )
                                .as_ref(),
                            )?;
                        } else {
                            writer.write_all(entry.to_string().as_ref())?;
                        }
                    }

                    writer.write_all("\n".as_ref())?;
                }
            }
        }
    }

    Ok(())
}
