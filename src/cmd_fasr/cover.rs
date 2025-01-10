use clap::*;
use std::collections::BTreeMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("cover")
        .about("Output covered regions on chromosomes")
        .after_help(
            r###"
This subcommand outputs the coverage of sequences on chromosomes from block FA files.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- The output is in JSON format, showing the coverage of sequences on chromosomes.
- Optionally, you can specify a species name to limit the output to that species.
- For lastz results, use --trim 10

Examples:
1. Calculate coverage for all species:
   fasr cover tests/fasr/example.fas

2. Calculate coverage for a specific species:
   fasr cover tests/fasr/example.fas --name S288c

3. Trim alignment borders to avoid overlaps:
   fasr cover tests/fasr/example.fas --trim 10

4. Output results to a file:
   fasr cover tests/fasr/example.fas -o output.json

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
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("Only output regions for this species"),
        )
        .arg(
            Arg::new("trim")
                .long("trim")
                .num_args(1)
                .value_parser(value_parser!(i32))
                .default_value("0")
                .help("Trim align borders to avoid overlaps"),
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
    let opt_trim = *args.get_one::<i32>("trim").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();

    //----------------------------
    // Ops
    //----------------------------
    let mut res_of: BTreeMap<String, BTreeMap<String, intspan::IntSpan>> = BTreeMap::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let block_names = block.names;

            if !opt_name.is_empty() {
                if !res_of.contains_key(opt_name) {
                    res_of.insert(opt_name.to_string(), BTreeMap::new());
                }
            } else {
                for name in &block_names {
                    if !res_of.contains_key(name) {
                        res_of.insert(name.to_string(), BTreeMap::new());
                    }
                }
            }

            for entry in &block.entries {
                let range = entry.range();
                if !range.is_valid() {
                    continue;
                }

                if !opt_name.is_empty() {
                    if opt_name != range.name() {
                        continue;
                    }
                }

                let res = res_of.get_mut(entry.range().name()).unwrap();

                if !res.contains_key(entry.range().chr()) {
                    res.insert(entry.range().chr().to_string(), intspan::IntSpan::new());
                }

                let intspan = range.intspan().clone().trim(opt_trim);
                res.get_mut(entry.range().chr()).unwrap().merge(&intspan);
            }
        }
    }

    //----------------------------
    // Output
    //----------------------------
    let out_json = if !opt_name.is_empty() {
        // Output coverage for a single species
        intspan::set2json(res_of.first_key_value().unwrap().1)
    } else {
        // Output coverage for all species
        intspan::set2json_m(&res_of)
    };
    // Write the JSON output to the specified file or stdout
    intspan::write_json(args.get_one::<String>("outfile").unwrap(), &out_json)?;

    Ok(())
}
