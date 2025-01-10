use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filter blocks and optionally format sequences")
        .after_help(
            r###"
This subcommand filters blocks in block FA files based on species name and sequence length.
It can also format sequences by converting them to uppercase or removing dashes.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- If `--name` is not specified, the first species in each block is used as the default.
- Sequences can be filtered based on length using `--ge` (greater than or equal) and `--le` (less than or equal).
- Sequences can be formatted using `--upper` (convert to uppercase) and `--dash` (remove dashes).

Examples:
1. Filter blocks for a specific species:
   fasr filter tests/fasr/example.fas --name S288c

2. Filter blocks with sequences >= 100 bp:
   fasr filter tests/fasr/example.fas --ge 100

3. Filter blocks with sequences <= 200 bp:
   fasr filter tests/fasr/example.fas --le 200

4. Convert sequences to uppercase and remove dashes:
   fasr filter tests/fasr/example.fas --upper --dash

5. Output results to a file:
   fasr filter tests/fasr/example.fas -o output.fas

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
                .help("Filter blocks based on this species"),
        )
        .arg(
            Arg::new("ge")
                .long("ge")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .help("Filter sequences with length >= this value"),
        )
        .arg(
            Arg::new("le")
                .long("le")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .help("Filter sequences with length <= this value"),
        )
        .arg(
            Arg::new("upper")
                .long("upper")
                .action(ArgAction::SetTrue)
                .help("Convert sequences to uppercase"),
        )
        .arg(
            Arg::new("dash")
                .long("dash")
                .action(ArgAction::SetTrue)
                .help("Remove dashes '-'"),
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
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();
    let opt_ge = args.get_one::<usize>("ge").copied().unwrap_or(usize::MAX);
    let opt_le = args.get_one::<usize>("le").copied().unwrap_or(usize::MAX);

    let is_upper = args.get_flag("upper");
    let is_dash = args.get_flag("dash");

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        'BLOCK: while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            // Determine the index of the species
            let idx = if !opt_name.is_empty() {
                if !block.names.contains(opt_name) {
                    continue 'BLOCK;
                }
                block.names.iter().position(|x| x == opt_name).unwrap()
            } else {
                0
            };

            let idx_seq = block.entries[idx].seq();

            // --ge
            if opt_ge != usize::MAX && idx_seq.len() < opt_ge {
                continue 'BLOCK;
            }

            // --le
            if opt_le != usize::MAX && idx_seq.len() > opt_le {
                continue 'BLOCK;
            }

            for entry in &block.entries {
                let mut out_seq: Vec<u8> = vec![];

                for char in entry.seq() {
                    if is_dash && *char == b'-' {
                        continue;
                    }
                    out_seq.push(*char);
                }

                let out_seq = if is_upper {
                    out_seq.to_ascii_uppercase()
                } else {
                    out_seq
                };

                //----------------------------
                // Output
                //----------------------------
                let out_entry = hnsm::FasEntry::from(entry.range(), &out_seq);
                writer.write_all(out_entry.to_string().as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
