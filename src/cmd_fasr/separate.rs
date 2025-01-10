use clap::*;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("separate")
        .about("Separate block FA files by species")
        .after_help(
            r###"
This subcommand separates block FA files by species, creating individual output files for each species.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- Dashes ('-') in sequences are removed.
- If the target file already exists, it will be overwritten.
- Optionally, sequences can be reverse-complemented if the chromosome strand is '-'.

Examples:
1. Separate block FA files by species:
   fasr separate tests/fasr/example.fas -o output_dir

2. Separate block FA files and reverse-complement sequences:
   fasr separate tests/fasr/example.fas -o output_dir --rc

3. Use a custom suffix for output files:
   fasr separate tests/fasr/example.fas -o output_dir --suffix .fa

4. Output to stdout:
   fasr separate tests/fasr/example.fas

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
            Arg::new("suffix")
                .long("suffix")
                .short('s')
                .num_args(1)
                .default_value(".fasta")
                .help("File extension for output files"),
        )
        .arg(
            Arg::new("rc")
                .long("rc")
                .action(ArgAction::SetTrue)
                .help("Reverse-complement sequences when chromosome strand is '-'"),
        )
        .arg(
            Arg::new("outdir")
                .short('o')
                .long("outdir")
                .num_args(1)
                .default_value("stdout")
                .help("Output location. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_suffix = args.get_one::<String>("suffix").unwrap();
    let is_rc = args.get_flag("rc");

    //----------------------------
    // Ops
    //----------------------------
    let mut file_of: BTreeMap<String, std::fs::File> = BTreeMap::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            for entry in &block.entries {
                let entry_name = entry.range().name(); // Don't borrow the following `range`
                let mut range = entry.range().clone();

                // Reverse-complement the sequence if needed
                let seq = if is_rc && range.strand() == "-" {
                    *range.strand_mut() = "+".to_string();
                    bio::alphabets::dna::revcomp(entry.seq())
                } else {
                    entry.seq().to_vec()
                };

                // Remove dashes from the sequence
                let seq = std::str::from_utf8(&seq)?.to_string().replace('-', "");

                //----------------------------
                // Output
                //----------------------------
                if outdir == "stdout" {
                    print!(">{}\n{}\n", range, seq);
                } else {
                    if !file_of.contains_key(entry_name) {
                        let path =
                            std::path::Path::new(outdir).join(range.name().to_owned() + opt_suffix);
                        let file = std::fs::OpenOptions::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(path)?;
                        file_of.insert(entry_name.to_string(), file);
                    }
                    write!(file_of.get(entry_name).unwrap(), ">{}\n{}\n", range, seq)?;
                }
            }
        }
    }

    Ok(())
}
