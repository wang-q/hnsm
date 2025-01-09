use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("check")
        .about("Check genome locations in block FA headers")
        .after_help(
            r###"
This tool verifies that the sequences in block FA files match the corresponding locations in a reference genome.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- The reference genome must be provided as a multi-FASTA file.
- `samtools` must be installed and available in $PATH.

Examples:
1. Check all sequences in a block FA file:
   fasr check tests/fasr/NC_000932.fa tests/fasr/A_tha.pair.fas

2. Check sequences for a specific species:
   fasr check tests/fasr/NC_000932.fa tests/fasr/A_tha.pair.fas --name A_tha

"###,
        )
        .arg(
            Arg::new("genome.fa")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to the reference genome FA file"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Input block FA file(s) to check"),
        )
        .arg(
            Arg::new("name")
                .long("name")
                .num_args(1)
                .help("Check sequences for a specific species"),
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
    let opt_genome = args.get_one::<String>("genome.fa").unwrap();
    let opt_name = &args
        .get_one::<String>("name")
        .map(|s| s.as_str())
        .unwrap_or("")
        .to_string();

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let block_names = block.names;

            // Check if a specific species is requested
            if !opt_name.is_empty() && block_names.contains(opt_name) {
                for entry in &block.entries {
                    let entry_name = entry.range().name();
                    if entry_name == opt_name {
                        let status = check_seq(entry, opt_genome)?;
                        writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                    }
                }
            } else if opt_name.is_empty() {
                // Check all sequences in the block
                for entry in &block.entries {
                    let status = check_seq(entry, opt_genome)?;
                    writer.write_all(format!("{}\t{}\n", entry.range(), status).as_ref())?;
                }
            }
        }
    }

    Ok(())
}

fn check_seq(entry: &hnsm::FasEntry, genome: &str) -> anyhow::Result<String> {
    let range = entry.range();
    let seq = if range.strand() == "-" {
        bio::alphabets::dna::revcomp(entry.seq())
    } else {
        entry.seq().to_vec()
    };
    let seq = std::str::from_utf8(&seq)
        .unwrap()
        .to_string()
        .to_ascii_uppercase()
        .replace('-', "");

    let pos = format!("{}:{}-{}", range.chr(), range.start(), range.end());
    let gseq = intspan::get_seq_faidx(genome, &pos)?.to_ascii_uppercase();

    let status = if seq == gseq { "OK" } else { "FAILED" };

    Ok(status.to_string())
}
