use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("axt2fas")
        .about("Convert AXT format files to block FA format")
        .after_help(
            r###"
AXT is a format for representing pairwise genomic alignments.
This subcommand converts AXT files into block FA format for further analysis.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- A chromosome sizes file (chr.sizes) for the query genome is required to correctly handle
  coordinates on the negative strand.
- The output file defaults to standard output (stdout). Use the -o option to specify an output file.

Examples:
1. Convert from a file and output to stdout:
   fasr axt2fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt

2. Read from stdin and output to a file:
   cat tests/fasr/example.axt | fasr axt2fas tests/fasr/RM11_1a.chr.sizes stdin -o output.fas

3. Specify target and query names:
   fasr axt2fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt --tname S288c --qname RM11_1a

"###,
        )
        .arg(
            Arg::new("chr.sizes")
                .required(true)
                .index(1)
                .num_args(1)
                .help("Chromosome sizes file for the query genome"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Input AXT file(s) to process"),
        )
        .arg(
            Arg::new("tname")
                .long("tname")
                .num_args(1)
                .default_value("target")
                .help("Target name"),
        )
        .arg(
            Arg::new("qname")
                .long("qname")
                .num_args(1)
                .default_value("query")
                .help("Query name"),
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
    let sizes = intspan::read_sizes(args.get_one::<String>("chr.sizes").unwrap());

    let opt_tname = args.get_one::<String>("tname").unwrap();
    let opt_qname = args.get_one::<String>("qname").unwrap();

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        // Parse each AXT block
        while let Ok(block) = hnsm::next_axt_block(&mut reader, &sizes, opt_tname, opt_qname) {
            for entry in block.entries {
                //----------------------------
                // Output
                //----------------------------
                writer.write_all(entry.to_string().as_ref())?;
            }

            // Add a newline to separate blocks
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
