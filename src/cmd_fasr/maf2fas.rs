use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("maf2fas")
        .about("Convert maf to block fasta")
        .after_help(
            r###"
* <infiles> are paths to maf files, .maf.gz is supported
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

    //----------------------------
    // Operating
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_maf_block(&mut reader) {
            // Can't use reference as entry.alignment does not Copy
            for entry in block.entries {
                let range = entry.to_range();
                let seq = String::from_utf8(entry.alignment).unwrap();

                //----------------------------
                // Output
                //----------------------------
                writer.write_all(format!(">{}\n{}\n", range, seq).as_ref())?;
            }

            // end of a block
            writer.write_all("\n".as_ref())?;
        }
    }

    Ok(())
}
