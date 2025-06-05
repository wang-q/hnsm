use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("fq2fa")
        .about("Convert FASTQ to FASTA format")
        .after_help(
            r###"
This command converts FASTQ format sequences to FASTA format.

Features:
* Automatic format detection
* Preserves sequence names
* Supports compressed input/output
* Processes multiple input files

Examples:
1. Convert a FASTQ file to FASTA:
   hnsm fq2fa input.fq -o output.fa

2. Convert multiple FASTQ files to a single FASTA:
   hnsm fq2fa input1.fq input2.fq -o output.fa

3. Convert and write to stdout:
   hnsm fq2fa input.fq
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTQ file(s)"),
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
    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut seq_in = noodles_fastq::io::Reader::new(reader);

        for result in seq_in.records() {
            // obtain record or fail with error
            let record = result?;

            // Output FASTA format
            let name = String::from_utf8(record.name().to_vec())?;
            let definition = noodles_fasta::record::Definition::new(name, None);
            let sequence = noodles_fasta::record::Sequence::from(record.sequence().to_vec());
            let record_out = noodles_fasta::Record::new(definition, sequence);
            fa_out.write_record(&record_out)?;
        }
    }

    Ok(())
}
