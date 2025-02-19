use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mask")
        .about("Mask regions in FA file(s)")
        .after_help(
            r###"
This command masks specified regions in FASTA sequences.

Masking modes:
* Soft-masking (default): Convert to lowercase
* Hard-masking (--hard): Replace with N's

Input format (runlist.json):
{
    "seq1": "1-100,200-300",    # Mask positions 1-100 and 200-300
    "seq2": "50-150",           # Mask positions 50-150
    "seq3": "1-50,90-100,..."   # Multiple regions allowed
}

Notes:
* 1-based coordinates
* Inclusive ranges
* Sequences not in runlist remain unchanged
* Supports both plain text and gzipped (.gz) files
* Invalid ranges are silently ignored

Examples:
1. Soft-mask regions:
   hnsm mask input.fa regions.json -o output.fa

2. Hard-mask regions:
   hnsm mask input.fa regions.json --hard -o output.fa

3. Process gzipped files:
   hnsm mask input.fa.gz regions.json -o output.fa.gz

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input FA file to process"),
        )
        .arg(
            Arg::new("runlist")
                .required(true)
                .num_args(1)
                .index(2)
                .help("JSON file specifying regions to mask"),
        )
        .arg(
            Arg::new("hard")
                .long("hard")
                .action(ArgAction::SetTrue)
                .help("Hard-mask regions (replace with N)"),
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
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = noodles_fasta::io::Reader::new(reader);

    let json = intspan::read_json(args.get_one::<String>("runlist").unwrap());
    let runlists = intspan::json2set(&json);

    let is_hard = args.get_flag("hard");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    //----------------------------
    // Process
    //----------------------------
    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        if !runlists.contains_key(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        // Get the regions to mask for this sequence
        let ints = runlists.get(&name).unwrap();
        let mut seq_out = String::from_utf8(seq[..].into())?;

        for (lower, upper) in ints.spans().iter() {
            let offset = (lower - 1) as usize;
            let length = (upper - lower + 1) as usize;

            let str = if is_hard {
                "N".repeat(length)
            } else {
                seq_out[offset..offset + length].to_lowercase()
            };
            seq_out.replace_range(offset..offset + length, &str);
        }

        //----------------------------
        // Output
        //----------------------------
        let definition = noodles_fasta::record::Definition::new(&*name, None);
        let seq_out = noodles_fasta::record::Sequence::from(seq_out.as_bytes().to_vec());
        let record_out = noodles_fasta::Record::new(definition, seq_out);
        fa_out.write_record(&record_out)?;
    }

    Ok(())
}
