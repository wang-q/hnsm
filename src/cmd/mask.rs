use clap::*;
use noodles_fasta as fasta;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mask")
        .about("Soft/hard-masking regions in FA file(s)")
        .after_help(
            r###"
This command masks regions in a FA file based on a runlist JSON file. The runlist specifies
regions to be masked, and the masking can be either soft (lowercase) or hard (replace with N).

The runlist JSON file should have the following format:
    {
        "seq_name": "start1-end1,start2-end2,...",
        ...
    }

Examples:
    1. Soft-mask regions specified in runlist.json:
       hnsm mask input.fa runlist.json -o masked.fa

    2. Hard-mask regions (replace with N):
       hnsm mask input.fa runlist.json --hard -o masked.fa

    3. Output to stdout:
       hnsm mask input.fa runlist.json

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
    let mut fa_in = fasta::io::Reader::new(reader);

    let json = intspan::read_json(args.get_one::<String>("runlist").unwrap());
    let runlists = intspan::json2set(&json);

    let is_hard = args.get_flag("hard");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    //----------------------------
    // Ops
    //----------------------------
    for result in fa_in.records() {
        // obtain record or fail with error
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

            let mut str = seq_out[offset..offset + length].to_string();
            if is_hard {
                str = "N".repeat(length); // Hard-mask with N
            } else {
                str = str.to_lowercase(); // Soft-mask with lowercase
            }
            seq_out.replace_range(offset..offset + length, &str);
        }

        //----------------------------
        // Output
        //----------------------------
        let definition = fasta::record::Definition::new(&*name, None);
        let seq_out = fasta::record::Sequence::from(seq_out.as_bytes().to_vec());
        let record_out = fasta::Record::new(definition, seq_out);
        fa_out.write_record(&record_out)?;
    }

    Ok(())
}
