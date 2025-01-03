use clap::*;
use noodles_fasta as fasta;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mask")
        .about("Soft/hard-masking regions in FA file(s)")
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("runlist")
                .required(true)
                .num_args(1)
                .index(2)
                .help("The runlist json file"),
        )
        .arg(
            Arg::new("hard")
                .long("hard")
                .action(ArgAction::SetTrue)
                .help("Change masked regions to N"),
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

        let name = String::from_utf8(record.name().into()).unwrap();
        let seq = record.sequence();

        if !runlists.contains_key(&name) {
            fa_out.write_record(&record)?;
            continue;
        }

        let ints = runlists.get(&name).unwrap();
        let mut seq_out = String::from_utf8(seq[..].into()).unwrap();

        for (lower, upper) in ints.spans().iter() {
            let offset = (lower - 1) as usize;
            let length = (upper - lower + 1) as usize;

            let mut str = seq_out[offset..offset + length].to_string();
            if is_hard {
                str = "N".repeat(length);
            } else {
                str = str.to_lowercase();
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
