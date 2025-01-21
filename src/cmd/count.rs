use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("count")
        .about("Count base statistics in FA file(s)")
        .after_help(
            r###"
This command calculates the base statistics (A, C, G, T, N) for each sequence in one or more FA files.
It outputs a table with the sequence name, length, and counts of each base.

Examples:
1. Count base statistics for a single FA file:
   hnsm count input.fa

2. Count base statistics for multiple FA files:
   hnsm count input1.fa input2.fa.gz

3. Save the output to a file:
   hnsm count input.fa -o output.tsv

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTA file(s) to process"),
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
    // Ops
    //----------------------------
    // Initialize counters
    let mut total_len = 0usize;
    let mut total_base_cnt = [0usize; 5]; // A, C, G, T, N

    // Write the header
    writer.write_fmt(format_args!("#seq\tlen\tA\tC\tG\tT\tN\n"))?;

    // Process each input file
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

        // Process each record
        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq = record.sequence();

            // Count bases in the sequence
            let (len, base_cnt) = count_bases(seq.get(..).unwrap());

            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                name,
                len,
                base_cnt[hnsm::Nt::A as usize],
                base_cnt[hnsm::Nt::C as usize],
                base_cnt[hnsm::Nt::G as usize],
                base_cnt[hnsm::Nt::T as usize],
                base_cnt[hnsm::Nt::N as usize],
            ))?;

            // Update total statistics
            total_len += len;
            for &nt in &[
                hnsm::Nt::A,
                hnsm::Nt::C,
                hnsm::Nt::G,
                hnsm::Nt::T,
                hnsm::Nt::N,
            ] {
                total_base_cnt[nt as usize] += base_cnt[nt as usize];
            }
        }
    }

    writer.write_fmt(format_args!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        "total",
        total_len,
        total_base_cnt[hnsm::Nt::A as usize],
        total_base_cnt[hnsm::Nt::C as usize],
        total_base_cnt[hnsm::Nt::G as usize],
        total_base_cnt[hnsm::Nt::T as usize],
        total_base_cnt[hnsm::Nt::N as usize],
    ))?;

    Ok(())
}

// Count bases in a sequence
fn count_bases(seq: &[u8]) -> (usize, [usize; 5]) {
    let mut len = 0usize;
    let mut base_cnt = [0usize; 5]; // A, C, G, T, N

    for &el in seq {
        let nt = hnsm::to_nt(el);
        if !matches!(nt, hnsm::Nt::Invalid) {
            len += 1;
            base_cnt[nt as usize] += 1;
        }
    }

    (len, base_cnt)
}
