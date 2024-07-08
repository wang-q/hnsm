use clap::*;
use hnsm::Nt;
use noodles_fasta as fasta;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("count")
        .about("Count base statistics in FA file(s)")
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input file to use"),
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
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let mut total_len = 0usize;
    let mut total_base_cnt = [0usize; 5];

    writer.write_fmt(format_args!("#seq\tlen\tA\tC\tG\tT\tN\n"))?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;
            let name = String::from_utf8(record.name().into()).unwrap();
            let seq = record.sequence();

            let mut len = 0usize;
            let mut base_cnt = [0usize; 5];
            for el in seq.get(..).unwrap().iter() {
                let nt = hnsm::to_nt(*el);
                if !matches!(nt, Nt::Invalid) {
                    len += 1;
                    base_cnt[nt as usize] += 1;
                }
            }

            writer.write_fmt(format_args!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                name,
                len,
                base_cnt[Nt::A as usize],
                base_cnt[Nt::C as usize],
                base_cnt[Nt::G as usize],
                base_cnt[Nt::T as usize],
                base_cnt[Nt::N as usize],
            ))?;

            total_len += len;
            for &nt in &[Nt::A, Nt::C, Nt::G, Nt::T, Nt::N] {
                total_base_cnt[nt as usize] += base_cnt[nt as usize];
            }
        }
    }

    writer.write_fmt(format_args!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        "total",
        total_len,
        total_base_cnt[Nt::A as usize],
        total_base_cnt[Nt::C as usize],
        total_base_cnt[Nt::G as usize],
        total_base_cnt[Nt::T as usize],
        total_base_cnt[Nt::N as usize],
    ))?;

    Ok(())
}
