use clap::*;
use noodles_fastq::Record;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("interleave")
        .about("Interleave two PE files")
        .after_help(
            r###"
* The number of input files can be one or two.
  When the number is one, output an `N` for the second sequence.
* When using the `--fq` option, all arbitrary quality values will be set to ! (ASCII 33)
* Multiple reads of the input file are required, so reading from stdin is not supported

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..=2)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("fq")
                .long("fq")
                .action(ArgAction::SetTrue)
                .help("Write FQ"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .num_args(1)
                .default_value("read")
                .help("Prefix of record names"),
        )
        .arg(
            Arg::new("start")
                .long("start")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Starting index"),
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

    let is_out_fq = args.get_flag("fq");
    let opt_prefix = args.get_one::<String>("prefix").unwrap();
    let mut opt_start = *args.get_one::<usize>("start").unwrap();

    let infiles = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>();
    let is_in_fq = hnsm::is_fq(infiles[0]);

    //----------------------------
    // Ops
    //----------------------------
    if infiles.len() == 1 {
        if is_in_fq {
            let reader = intspan::reader(infiles[0]);
            let mut seq_in = noodles_fastq::io::Reader::new(reader);
            for result in seq_in.records() {
                // obtain record or fail with error
                let record = result?;

                if is_out_fq {
                    // Output FASTQ format
                    // R1
                    write_fq(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        record.sequence(),
                        record.quality_scores(),
                    )?;

                    // R2
                    write_fq(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        b"N",
                        b"!",
                    )?;
                } else {
                    // Output FASTA format
                    // R1
                    write_fa(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        record.sequence(),
                    )?;

                    // R2
                    write_fa(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        b"\n",
                    )?;
                }

                opt_start += 1;
            }
        } else {
            let reader = intspan::reader(infiles[0]);
            let mut seq_in = noodles_fasta::io::Reader::new(reader);
            for result in seq_in.records() {
                // obtain record or fail with error
                let record = result?;

                if is_out_fq {
                    // Output FASTQ format
                    // R1
                    write_fq(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        &record.sequence()[..],
                        &vec![b'!'; record.sequence().len()],
                    )?;

                    // R2
                    write_fq(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        b"N",
                        b"!",
                    )?;
                } else {
                    // Output FASTA format
                    // R1
                    write_fa(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        &record.sequence()[..],
                    )?;

                    // R2
                    write_fa(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        b"\n",
                    )?;
                }

                opt_start += 1;
            }
        }
    } else {
        if is_in_fq {
            let reader = intspan::reader(infiles[0]);
            let mut seq1_in = noodles_fastq::io::Reader::new(reader);
            let reader = intspan::reader(infiles[1]);
            let mut seq2_in = noodles_fastq::io::Reader::new(reader);

            let mut zipped = std::iter::zip(seq1_in.records(), seq2_in.records());

            for (result1, result2) in zipped {
                // obtain record or fail with error
                let record1 = result1?;
                let record2 = result2?;

                if is_out_fq {
                    // Output FASTQ format
                    // R1
                    write_fq(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        record1.sequence(),
                        record1.quality_scores(),
                    )?;

                    // R2
                    write_fq(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        record2.sequence(),
                        record2.quality_scores(),
                    )?;
                } else {
                    // Output FASTA format
                    // R1
                    write_fa(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        record1.sequence(),
                    )?;

                    // R2
                    write_fa(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        record2.sequence(),
                    )?;
                }

                opt_start += 1;
            }
        } else {
            let reader = intspan::reader(infiles[0]);
            let mut seq1_in = noodles_fasta::io::Reader::new(reader);
            let reader = intspan::reader(infiles[1]);
            let mut seq2_in = noodles_fasta::io::Reader::new(reader);

            let mut zipped = std::iter::zip(seq1_in.records(), seq2_in.records());

            for (result1, result2) in zipped {
                // obtain record or fail with error
                let record1 = result1?;
                let record2 = result2?;

                if is_out_fq {
                    // Output FASTQ format
                    // R1
                    write_fq(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        &record1.sequence()[..],
                        &vec![b'!'; record1.sequence().len()],
                    )?;

                    // R2
                    write_fq(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        &record2.sequence()[..],
                        &vec![b'!'; record2.sequence().len()],
                    )?;
                } else {
                    // Output FASTA format
                    // R1
                    write_fa(
                        &mut writer,
                        &format!("{}{}/1", opt_prefix, opt_start),
                        &record1.sequence()[..],
                    )?;

                    // R2
                    write_fa(
                        &mut writer,
                        &format!("{}{}/2", opt_prefix, opt_start),
                        &record2.sequence()[..],
                    )?;
                }

                opt_start += 1;
            }
        }
    }

    Ok(())
}

fn write_fq(
    writer: &mut Box<dyn Write>,
    seq_name: &str,
    seq: &[u8],
    qual: &[u8],
) -> Result<(), Error> {
    writer.write_fmt(format_args!("@{}\n", seq_name))?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    writer.write_all(b"+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_fa(writer: &mut Box<dyn Write>, seq_name: &str, seq: &[u8]) -> Result<(), Error> {
    writer.write_fmt(format_args!(">{}\n", seq_name))?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}
