use clap::*;
use noodles_fasta as fasta;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("masked")
        .about("Masked regions in FA file(s)")
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("gap")
                .long("gap")
                .short('g')
                .action(ArgAction::SetTrue)
                .help("Only regions of N/n"),
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
    let is_gap = args.get_flag("gap");

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let name = String::from_utf8(record.name().into()).unwrap();
            let seq = record.sequence();

            let mut begin = usize::MAX;
            let mut end = usize::MAX;

            for (i, el) in seq.get(..).unwrap().iter().enumerate() {
                let is_masked = if is_gap {
                    hnsm::is_n(*el)
                } else {
                    hnsm::is_n(*el) || hnsm::is_lower(*el)
                };
                if is_masked {
                    if begin == usize::MAX {
                        begin = i;
                        end = i;
                    } else {
                        end = i;
                    }
                } else if begin != usize::MAX {
                    writer.write_all(out_line(&name, begin, end).as_ref())?;

                    // reset
                    begin = usize::MAX;
                    end = usize::MAX;
                }
            }

            // last region
            if begin != usize::MAX {
                writer.write_all(out_line(&name, begin, end).as_ref())?;
            }
        }
    }

    Ok(())
}

fn out_line(name: &str, begin: usize, end: usize) -> String {
    if begin == end {
        format!("{}:{}\n", name, begin + 1)
    } else {
        format!("{}:{}-{}\n", name, begin + 1, end + 1)
    }
}
