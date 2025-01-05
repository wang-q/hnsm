use clap::*;
use noodles_fasta as fasta;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("masked")
        .about("Identify masked regions in FA file(s")
        .after_help(
            r###"
This command identifies masked regions in one or more FA files. Masked regions can be:
- Lowercase letters
- Regions of N/n

The output is a list of regions in the format:
    seq_name:start-end
    seq_name:position (if start == end)

Examples:
    1. Identify masked regions (lowercase and N/n):
       hnsm masked input.fa -o masked_regions.txt

    2. Identify only N/n gap regions:
       hnsm masked input.fa --gap -o gap_regions.txt

    3. Process multiple input files:
       hnsm masked input1.fa input2.fa -o masked_regions.txt

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FA file(s) to process"),
        )
        .arg(
            Arg::new("gap")
                .long("gap")
                .short('g')
                .action(ArgAction::SetTrue)
                .help("Only identify regions of N/n (gaps)"),
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
    let is_gap = args.get_flag("gap");

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
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

            // Write the last masked region (if any)
            if begin != usize::MAX {
                writer.write_all(out_line(&name, begin, end).as_ref())?;
            }
        }
    }

    Ok(())
}

// Generate the output line for a masked region
fn out_line(name: &str, begin: usize, end: usize) -> String {
    if begin == end {
        format!("{}:{}\n", name, begin + 1)
    } else {
        format!("{}:{}-{}\n", name, begin + 1, end + 1)
    }
}
