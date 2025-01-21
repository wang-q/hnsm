use clap::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sixframe")
        .about("Six-Frame Translation")
        .after_help(
            r###"
This command performs six-frame translation on DNA sequences in a FA file.
It translates each sequence into all six possible reading frames (three forward and three reverse)
and identifies open reading frames (ORFs) in the translated protein sequences.

The output includes the ORF sequences along with their positions and frames.

Examples:
1. Perform six-frame translation on a FA file and output to stdout:
   hnsm sixframe input.fa

2. Perform six-frame translation and save the results to a file:
   hnsm sixframe input.fa -o output.fa

3. Combine filters (e.g., ORFs starting with M, ending with *, and at least 50 amino acids):
   hnsm sixframe input.fa --len 50 --start --end

For more information on six-frame translation, visit: https://web.expasy.org/translate/

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input FA file containing DNA sequences"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Minimum length of the amino acid sequence to consider"),
        )
        .arg(
            Arg::new("start")
                .long("start")
                .action(ArgAction::SetTrue)
                .help("Only consider ORFs that start with Methionine (M)"),
        )
        .arg(
            Arg::new("end")
                .long("end")
                .action(ArgAction::SetTrue)
                .help("Only consider ORFs that end with a stop codon (*)"),
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

    let opt_len = *args.get_one::<usize>("len").unwrap();
    let is_start = args.get_flag("start");
    let is_end = args.get_flag("end");

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        // Perform six-frame translation
        let translations = six_frame_translation(&seq[..]);

        // Iterate over each translation frame
        for (protein, frame, is_reverse) in translations {
            // Detect ORFs in the translated protein sequence
            let orfs = hnsm::find_orfs(&protein);

            // Calculate the starting position in the DNA sequence
            let dna_start = if is_reverse {
                seq.len() - frame // Starting position for reverse strand
            } else {
                frame // Starting position for forward strand
            };

            // Adjust dna positions and write each ORF to the output file
            for (orf_seq, start, end) in orfs {
                // Filter ORFs based on the provided options
                if orf_seq.len() < opt_len {
                    continue;
                }
                if is_start && !orf_seq.starts_with('M') {
                    continue;
                }
                if is_end && !orf_seq.ends_with('*') {
                    continue;
                }

                // 1-based
                let orf_start = if is_reverse {
                    dna_start - end * 3 + 1
                } else {
                    dna_start + start * 3 + 1
                };
                let orf_end = if is_reverse {
                    dna_start - start * 3
                } else {
                    dna_start + end * 3
                };

                let header = format!(
                    "{}({}):{}-{}|frame={}",
                    name,
                    if is_reverse { "-" } else { "+" },
                    orf_start,
                    orf_end,
                    frame,
                );
                writer.write_fmt(format_args!(">{}\n{}\n", header, orf_seq))?;
            }
        }
    }

    Ok(())
}

fn six_frame_translation(dna: &[u8]) -> Vec<(String, usize, bool)> {
    let mut translations = Vec::new();

    // Translate the three forward frames
    for frame in 0..3 {
        let frame_dna = &dna[frame..];
        let protein = hnsm::translate(frame_dna);
        translations.push((protein, frame, false)); // false indicates forward strand
    }

    // Translate the three forward frames
    let dna_rc = hnsm::rev_comp(dna).collect::<Vec<_>>();
    for frame in 0..3 {
        let frame_dna = &dna_rc[frame..];
        let protein = hnsm::translate(frame_dna);
        translations.push((protein, frame, true)); // true indicates reverse strand
    }

    translations
}
