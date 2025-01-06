use clap::*;
use noodles::fasta;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sixframe")
        .about("Six-Frame Translation")
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
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
    //----------------------------
    // Args
    //----------------------------
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = fasta::io::Reader::new(reader);

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
