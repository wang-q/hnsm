use clap::*;
use itertools::Itertools;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("link")
        .about("Output bi/multi-lateral range links from block FA files")
        .after_help(
            r###"
This subcommand extracts bi/multi-lateral range links from block FA files.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- By default, the tool outputs multi-lateral links (all ranges in a block).
- Use `--pair` to output bilateral (pairwise) links.
- Use `--best` to output best-to-best bilateral links based on sequence similarity.

Examples:
1. Output multi-lateral links:
   fasr link tests/fasr/example.fas

2. Output bilateral (pairwise) links:
   fasr link tests/fasr/example.fas --pair

3. Output best-to-best bilateral links:
   fasr link tests/fasr/example.fas --best

4. Output results to a file:
   fasr link tests/fasr/example.fas -o output.tsv

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("pair")
                .long("pair")
                .action(ArgAction::SetTrue)
                .help("Output bilateral (pairwise) links"),
        )
        .arg(
            Arg::new("best")
                .long("best")
                .action(ArgAction::SetTrue)
                .help("Output best-to-best bilateral links"),
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
    let is_pair = args.get_flag("pair");
    let is_best = args.get_flag("best");

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let headers: Vec<String> = block
                .entries
                .iter()
                .map(|entry| entry.range().to_string())
                .collect();

            //----------------------------
            // Output
            //----------------------------
            if is_pair {
                // Output bilateral (pairwise) links
                for (i, j) in (0..headers.len()).tuple_combinations() {
                    writer.write_all(format!("{}\t{}\n", headers[i], headers[j]).as_ref())?;
                }
            } else if is_best {
                // Output best-to-best bilateral links
                let mut best_pair: Vec<(usize, usize)> = vec![];
                for i in 0..headers.len() {
                    let mut dist_idx: (f32, usize) = (1.0, headers.len() - 1);
                    for j in 0..headers.len() {
                        if i == j {
                            continue;
                        }
                        let dist = hnsm::pair_d(block.entries[i].seq(), block.entries[j].seq());
                        if dist < dist_idx.0 {
                            dist_idx = (dist, j);
                        }
                    }
                    if i < dist_idx.1 {
                        best_pair.push((i, dist_idx.1));
                    } else {
                        best_pair.push((dist_idx.1, i));
                    }
                }
                // // Deduplicate pairs using itertools
                let best_pair: Vec<(usize, usize)> = best_pair.into_iter().unique().collect();

                for (i, j) in best_pair {
                    writer.write_all(format!("{}\t{}\n", headers[i], headers[j]).as_ref())?;
                }
            } else {
                // Output multi-lateral links
                writer.write_all(format!("{}\n", headers.join("\t")).as_ref())?;
            }
        }
    }

    Ok(())
}
