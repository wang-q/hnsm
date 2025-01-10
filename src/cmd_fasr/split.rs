use clap::*;
use std::collections::BTreeMap;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Split block FA files into per-alignment or per-chromosome FA files")
        .after_help(
            r###"
This subcommand splits block FA files into individual FA files, either per alignment or per chromosome.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- By default, each alignment block is written to a separate file.
- Use `--chr` to split files by chromosome.
- Use `--simple` to simplify headers by keeping only species names.

Examples:
1. Split block FA files into per-alignment files:
   fasr split tests/fasr/example.fas -o output_dir

2. Split block FA files into per-chromosome files:
   fasr split tests/fasr/example.fas -o output_dir --chr

3. Simplify headers in output files:
   fasr split tests/fasr/example.fas -o output_dir --simple

4. Use a custom suffix for output files:
   fasr split tests/fasr/example.fas -o output_dir --suffix .fa

5. Output to stdout:
   fasr split tests/fasr/example.fas

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
            Arg::new("suffix")
                .long("suffix")
                .short('s')
                .num_args(1)
                .default_value(".fas")
                .help("File extension for output files"),
        )
        .arg(
            Arg::new("chr")
                .long("chr")
                .action(ArgAction::SetTrue)
                .help("Split files by chromosomes"),
        )
        .arg(
            Arg::new("simple")
                .long("simple")
                .action(ArgAction::SetTrue)
                .help("Simplify headers by keeping only species names"),
        )
        .arg(
            Arg::new("outdir")
                .short('o')
                .long("outdir")
                .num_args(1)
                .default_value("stdout")
                .help("Output location. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_suffix = args.get_one::<String>("suffix").unwrap();
    let is_chr = args.get_flag("chr");
    let is_simple = args.get_flag("simple");

    let mut file_of: BTreeMap<String, std::fs::File> = BTreeMap::new();

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = intspan::reader(infile);

        while let Ok(block) = hnsm::next_fas_block(&mut reader) {
            let filename = if is_chr {
                let tname = block.entries.first().unwrap().range().name();
                let tchr = block.entries.first().unwrap().range().chr();
                format!("{}.{}", tname, tchr)
            } else {
                let trange = &block.entries.first().unwrap().range().clone();
                trange.to_string()
            }
            .replace(['(', ')', ':'], "_")
            .replace("__", "_");

            for entry in &block.entries {
                let range = entry.range().clone();
                let seq = std::str::from_utf8(entry.seq()).unwrap();

                //----------------------------
                // Output
                //----------------------------
                if outdir == "stdout" {
                    let header = if is_simple {
                        range.name().to_string()
                    } else {
                        range.to_string()
                    };
                    print!(">{}\n{}\n", header, seq);
                } else {
                    if !file_of.contains_key(&filename) {
                        let path = std::path::Path::new(outdir).join(filename.clone() + opt_suffix);
                        let file = std::fs::OpenOptions::new()
                            .create(true)
                            .write(true)
                            .truncate(true)
                            .open(path)?;
                        file_of.insert(filename.clone(), file);
                    }
                    write!(file_of.get(&filename).unwrap(), ">{}\n{}\n", range, seq)?;
                }
            }

            // end of a block
            if outdir == "stdout" {
                println!();
            } else {
                writeln!(file_of.get(&filename).unwrap())?;
            }
        }
    }

    Ok(())
}
