use clap::*;
use noodles_fasta as fasta;
use std::collections::BTreeSet;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filter records in FA file(s)")
        .after_help(
            r###"
* Not all faFilter options have been implemented
  Wildcards for names can be easily implemented with `hnsm some`
* This subcommand is also a formatter
    * -l is used to set the number of bases per line
    * -b/--block is not implemented here

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("minsize")
                .long("minsize")
                .short('a')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences at least this big ('a'-smallest)"),
        )
        .arg(
            Arg::new("maxsize")
                .long("maxsize")
                .short('z')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences this size or smaller ('z'-biggest)"),
        )
        .arg(
            Arg::new("maxn")
                .long("maxn")
                .short('n')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Pass sequences with fewer than this number of Ns"),
        )
        .arg(
            Arg::new("uniq")
                .long("uniq")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Unique, removes duplicated ids, keeping the first"),
        )
        .arg(
            Arg::new("upper")
                .long("upper")
                .short('U')
                .action(ArgAction::SetTrue)
                .help("Convert all sequences to upper cases"),
        )
        .arg(
            Arg::new("iupac")
                .long("iupac")
                .short('N')
                .action(ArgAction::SetTrue)
                .help("Convert IUPAC ambiguous codes to 'N'"),
        )
        .arg(
            Arg::new("dash")
                .long("dash")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("Remove dashes '-'"),
        )
        .arg(
            Arg::new("simplify")
                .long("simplify")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Simplify sequence names"),
        )
        .arg(
            Arg::new("line")
                .long("line")
                .short('l')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Sequence line length"),
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
    let opt_minsize = if args.contains_id("minsize") {
        *args.get_one::<usize>("minsize").unwrap()
    } else {
        usize::MAX
    };
    let opt_maxsize = if args.contains_id("maxsize") {
        *args.get_one::<usize>("maxsize").unwrap()
    } else {
        usize::MAX
    };
    let opt_maxn = if args.contains_id("maxn") {
        *args.get_one::<usize>("maxn").unwrap()
    } else {
        usize::MAX
    };
    let opt_line = if args.contains_id("line") {
        *args.get_one::<usize>("line").unwrap()
    } else {
        usize::MAX
    };

    let is_uniq = args.get_flag("uniq");
    let is_upper = args.get_flag("upper");
    let is_iupac = args.get_flag("iupac");
    let is_dash = args.get_flag("dash");
    let is_simplify = args.get_flag("simplify");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(opt_line)
        .build_with_writer(writer);

    let mut set_list: BTreeSet<String> = BTreeSet::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let mut name = String::from_utf8(record.name().into()).unwrap();
            if is_simplify {
                if let Some(i) = name.find(&[' ', '.', ',', '-'][..]) {
                    name = name[..i].to_string();
                }
            }
            let seq = record.sequence();

            // filters
            let mut flag_pass = true;
            if opt_minsize != usize::MAX && seq.len() < opt_minsize {
                flag_pass = false;
            } else if opt_maxsize != usize::MAX && seq.len() > opt_maxsize {
                flag_pass = false;
            } else if opt_maxn != usize::MAX {
                if hnsm::count_n(seq.get(..).unwrap()) > opt_maxn {
                    flag_pass = false;
                }
            } else if is_uniq {
                // If the set did not previously contain an equal value, true is returned.
                let seen = !set_list.insert(name.clone());
                if seen {
                    flag_pass = false;
                }
            }

            if !flag_pass {
                continue;
            }

            // formatters
            let mut seq_out = String::new();
            for nt in seq.get(..).unwrap().iter() {
                if is_dash && *nt == b'-' {
                        continue;

                }
                if is_iupac {
                    if is_upper {
                        seq_out.push(char::from(hnsm::to_n(*nt)).to_ascii_uppercase());
                    } else {
                        seq_out.push(char::from(hnsm::to_n(*nt)));
                    }
                } else {
                    if is_upper {
                        seq_out.push(char::from(*nt).to_ascii_uppercase());
                    } else {
                        seq_out.push(char::from(*nt));
                    }
                }
            } // end of each nt

            let definition = fasta::record::Definition::new(&*name, None);
            let seq_out = fasta::record::Sequence::from(seq_out.as_bytes().to_vec());
            let record_out = fasta::Record::new(definition, seq_out);
            fa_out.write_record(&record_out)?;
        }
    }

    Ok(())
}
