use clap::*;
use std::collections::BTreeSet;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("filter")
        .about("Filter and format sequences in FA file(s)")
        .after_help(
            r###"
This command filters and formats sequences in FA files.

Filters:
* --minsize N: Keep sequences >= N bp
* --maxsize N: Keep sequences <= N bp
* --maxn N: Keep sequences with < N ambiguous bases
* --uniq: Remove duplicate sequence IDs

Formatters:
* --upper: Convert sequences to uppercase
* --iupac: Convert ambiguous codes to 'N'
* --dash: Remove dashes from sequences
* --simplify: Simplify sequence names (truncate at first space/./,/-)
* --line N: Set sequence line length

Notes:
* Multiple filters can be combined
* Supports both plain text and gzipped (.gz) files
* For duplicate IDs, keeps the first occurrence
* Not all faFilter options have been implemented
  Wildcards for names can be easily implemented with `hnsm some`

Examples:
1. Filter by size:
   hnsm filter input.fa --minsize 100 --maxsize 1000

2. Format sequences:
   hnsm filter input.fa --upper --iupac --line 80

3. Process multiple files:
   hnsm filter *.fa --uniq --simplify -o output.fa

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
    //----------------------------
    // Args
    //----------------------------
    let opt_minsize = args
        .get_one::<usize>("minsize")
        .copied()
        .unwrap_or(usize::MAX);
    let opt_maxsize = args
        .get_one::<usize>("maxsize")
        .copied()
        .unwrap_or(usize::MAX);
    let opt_maxn = args.get_one::<usize>("maxn").copied().unwrap_or(usize::MAX);
    let opt_line = args.get_one::<usize>("line").copied().unwrap_or(usize::MAX);

    let is_uniq = args.get_flag("uniq");
    let is_upper = args.get_flag("upper");
    let is_iupac = args.get_flag("iupac");
    let is_dash = args.get_flag("dash");
    let is_simplify = args.get_flag("simplify");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = noodles_fasta::io::writer::Builder::default()
        .set_line_base_count(opt_line)
        .build_from_writer(writer);

    //----------------------------
    // Process
    //----------------------------
    let mut set_list: BTreeSet<String> = BTreeSet::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let mut name = String::from_utf8(record.name().into())?;
            if is_simplify {
                if let Some(i) = name.find(&[' ', '.', ',', '-'][..]) {
                    name = name[..i].to_string();
                }
            }
            let seq = record.sequence();

            // Apply filters
            if !pass_filters(
                seq,
                opt_minsize,
                opt_maxsize,
                opt_maxn,
                is_uniq,
                &mut set_list,
                &name,
            ) {
                continue;
            }

            // Apply formatters
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

            let definition = noodles_fasta::record::Definition::new(&*name, None);
            let seq_out = noodles_fasta::record::Sequence::from(seq_out.as_bytes().to_vec());
            let record_out = noodles_fasta::Record::new(definition, seq_out);
            fa_out.write_record(&record_out)?;
        }
    }

    Ok(())
}

// Check if a sequence passes all filters
fn pass_filters(
    seq: &noodles_fasta::record::Sequence,
    minsize: usize,
    maxsize: usize,
    maxn: usize,
    is_uniq: bool,
    set_list: &mut BTreeSet<String>,
    name: &str,
) -> bool {
    if minsize != usize::MAX && seq.len() < minsize {
        return false;
    }
    if maxsize != usize::MAX && seq.len() > maxsize {
        return false;
    }
    if maxn != usize::MAX && hnsm::count_n(seq.get(..).unwrap()) > maxn {
        return false;
    }
    // If the set did not previously contain an equal value, true is returned
    if is_uniq && !set_list.insert(name.to_string()) {
        return false;
    }
    true
}
