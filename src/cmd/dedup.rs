use clap::*;
use noodles_fasta as fasta;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dedup")
        .about("Deduplicate records in FA file(s)")
        .after_help(
            r###"
* The default behavior is the same as `hnsm filter -u`
* By default, only the forward strand is compared, setting `-b` compares both strands
* `-b` implies `-c`

TODO:
* Remove fully contained sequences

 sequence name
 | |
>sq0 LN:13
     |   |
     description

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("desc")
                .long("desc")
                .short('d')
                .action(ArgAction::SetTrue)
                .help("By name and description"),
        )
        .arg(
            Arg::new("seq")
                .long("seq")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("By sequence"),
        )
        .arg(
            Arg::new("both")
                .long("both")
                .short('b')
                .action(ArgAction::SetTrue)
                .help("Compares both strands"),
        )
        .arg(
            Arg::new("case")
                .long("case")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("Case insensitive"),
        )
        .arg(
            Arg::new("file")
                .long("file")
                .short('f')
                .num_args(1)
                .help("File to save duplicated names"),
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
    let is_desc = args.get_flag("desc");
    let is_seq = args.get_flag("seq");
    let is_both = args.get_flag("both");
    let is_insensitive = args.get_flag("case");

    let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
    let mut fa_out = fasta::io::writer::Builder::default()
        .set_line_base_count(usize::MAX)
        .build_from_writer(writer);

    let mut subject_map: HashMap<u64, Vec<String>> = HashMap::new();
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let name = record.name();
            let desc = record.description();
            let seq = record.sequence();

            let name_str = String::from_utf8(record.name().into()).unwrap();

            // filters
            let mut flag_pass = true;

            // name/desc/sequence to u64 signatures
            let subject = if is_seq {
                if is_both {
                    let fwd = xxhash_rust::xxh3::xxh3_64(&seq[..].to_ascii_uppercase());
                    let rc: fasta::record::Sequence =
                        seq.complement().rev().collect::<Result<_, _>>()?;
                    let rev = xxhash_rust::xxh3::xxh3_64(&rc[..].to_ascii_uppercase());
                    fwd.min(rev)
                } else if is_insensitive {
                    xxhash_rust::xxh3::xxh3_64(&seq[..].to_ascii_uppercase())
                } else {
                    xxhash_rust::xxh3::xxh3_64(&seq[..])
                }
            } else if is_desc && desc.is_some() {
                let full = [name, desc.unwrap()].concat();
                if is_insensitive {
                    xxhash_rust::xxh3::xxh3_64(&full.to_ascii_uppercase())
                } else {
                    xxhash_rust::xxh3::xxh3_64(&full)
                }
            } else {
                if is_insensitive {
                    xxhash_rust::xxh3::xxh3_64(&name.to_ascii_uppercase())
                } else {
                    xxhash_rust::xxh3::xxh3_64(name)
                }
            };

            if subject_map.contains_key(&subject) {
                flag_pass = false;
                subject_map.get_mut(&subject).unwrap().push(name_str);
            } else {
                subject_map.insert(subject, vec![name_str]);
            }

            if !flag_pass {
                continue;
            }
            fa_out.write_record(&record)?;
        }
    }

    if args.contains_id("file") {
        let opt_file = args.get_one::<String>("file").unwrap();
        let mut writer = intspan::writer(opt_file);

        for (_, v) in &subject_map {
            if v.len() < 2 {
                continue;
            }

            for i in 1..v.len() {
                writer.write_fmt(format_args!("{}\t{}\n", v[0], v[i]))?;
            }
        }
    }

    Ok(())
}
