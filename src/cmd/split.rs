use clap::*;
use noodles_fasta as fasta;
use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Split FA file(s) into several files")
        .after_help(
            r#"
Modes

* name  - using sequence names as file names
* about - about `count` bytes each by record
    * -c, -e, -m

"#,
        )
        .arg(
            Arg::new("mode")
                .required(true)
                .index(1)
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("name"),
                    builder::PossibleValue::new("about"),
                ])
                .help("Set the mode"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("count")
                .long("count")
                .short('c')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("bytes "),
        )
        .arg(
            Arg::new("even")
                .long("even")
                .short('e')
                .action(ArgAction::SetTrue)
                .help("Record number in one file should be EVEN"),
        )
        .arg(
            Arg::new("maxpart")
                .long("maxpart")
                .short('m')
                .num_args(1)
                .default_value("999")
                .value_parser(value_parser!(usize))
                .help("Max parts"),
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
    let mode = args.get_one::<String>("mode").unwrap();

    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        fs::create_dir_all(outdir)?;
    }

    let mut fh_of: BTreeMap<String, File> = BTreeMap::new();

    //----------------------------
    // Operating
    //----------------------------
    if mode == "name" {
        for infile in args.get_many::<String>("infiles").unwrap() {
            let reader = intspan::reader(infile);
            let mut fa_in = fasta::io::Reader::new(reader);

            for result in fa_in.records() {
                // obtain record or fail with error
                let record = result?;

                let name = String::from_utf8(record.name().into()).unwrap();
                let seq = record.sequence();
                let seq_str = String::from_utf8(seq.get(..).unwrap().to_vec()).unwrap();

                //----------------------------
                // Output
                //----------------------------
                if outdir == "stdout" {
                    print!(">{}\n{}\n", name, seq_str);
                } else {
                    let filename = name
                        .clone()
                        .replace(['(', ')', ':'], "_")
                        .replace("__", "_");
                    gen_fh(outdir, &mut fh_of, &filename)?;
                    write!(fh_of.get(&filename).unwrap(), ">{}\n{}\n", name, seq_str)?;
                }
            }
        }
    } else if mode == "about" {
        let opt_count = if args.contains_id("count") {
            *args.get_one::<usize>("count").unwrap()
        } else {
            usize::MAX
        };
        let is_even = args.get_flag("even");
        let opt_maxpart = *args.get_one::<usize>("maxpart").unwrap();

        let mut cur_cnt = 0;
        let mut record_sn = 0;
        let mut file_sn = 0;
        let part_width = (opt_maxpart.checked_ilog10().unwrap_or(0) + 1) as usize;

        'outer: for infile in args.get_many::<String>("infiles").unwrap() {
            let reader = intspan::reader(infile);
            let mut fa_in = fasta::io::Reader::new(reader);

            for result in fa_in.records() {
                if file_sn > opt_maxpart {
                    break 'outer;
                }

                // obtain record or fail with error
                let record = result?;

                let name = String::from_utf8(record.name().into()).unwrap();

                let seq = record.sequence();
                let seq_str = String::from_utf8(seq.get(..).unwrap().to_vec()).unwrap();

                //----------------------------
                // Output
                //----------------------------
                if outdir == "stdout" {
                    print!(">{}\n{}\n", name, seq_str);
                } else {
                    let filename = format!("{:0width$}", file_sn, width = part_width);
                    gen_fh(outdir, &mut fh_of, &filename)?;
                    write!(fh_of.get(&filename).unwrap(), ">{}\n{}\n", name, seq_str)?;
                }
                cur_cnt += seq.len();
                record_sn += 1;

                if is_even {
                    if record_sn % 2 != 0 {
                        continue;
                    }
                } else if cur_cnt > opt_count {
                    cur_cnt = 0;
                    record_sn = 0;
                    file_sn += 1;
                }
            } // record
        } // file
    }

    Ok(())
}

fn gen_fh(
    outdir: &String,
    fh_of: &mut BTreeMap<String, File>,
    filename: &String,
) -> Result<(), Error> {
    if !fh_of.contains_key(filename) {
        let path = Path::new(outdir).join(filename.clone() + ".fa");
        let file = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        fh_of.insert(filename.clone(), file);
    }
    Ok(())
}
