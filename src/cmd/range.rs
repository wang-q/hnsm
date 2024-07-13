use clap::*;
use std::collections::HashMap;
use std::{ffi, fs, path};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_fasta as fasta;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extract sequences defined by the range(s)")
        .after_help(
            r###"
* <infile> can be plain text or bgzf but not stdin or gzip
* Range - seq_name(strand):start-end
          ^^^^^^^^ required

    Mito
    I:1-100
    I(+):90-150
    S288c.I(-):190-200
    II:21294-22075
    II:23537-24097

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("ranges")
                .required(false)
                .index(2)
                .num_args(0..)
                .help("Ranges of interest"),
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
    let infile = args.get_one::<String>("infile").unwrap();

    let mut fa_out = {
        let writer = intspan::writer(args.get_one::<String>("outfile").unwrap());
        fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_with_writer(writer)
    };

    let is_bgzf = {
        let path = path::Path::new(infile);
        path.extension() == Some(ffi::OsStr::new("gz"))
    };

    let ranges = if args.contains_id("ranges") {
        args.get_many::<String>("ranges")
            .unwrap()
            .cloned()
            .collect()
    } else {
        vec![]
    };

    //----------------------------
    // Open files
    //----------------------------
    let loc_file = format!("{}.loc", infile);
    if !path::Path::new(&loc_file).is_file() {
        hnsm::create_loc(infile, &loc_file, is_bgzf)?;
    }
    let loc_of: HashMap<String, (u64, usize)> = hnsm::load_loc(&loc_file).unwrap();

    let mut reader = if is_bgzf {
        hnsm::Input::Bgzf(
            bgzf::indexed_reader::Builder::default()
                .build_from_path(infile)
                .unwrap(),
        )
    } else {
        hnsm::Input::File(fs::File::open(path::Path::new(infile))?)
    };

    //----------------------------
    // Output
    //----------------------------
    for el in ranges.iter() {
        let rg = intspan::Range::from_str(el);
        let seq_id = rg.chr().to_string();
        if !loc_of.contains_key(&seq_id) {
            eprintln!("{} for [{}] not found in the .loc index file\n", seq_id, el);
            continue;
        }

        let record = hnsm::record_loc(&mut reader, &loc_of, &seq_id)?;

        // name only
        if *rg.start() == 0 {
            fa_out.write_record(&record)?;
            continue;
        }

        let definition = fasta::record::Definition::new(rg.to_string(), None);

        // slice here is 1-based
        let start = Position::new(*rg.start() as usize).unwrap();
        let end = Position::new(*rg.end() as usize).unwrap();

        let record_rg = if rg.strand() == "-" {
            let seq_rc: fasta::record::Sequence = record
                .sequence()
                .complement()
                .rev()
                .collect::<Result<_, _>>()?;
            let slice = seq_rc.slice(start..=end).unwrap();
            fasta::Record::new(definition, slice)
        } else {
            let slice = record.sequence().slice(start..=end).unwrap();
            fasta::Record::new(definition, slice)
        };
        fa_out.write_record(&record_rg)?;
    }

    Ok(())
}

// fn print_type_of<T: ?Sized>(_: &T) {
//     println!("{}", std::any::type_name::<T>())
// }
