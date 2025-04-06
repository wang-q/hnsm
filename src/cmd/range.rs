use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("range")
        .about("Extract sequence regions by coordinates")
        .after_help(
            r###"
This command extracts sequence regions from FA files using genomic coordinates.

Range format:
    seq_name(strand):start-end

* seq_name: Required, sequence identifier
* strand: Optional, + (default) or -
* start-end: Required, 1-based coordinates

Examples:
    Mito
    I:1-100
    I(+):90-150
    S288c.I(-):190-200
    II:21294-22075
    II:23537-24097

Input methods:
* Command line: hnsm range input.fa "chr1:1-1000"
* Range file: hnsm range input.fa -r ranges.txt

Features:
* Supports BGZF compressed files (.gz)
* Automatic index creation (.loc)
* LRU caching for better performance
* Reverse complement for negative strand

Notes:
* Cannot read from stdin or gzip
* All coordinates (<start> and <end>) are based on the positive strand, regardless of the specified strand.
* Sort range file for better performance
* Cache size affects memory usage

Examples:
1. Single range:
   hnsm range input.fa "chr1:1-1000"

2. Multiple ranges:
   hnsm range input.fa "chr1:1-1000" "chr2(-):2000-3000"

3. From range file with larger cache:
   hnsm range input.fa -r ranges.txt -c 10

4. Force update the index file:
   hnsm range input.fa "chr1:1-1000" --update

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
            Arg::new("rgfile")
                .long("rgfile")
                .short('r')
                .num_args(1)
                .help("File of regions, one per line"),
        )
        .arg(
            Arg::new("cache")
                .long("cache")
                .short('c')
                .value_parser(value_parser!(std::num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Set the capacity of the LRU cache"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
        .arg(
            Arg::new("update")
                .long("update")
                .short('u')
                .action(ArgAction::SetTrue)
                .help("Force update the .loc index file"),
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
        noodles_fasta::io::writer::Builder::default()
            .set_line_base_count(usize::MAX)
            .build_from_writer(writer)
    };

    let is_bgzf = {
        let path = std::path::Path::new(infile);
        path.extension() == Some(std::ffi::OsStr::new("gz"))
    };

    let mut ranges = if args.contains_id("ranges") {
        args.get_many::<String>("ranges")
            .unwrap()
            .cloned()
            .collect()
    } else {
        vec![]
    };

    if args.contains_id("rgfile") {
        let mut rgs = intspan::read_first_column(args.get_one::<String>("rgfile").unwrap());
        ranges.append(&mut rgs);
    }

    let opt_cache = *args.get_one::<std::num::NonZeroUsize>("cache").unwrap();
    let mut cache: lru::LruCache<String, noodles_fasta::Record> = lru::LruCache::new(opt_cache);

    //----------------------------
    // Open files
    //----------------------------
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file() || args.get_flag("update") {
        hnsm::create_loc(infile, &loc_file, is_bgzf)?;
    }
    let loc_of: indexmap::IndexMap<String, (u64, usize)> = hnsm::load_loc(&loc_file)?;

    let mut reader = if is_bgzf {
        hnsm::Input::Bgzf(noodles_bgzf::indexed_reader::Builder::default().build_from_path(infile)?)
    } else {
        hnsm::Input::File(std::fs::File::open(std::path::Path::new(infile))?)
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

        if !cache.contains(&seq_id) {
            let record = hnsm::record_rg(&mut reader, &loc_of, &seq_id)?;
            cache.put(seq_id.clone(), record);
        }

        let record: &noodles_fasta::Record = cache.get(&seq_id).unwrap();

        // name only
        if *rg.start() == 0 {
            fa_out.write_record(record)?;
            continue;
        }

        let definition = noodles_fasta::record::Definition::new(rg.to_string(), None);

        // slice here is 1-based
        let start = noodles_core::Position::new(*rg.start() as usize).unwrap();
        let end = noodles_core::Position::new(*rg.end() as usize).unwrap();

        let mut slice = record.sequence().slice(start..=end).unwrap();
        if rg.strand() == "-" {
            slice = slice.complement().rev().collect::<Result<_, _>>()?;
        }
        let record_rg = noodles_fasta::Record::new(definition, slice);

        fa_out.write_record(&record_rg)?;
    }

    Ok(())
}

// fn print_type_of<T: ?Sized>(_: &T) {
//     println!("{}", std::any::type_name::<T>())
// }
