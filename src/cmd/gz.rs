use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("gz")
        .about("Compressing a file using the BGZF format")
        .after_help(
            r###"
This command compresses a file using BGZF (Blocked Gzip Format).

Features:
* Parallel compression with multiple threads
* Creates index file (.gzi) for random access
* Supports stdin as input
* Preserves original file

Output files:
* <infile>.gz: Compressed file
* <infile>.gz.gzi: Index file

Notes:
* Cannot compress already gzipped files
* Requires bgzip in PATH for indexing
* Default thread count is 1
* Index creation is automatic

Examples:
1. Compress a file with default settings, and the outfile is input.fa.gz:
   hnsm gz input.fa

2. Multi-threaded compression:
   hnsm gz input.fa -p 4

3. From stdin with custom output:
   cat input.fa | hnsm gz stdin -o output.fa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input FA file to compress"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(std::num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Number of threads for parallel compression"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .help("Output filename (default: <infile>.gz)"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    let opt_parallel = *args.get_one::<std::num::NonZeroUsize>("parallel").unwrap();

    let outfile = if args.contains_id("outfile") {
        format!("{}.gz", args.get_one::<String>("outfile").unwrap())
    } else {
        format!("{}.gz", infile)
    };

    //----------------------------
    // Input
    //----------------------------
    let mut reader: Box<dyn std::io::BufRead> = if infile == "stdin" {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        let path = std::path::Path::new(infile);
        let file = match std::fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        Box::new(std::io::BufReader::new(file))
    };

    let mut writer = noodles_bgzf::MultithreadedWriter::with_worker_count(
        opt_parallel,
        Box::new(std::io::BufWriter::new(
            std::fs::File::create(&outfile).unwrap(),
        )),
    );

    //----------------------------
    // Output
    //----------------------------
    std::io::copy(&mut reader, &mut writer)?;
    writer.finish()?;

    let bin = if let Ok(pth) = which::which("bgzip") {
        pth.to_string_lossy().to_string()
    } else {
        "".to_string()
    };

    if bin.is_empty() {
        return Err(anyhow::anyhow!(
            "Can't find `bgzip` in $PATH. .gzi not created"
        ));
    }

    let res = std::process::Command::new(bin)
        .arg("-r")
        .arg(&outfile)
        .output()?;
    if !res.status.success() {
        return Err(anyhow::anyhow!("Command executed with failing error code"));
    }

    Ok(())
}
