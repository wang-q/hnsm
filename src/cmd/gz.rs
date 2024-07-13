use clap::*;

use std::io::{self};
use std::{fs, num, path, process};

use noodles_bgzf as bgzf;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("gz")
        .about("Compressing a file using the blocked gzip format (BGZF)")
        .after_help(
            r###"
* <infile> can be stdin but not gzipped
* The upstream library `noodles` not yet implemented write .gzi
* `bgzip -r` will be called to create .gzi
* The output is hardcoded as "outfile.gz" and "outfile.gz.gzi", for more flexible output use `bgzip`

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(num::NonZeroUsize))
                .num_args(1)
                .default_value("1")
                .help("Running in parallel mode, the number of threads"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .help("Output filename. Default is infile.gz"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    let opt_parallel = *args.get_one::<num::NonZeroUsize>("parallel").unwrap();

    let outfile = format!(
        "{}.gz",
        if args.contains_id("outfile") {
            args.get_one::<String>("outfile").unwrap()
        } else {
            infile
        }
    );

    //----------------------------
    // Open files
    //----------------------------
    let mut reader: Box<dyn io::BufRead> = if infile == "stdin" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        let path = path::Path::new(infile);
        let file = match fs::File::open(path) {
            Err(why) => panic!("could not open {}: {}", path.display(), why),
            Ok(file) => file,
        };

        Box::new(io::BufReader::new(file))
    };

    let mut writer = bgzf::MultithreadedWriter::with_worker_count(
        opt_parallel,
        Box::new(io::BufWriter::new(fs::File::create(&outfile).unwrap())),
    );

    //----------------------------
    // Output
    //----------------------------
    io::copy(&mut reader, &mut writer)?;
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

    let res = process::Command::new(bin)
        .arg("-r")
        .arg(&outfile)
        .output()?;
    if !res.status.success() {
        return Err(anyhow::anyhow!("Command executed with failing error code"));
    }

    Ok(())
}
