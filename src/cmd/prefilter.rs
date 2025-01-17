use clap::*;
use cmd_lib::*;
use noodles_bgzf as bgzf;
use noodles_fasta as fasta;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("prefilter")
        .about("Prefilter genome/metagenome assembly by hypervectors")
        .after_help(
            r###"
* <infile> can be plain text or bgzf but not stdin or gzip

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("match")
                .required(true)
                .index(2)
                .help("The match file"),
        )
        .arg(
            Arg::new("chunk")
                .long("chunk")
                .short('c')
                .num_args(1)
                .default_value("100000")
                .value_parser(value_parser!(usize))
                .help("Size of each chunk in bytes"),
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
    let match_file = args.get_one::<String>("match").unwrap();

    let opt_chunk = *args.get_one::<usize>("chunk").unwrap();

    let is_bgzf = {
        let path = std::path::Path::new(infile);
        path.extension() == Some(std::ffi::OsStr::new("gz"))
    };

    //----------------------------
    // Open files
    //----------------------------
    let loc_file = format!("{}.loc", infile);
    if !std::path::Path::new(&loc_file).is_file() {
        hnsm::create_loc(infile, &loc_file, is_bgzf)?;
    }

    let mut reader = if is_bgzf {
        hnsm::Input::Bgzf(bgzf::indexed_reader::Builder::default().build_from_path(infile)?)
    } else {
        hnsm::Input::File(std::fs::File::open(std::path::Path::new(infile))?)
    };

    // Split .loc file into chunks
    let chunks = split_loc_file(&loc_file, opt_chunk)?;

    let hnsm = std::env::current_exe()?.display().to_string();

    for (first, offset, size) in chunks.iter() {
        let chunk = hnsm::read_offset(&mut reader, *offset, *size)?;

        let mut temp_file = tempfile::NamedTempFile::new()?;
        temp_file.write_all(&chunk)?;
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        run_cmd!(
            ${hnsm} sixframe ${temp_path} --len 15 |
                ${hnsm} distance stdin ${match_file}
        )?;

        // eprintln!(
        //     "Processed chunk {}: first sequence = {}, temp file = {:?}",
        //     i, first, temp_path
        // );
        // hnsm::pause();
    }

    Ok(())
}

// Split .loc file into chunks
fn split_loc_file(loc_file: &str, chunk_size: usize) -> anyhow::Result<Vec<(String, u64, usize)>> {
    // Load .loc file
    let loc_of: indexmap::IndexMap<String, (u64, usize)> = hnsm::load_loc(loc_file)?;

    let mut chunks: Vec<(String, u64, usize)> = Vec::new();
    let mut cur_size = 0;
    let mut cur_start_offset = 0;
    let mut cur_first_seq = String::new();

    // Iterate over each sequence in the .loc file
    for (seq_id, &(offset, size)) in &loc_of {
        // If the current chunk size exceeds the specified size,
        //   record the current chunk and start a new one
        if cur_size + size > chunk_size && !cur_first_seq.is_empty() {
            chunks.push((cur_first_seq.clone(), cur_start_offset, cur_size));
            cur_size = 0;
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        // If it's the first sequence, record the start offset and sequence name
        if cur_size == 0 {
            cur_start_offset = offset;
            cur_first_seq = seq_id.clone();
        }

        // Update the current chunk size
        cur_size += size;
    }

    // Add the last chunk
    if !cur_first_seq.is_empty() {
        chunks.push((cur_first_seq, cur_start_offset, cur_size));
    }

    Ok(chunks)
}
