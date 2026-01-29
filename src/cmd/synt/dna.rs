use clap::*;
use hnsm::libs::synteny::algo::SyntenyFinder;
use std::collections::HashMap;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("dna")
        .about("Synteny detection using minimizer graphs")
        .arg(
            Arg::new("infiles")
                .help("Input FASTA files")
                .required(true)
                .num_args(1..)
                .index(1),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .long("kmer")
                .help("K-mer size")
                .default_value("24")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("divergence")
                .short('d')
                .long("divergence")
                .help("Approximate sequence divergence (%)")
                .value_parser(value_parser!(f64)),
        )
        .arg(
            Arg::new("rounds")
                .short('r')
                .long("rounds")
                .help("Window sizes for iterative refinement (comma-separated). Defaults depend on divergence.")
                .value_parser(value_parser!(String)),
        )
        .arg(
            Arg::new("block_size")
                .short('b')
                .long("block-size")
                .help("Minimum synteny block size (bp). Defaults depend on divergence.")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("chain_gap")
                .long("chain-gap")
                .help("Maximum gap size between chained minimizers (bp). Defaults depend on divergence.")
                .value_parser(value_parser!(u32)),
        )
        .arg(
            Arg::new("min_weight")
                .long("min-weight")
                .help("Minimum edge weight (number of supporting genomes)")
                .default_value("2")
                .value_parser(value_parser!(usize)),
        )
        .arg(
            Arg::new("max_freq")
                .long("max-freq")
                .help("Maximum k-mer frequency to consider (filter repeats)")
                .default_value("1000")
                .value_parser(value_parser!(u32)),
        )
        .arg(
            Arg::new("soft_mask")
                .long("soft-mask")
                .action(clap::ArgAction::SetTrue)
                .help("Ignore soft-masked repeats (lowercase bases)"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .help("Output filename. [stdout] for screen")
                .default_value("stdout"),
        )
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .short('v')
                .action(clap::ArgAction::SetTrue)
                .help("Verbose output"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = matches
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let k = *matches.get_one::<usize>("kmer").unwrap();
    let min_weight = *matches.get_one::<usize>("min_weight").unwrap();
    let max_freq = *matches.get_one::<u32>("max_freq").unwrap();
    let soft_mask = matches.get_flag("soft_mask");
    let outfile = matches.get_one::<String>("outfile").unwrap();
    let verbose = matches.get_flag("verbose");

    // Default parameters based on divergence
    let divergence = matches.get_one::<f64>("divergence");

    // Logic from ntSynt
    let (default_rounds, default_block_size, default_chain_gap) =
        if let Some(d) = divergence {
            if *d < 1.0 {
                ("100,10", 500, 10000)
            } else if *d <= 10.0 {
                ("250,100", 1000, 50000)
            } else {
                ("500,250", 10000, 50000)
            }
        } else {
            // Fallback if no divergence specified, use "medium" defaults or what was previously hardcoded
            ("1000,100,10", 0, 20000)
        };

    let rounds_str = matches
        .get_one::<String>("rounds")
        .map(|s| s.as_str())
        .unwrap_or(default_rounds);
    let block_size = *matches
        .get_one::<usize>("block_size")
        .unwrap_or(&default_block_size);
    let chain_gap = *matches
        .get_one::<u32>("chain_gap")
        .unwrap_or(&default_chain_gap);

    let rounds: Vec<usize> = rounds_str
        .split(',')
        .map(|s| s.trim().parse::<usize>().expect("Invalid window size"))
        .collect();

    if verbose {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Info)
            .init();
    } else {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Warn)
            .init();
    }

    let finder = SyntenyFinder::new(
        k,
        rounds,
        min_weight,
        max_freq,
        block_size,
        chain_gap,
        soft_mask,
    );

    // Pre-scan to build seq_names map
    // We could do this inside provider but we need names for output.
    // Let's do a quick pass or just build it during the first run?
    // But SyntenyFinder runs multiple rounds.
    // It's safer to build it once.
    let mut seq_names: HashMap<u32, String> = HashMap::new();
    let mut global_seq_id = 0;

    // We can just iterate files once to get names, or assume names are stable.
    // Since we rely on global_seq_id being consistent, we must process files in same order.
    for infile in &infiles {
        let reader = intspan::reader(infile);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);

        let path = std::path::Path::new(infile);
        let file_name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        let file_stem = file_name.split('.').next().unwrap_or("unknown");
        let use_prefix = infiles.len() > 1;

        for result in fa_in.records() {
            let record = result?;
            global_seq_id += 1;
            let name = String::from_utf8(record.name().into())?;

            let name = if use_prefix {
                format!("{}.{}", file_stem, name)
            } else {
                name
            };

            seq_names.insert(global_seq_id, name);
        }
    }

    let mut writer = intspan::writer(outfile);
    writeln!(
        writer,
        "# Block_ID\tRange\tCount\tRound"
    )?;
    let mut block_counter = 0;

    let provider = |emit: &mut dyn FnMut(&str, &[u8])| -> anyhow::Result<()> {
        for infile in &infiles {
            let reader = intspan::reader(infile);
            let mut fa_in = noodles_fasta::io::Reader::new(reader);
            for result in fa_in.records() {
                let record = result?;
                let name = String::from_utf8(record.name().into())?;
                let seq = record.sequence();
                emit(&name, seq.as_ref());
            }
        }
        Ok(())
    };

    finder.run(provider, |w, block| {
        let mut ranges: Vec<_> = block.ranges.values().collect();
        ranges.sort_by_key(|r| r.seq_id);

        let flip = if let Some(first) = ranges.first() {
            !first.strand
        } else {
            false
        };

        for range in ranges {
            let seq_name = seq_names
                .get(&range.seq_id)
                .cloned()
                .unwrap_or_else(|| format!("Seq_{}", range.seq_id));
            
            let current_strand = if flip { !range.strand } else { range.strand };
            let strand_char = if current_strand { '+' } else { '-' };
            
            let _ = writeln!(
                writer,
                "{}\t{}({}):{}-{}\t{}\t{}",
                block_counter, seq_name, strand_char, range.start, range.end, range.count, w
            );
        }
        block_counter += 1;
    })?;

    Ok(())
}
