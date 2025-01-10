use clap::*;
use rayon::prelude::*;
use std::io::BufRead;
use std::simd::prelude::*;

const LANES: usize = 8; // 32 * 8 = 256, AVX2

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("similarity")
        .about("Calculate similarity between vectors")
        .after_help(
            r###"
This command calculates pairwise similarity between vectors in input file(s).

modes:
    * euclidean distance
        * --mode euclid
    * euclidean distance to similarity
        * --mode euclid --sim
    * binary euclidean distance
        * --mode euclid --bin
    * binary euclidean distance to dissimilarity
        * --mode euclid --bin --sim --dis

    * cosine similarity, -1 -- 1
        * --mode cosine
    * cosine distance, 0 -- 2
        * --mode cosine --dis
    * binary cosine similarity
        * --mode cosine --bin
    * binary cosine similarity
        * --mode cosine --bin --dis

    * jaccard index
        * --mode jaccard --bin
    * weighted jaccard similarity
        * --mode jaccard

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..=2)
                .index(1)
                .required(true)
                .help("Input filenames. [stdin] for standard input"),
        )
        .arg(
            Arg::new("mode")
                .long("mode")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("euclid"),
                    builder::PossibleValue::new("cosine"),
                    builder::PossibleValue::new("jaccard"),
                ])
                .default_value("euclid")
                .help("Mode of calculation"),
        )
        .arg(
            Arg::new("bin")
                .long("bin")
                .action(ArgAction::SetTrue)
                .help("Treat values in list as binary (0 or 1)"),
        )
        .arg(
            Arg::new("sim")
                .long("sim")
                .action(ArgAction::SetTrue)
                .help("Convert distance to similarity"),
        )
        .arg(
            Arg::new("dis")
                .long("dis")
                .action(ArgAction::SetTrue)
                .help("Convert to dissimilarity"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads for parallel processing"),
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
    let opt_mode = args.get_one::<String>("mode").unwrap();

    let is_bin = args.get_flag("bin");
    let is_sim = args.get_flag("sim");
    let is_dis = args.get_flag("dis");

    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>();

    // Create a channel for sending results to the writer thread
    let (sender, receiver) = crossbeam::channel::bounded::<String>(256);

    // Spawn a writer thread
    let output = args.get_one::<String>("outfile").unwrap().to_string();
    let writer_thread = std::thread::spawn(move || {
        let mut writer = intspan::writer(&output);
        for result in receiver {
            writer.write_all(result.as_bytes()).unwrap();
        }
    });

    // Set the number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt_parallel)
        .build_global()?;

    //----------------------------
    // Ops
    //----------------------------
    let entries = load_file(infiles.get(0).unwrap(), is_bin);
    let others = if infiles.len() == 2 {
        load_file(infiles.get(1).unwrap(), is_bin)
    } else {
        entries.clone()
    };

    // Use rayon to parallelize the outer loop
    entries.par_iter().for_each(|e1| {
        let mut lines = "".to_string();
        for (i, e2) in others.iter().enumerate() {
            let score = calc(e1.list(), e2.list(), opt_mode, is_sim, is_dis);
            let out_string = format!("{}\t{}\t{:.4}\n", e1.name(), e2.name(), score);

            lines.push_str(&out_string);
            if i > 1 && i % 1000 == 0 {
                sender.send(lines.clone()).unwrap();
                lines.clear();
            }
        }
        if !lines.is_empty() {
            sender.send(lines).unwrap();
        }
    });

    // Drop the sender to signal the writer thread to exit
    drop(sender);
    // Wait for the writer thread to finish
    writer_thread.join().unwrap();

    Ok(())
}

fn load_file(infile: &str, is_bin: bool) -> Vec<hnsm::AsmEntry> {
    let mut entries = vec![];
    let reader = intspan::reader(infile);
    'LINE: for line in reader.lines().map_while(Result::ok) {
        let mut entry = hnsm::AsmEntry::parse(&line);
        if entry.name().is_empty() {
            continue 'LINE;
        }
        if is_bin {
            let bin_list = entry
                .list()
                .iter()
                .map(|e| if *e > 0.0 { 1.0 } else { 0.0 })
                .collect::<Vec<f32>>();
            entry = hnsm::AsmEntry::from(entry.name(), &bin_list);
        }
        entries.push(entry);
    }
    entries
}

fn calc(l1: &[f32], l2: &[f32], mode: &str, is_sim: bool, is_dis: bool) -> f32 {
    let mut score = match mode {
        "euclid" => euclidean_distance(l1, l2),
        "cosine" => cosine_similarity(l1, l2),
        "jaccard" => weighted_jaccard_similarity(l1, l2),
        _ => unreachable!(),
    };

    if is_sim {
        score = d2s(score);
    }
    if is_dis {
        score = dis(score);
    }

    score
}

// https://www.maartengrootendorst.com/blog/distances/
// https://crates.io/crates/semanticsimilarity_rs
fn euclidean_distance(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        let diff = x - y;
        *d = diff * diff;
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        let diff = f32x8::from_array(*x) - f32x8::from_array(*y);
        sums += diff * diff;
    });

    sums.reduce_sum().sqrt()
}

fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = x * y;
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::from_array(*x) * f32x8::from_array(*y);
    });

    sums.reduce_sum()
}

fn norm(a: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();

    let mut sums = [0.0; LANES];
    for (x, d) in std::iter::zip(a_extra, &mut sums) {
        *d = x * x;
    }

    let mut sums = f32x8::from_array(sums);
    a_chunks.into_iter().for_each(|x| {
        sums += f32x8::from_array(*x) * f32x8::from_array(*x);
    });

    sums.reduce_sum().sqrt()
}

fn cosine_similarity(a: &[f32], b: &[f32]) -> f32 {
    let dot_product = dot_product(a, b);
    let denominator = norm(a) * norm(b);

    if denominator == 0.0 {
        0.0
    } else {
        dot_product / denominator
    }
}

fn jaccard_intersection(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = f32::min(*x, *y);
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::simd_min(f32x8::from_array(*x), f32x8::from_array(*y));
    });

    sums.reduce_sum()
}

fn jaccard_union(a: &[f32], b: &[f32]) -> f32 {
    let (a_extra, a_chunks): (&[f32], &[[f32; LANES]]) = a.as_rchunks();
    let (b_extra, b_chunks): (&[f32], &[[f32; LANES]]) = b.as_rchunks();

    let mut sums = [0.0; LANES];
    for ((x, y), d) in std::iter::zip(a_extra, b_extra).zip(&mut sums) {
        *d = f32::max(*x, *y);
    }

    let mut sums = f32x8::from_array(sums);
    std::iter::zip(a_chunks, b_chunks).for_each(|(x, y)| {
        sums += f32x8::simd_max(f32x8::from_array(*x), f32x8::from_array(*y));
    });

    sums.reduce_sum()
}

fn weighted_jaccard_similarity(a: &[f32], b: &[f32]) -> f32 {
    let numerator = jaccard_intersection(a, b);
    let denominator = jaccard_union(a, b);

    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

// Schölkopf, B. (2000). The kernel trick for distances. In Neural Information Processing Systems, pages 301-307.
// https://stats.stackexchange.com/questions/146309/turn-a-distance-measure-into-a-kernel-function
// https://stats.stackexchange.com/questions/158279/how-i-can-convert-distance-euclidean-to-similarity-score
fn d2s(dist: f32) -> f32 {
    1.0 / dist.abs().exp()
}

fn dis(dist: f32) -> f32 {
    1.0 - dist
}
