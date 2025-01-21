use clap::*;
use noodles_fasta as fasta;
use rayon::prelude::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("hv")
        .about("Estimate distances between DNA/protein files using hypervectors")
        .after_help(
            r###"
This command calculates pairwise distances between files in FA file(s) using minimizers and hypervectors.

* The outputs are printed to stdout in the following format:
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Minimizers and Hash Algorithms are the same as `hnsm distance`

* Input Modes:
    * For a single sequence file: Merge all sequences within the file into a single hypervector.
      Note that comparing this set to itself (self-comparison) is not meaningful,
      as the distance will always be 0 and the similarity will always be 1.
    * For two sequence files: Merge all sequences within each file into a single hypervector,
      and calculate distances between the two hypervectors.
    * When --list is set:
      - For each file listed in the list file, merge all sequences within that file
        into a single hypervector, and calculate distances between these hypervectors.
      - The merging does not span across multiple files listed in the list file.

Examples:
1. Merge all sequences in a file and compare to another:
   hnsm hv file1.fa file2.fa

2. Use Mod-Minimizer for DNA sequences (canonical k-mers):
   hnsm hv file1.fa file2.fa --hasher mod -k 21 -w 5

3. Treat input as a list file and calculate distances:
   hnsm hv list.txt --list

4. Use 4 threads for parallel processing:
   hnsm hv input.fa --parallel 4

5. Perform six-frame translation on a FA file and match to another
    hnsm sixframe input.fa |
        hnsm hv stdin match.fa

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..=2)
                .index(1)
                .required(true)
                .help("Input FA/list file(s). [stdin] for standard input"),
        )
        .arg(
            Arg::new("hasher")
                .long("hasher")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("rapid"),
                    builder::PossibleValue::new("fx"),
                    builder::PossibleValue::new("murmur"),
                    builder::PossibleValue::new("mod"),
                ])
                .default_value("rapid")
                .help("Hash algorithm to use"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("K-mer size"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .short('w')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Window size for minimizers"),
        )
        .arg(
            Arg::new("dim")
                .long("dim")
                .short('d')
                .num_args(1)
                .default_value("4096")
                .value_parser(value_parser!(usize))
                .help("The dimension size should be a multiple of 32."),
        )
        .arg(
            Arg::new("sim")
                .long("sim")
                .action(ArgAction::SetTrue)
                .help("Convert distance to similarity (1 - distance)"),
        )
        .arg(
            Arg::new("list")
                .long("list")
                .action(ArgAction::SetTrue)
                .help("Treat infiles as list files, where each line is a path to a sequence file"),
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

#[derive(Debug, Default, Clone)]
struct HvEntry {
    name: String,
    set: Vec<i32>,
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();
    let opt_dim = *args.get_one::<usize>("dim").unwrap();

    let is_sim = args.get_flag("sim");
    let is_list = args.get_flag("list"); // Whether to treat infiles as list files
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

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

    let infiles = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>();

    //----------------------------
    // Ops
    //----------------------------
    // Load data based on the number of input files and the --list flag
    let (entries1, entries2) = if infiles.len() == 1 {
        // Single file
        let paths = if is_list {
            intspan::read_first_column(infiles[0])
        } else {
            vec![infiles[0].to_string()] // Treat the input as a sequence file
        };
        let entries = load_entries(&paths, opt_hasher, opt_kmer, opt_window, opt_dim)?;
        (entries.clone(), entries) // Calculate pairwise distances within the same set
    } else {
        // Two files
        let paths1 = if is_list {
            intspan::read_first_column(infiles[0])
        } else {
            vec![infiles[0].to_string()]
        };
        let paths2 = if is_list {
            intspan::read_first_column(infiles[1])
        } else {
            vec![infiles[1].to_string()]
        };
        let entries1 = load_entries(&paths1, opt_hasher, opt_kmer, opt_window, opt_dim)?;
        let entries2 = load_entries(&paths2, opt_hasher, opt_kmer, opt_window, opt_dim)?;
        (entries1, entries2) // Calculate pairwise distances between the two sets
    };

    // Use rayon to parallelize the outer loop
    entries1.par_iter().for_each(|e1| {
        let mut lines = String::with_capacity(1024);
        for (i, e2) in entries2.iter().enumerate() {
            let (total1, total2, inter, union, mash, jaccard, containment) =
                calc_distances(&e1.set, &e2.set, opt_kmer);

            let out_string = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name,
                e2.name,
                total1,
                total2,
                inter,
                union,
                if is_sim { 1.0 - mash } else { mash },
                jaccard,
                containment
            );

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

// Load entries from a list of paths
fn load_entries(
    paths: &[String],
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
    opt_dim: usize,
) -> anyhow::Result<Vec<HvEntry>> {
    let mut entries = Vec::new();

    for path in paths {
        let mut loaded = load_file(path, opt_hasher, opt_kmer, opt_window, opt_dim)?;
        entries.append(&mut loaded);
    }

    Ok(entries)
}

fn load_file(
    infile: &str,
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
    opt_dim: usize,
) -> anyhow::Result<Vec<HvEntry>> {
    let reader = intspan::reader(infile);
    let mut fa_in = fasta::io::Reader::new(reader);

    let mut file_set = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;
        let seq = record.sequence();

        let set: rapidhash::RapidHashSet<u64> =
            hnsm::seq_mins(&seq[..], opt_hasher, opt_kmer, opt_window)?;

        file_set.extend(set);
    }

    let hv: Vec<i32> = hnsm::hash_hv(&file_set, opt_dim);
    let entry = HvEntry {
        name: infile.to_string(),
        set: hv,
    };

    Ok(vec![entry])
}

// Calculate Jaccard, Containment, and Mash distance between two sets
fn calc_distances(
    s1: &[i32],
    s2: &[i32],
    opt_kmer: usize,
) -> (usize, usize, usize, usize, f32, f32, f32) {
    let card1 = hnsm::hv_cardinality(s1);
    let card2 = hnsm::hv_cardinality(s2);

    let inter = hnsm::hv_dot(s1, s2).min(card1 as f32).min(card2 as f32);
    let union = card1 as f32 + card2 as f32 - inter;

    let jaccard = inter / union;
    let containment = inter / card1 as f32;
    // https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation
    let mash = if jaccard == 0.0 {
        1.0
    } else {
        ((-1.0 / opt_kmer as f32) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()).abs()
    };

    (
        card1,
        card2,
        inter as usize,
        union as usize,
        mash,
        jaccard,
        containment,
    )
}
