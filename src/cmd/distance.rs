use clap::*;
use noodles_fasta as fasta;
use rayon::prelude::*;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("distance")
        .about("Estimate sequence distances using minimizers")
        .after_help(
            r###"
This command calculates pairwise distances between sequences in FA file(s) using minimizers.

* The outputs are printed to stdout in the following format:
    <sequence1> <sequence2> <mash_distance> <jaccard_index> <containment_index>
* With --merge
    <file1> <file2> <total1> <total2> <inter> <union> <mash_distance> <jaccard_index> <containment_index>

* Minimizers
    Given a $(k + w - 1)$-mer, consider the $w$ contained $k$-mers. The (rightmost) $k$-mer with
    minimal hash (for some given hash function) is the minimizer.

* We use minimizers here to sample kmers
    * For proteins, the length is short, so the window size can be set as: `-k 7 -w 2`
    * DNA: `-k 21 -w 5`
    * Increasing the window size speeds up processing

* Hash Algorithms (--hasher):
    * The `--hasher` parameter selects the hash algorithm used for minimizer calculation.
    * Available options:
        - `rapid`: RapidHash (default)
        - `fx`: FxHash
        - `murmur`: MurmurHash3
    * Note: The `mod` option is not a hash algorithm but a special mode for DNA sequences.

* Mod-Minimizer (--hasher mod):
    * It generates canonical k-mers, meaning that a sequence and its reverse complement
      are generating the same k-mer set.

* To get accurate pairwise sequence identities, use clustalo
  https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

* Input Modes:
    * By default (--list is false):
        * Single file: Treat the file as a sequence file and calculate pairwise distances
          for all sequences within it.
        * Two files: Treat both files as sequence files and calculate pairwise distances
          between sequences from the two files.
    * When --list is set:
        * Single file: Treat the file as a list file (each line is a path to a sequence file)
          and calculate pairwise distances for all sequences in the listed files.
        * Two files: Treat both files as list files and calculate pairwise distances
          between sequences from the two list files.

* --merge Behavior:
  - By default (--merge is false):
    * Distances are calculated between individual sequences.
  - When --merge is set:
    * For a single sequence file: Merge all sequences within the file into a single set
      of minimizers. Note that comparing this set to itself (self-comparison) is not
      meaningful, as the distance will always be 0 and the similarity will always be 1.
    * For two sequence files: Merge all sequences within each file into a single set,
      and calculate distances between the two sets.
    * When --list is set, --merge operates on each sequence file individually:
      - For each file listed in the list file, merge all sequences within that file
        into a single set, and calculate distances between these sets.
      - The merging does not span across multiple files listed in the list file.

Examples:
1. Calculate distances with default parameters:
   hnsm distance input.fa

2. Use Mod-Minimizer for DNA sequences (canonical k-mers):
   hnsm distance input.fa --hasher mod -k 21 -w 5

3. Compare two FA files:
   hnsm distance file1.fa file2.fa

4. Merge all sequences in a file and compare to another:
   hnsm distance file1.fa file2.fa --merge

5. Treat input as a list file and calculate distances:
   hnsm distance list.txt --list

6. Use 4 threads for parallel processing:
   hnsm distance input.fa --parallel 4

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
            Arg::new("sim")
                .long("sim")
                .action(ArgAction::SetTrue)
                .help("Convert distance to similarity (1 - distance)"),
        )
        .arg(
            Arg::new("zero")
                .long("zero")
                .action(ArgAction::SetTrue)
                .help("Also write results with zero Jaccard index"),
        )
        .arg(
            Arg::new("merge")
                .long("merge")
                .action(ArgAction::SetTrue)
                .help("Merge all sequences within a file into a single set for comparison"),
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
struct MinimizerEntry {
    name: String,
    set: rapidhash::RapidHashSet<u64>,
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();

    let is_sim = args.get_flag("sim");
    let is_zero = args.get_flag("zero");
    let is_merge = args.get_flag("merge"); // Whether to merge all sequences within a file
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
        let entries = load_entries(&paths, opt_hasher, opt_kmer, opt_window, is_merge)?;
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
        let entries1 = load_entries(&paths1, opt_hasher, opt_kmer, opt_window, is_merge)?;
        let entries2 = load_entries(&paths2, opt_hasher, opt_kmer, opt_window, is_merge)?;
        (entries1, entries2) // Calculate pairwise distances between the two sets
    };

    // Use rayon to parallelize the outer loop
    entries1.par_iter().for_each(|e1| {
        let mut lines = String::with_capacity(1024);
        for (i, e2) in entries2.iter().enumerate() {
            let (total1, total2, inter, union, mash, jaccard, containment) =
                calc_distances(&e1.set, &e2.set, opt_kmer);

            if !is_zero && jaccard == 0. {
                continue;
            }

            let out_string = if is_merge {
                format!(
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
                )
            } else {
                format!(
                    "{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                    e1.name,
                    e2.name,
                    if is_sim { 1.0 - mash } else { mash },
                    jaccard,
                    containment
                )
            };

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
    is_merge: bool,
) -> anyhow::Result<Vec<MinimizerEntry>> {
    let mut entries = Vec::new();

    for path in paths {
        let mut loaded = load_file(path, opt_hasher, opt_kmer, opt_window, is_merge)?;
        entries.append(&mut loaded);
    }

    Ok(entries)
}

fn load_file(
    infile: &str,
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
    is_merge: bool,
) -> anyhow::Result<Vec<MinimizerEntry>> {
    let reader = intspan::reader(infile);
    let mut fa_in = fasta::io::Reader::new(reader);

    let mut entries = vec![];
    // Set to merge all minimizers if --merge is true
    let mut all_set = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        let set: rapidhash::RapidHashSet<u64> =
            hnsm::seq_mins(&seq[..], opt_hasher, opt_kmer, opt_window)?;

        if is_merge {
            all_set.extend(set);
        } else {
            let entry = MinimizerEntry { name, set };
            entries.push(entry);
        }
    }

    if is_merge {
        let entry = MinimizerEntry {
            name: infile.to_string(),
            set: all_set,
        };
        entries.push(entry);
    }

    Ok(entries)
}

// Calculate Jaccard, Containment, and Mash distance between two sets
fn calc_distances(
    s1: &rapidhash::RapidHashSet<u64>,
    s2: &rapidhash::RapidHashSet<u64>,
    opt_kmer: usize,
) -> (usize, usize, usize, usize, f64, f64, f64) {
    let total1 = s1.len();
    let total2 = s2.len();

    let inter = s1.intersection(s2).cloned().count();
    let union = total1 + total2 - inter;

    let jaccard = inter as f64 / union as f64;
    let containment = inter as f64 / total1 as f64;
    // https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation
    let mash = if jaccard == 0.0 {
        1.0
    } else {
        ((-1.0 / opt_kmer as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()).abs()
    };

    (total1, total2, inter, union, mash, jaccard, containment)
}
