use clap::*;
use hnsm::Minimizer;
use noodles_fasta as fasta;
use rayon::prelude::*;
use std::iter::FromIterator;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("distance")
        .about("Estimate distances between DNA/protein sequences using minimizers")
        .after_help(
            r###"
This command calculates pairwise distances between sequences in FA file(s) using minimizers.

* The outputs are printed to stdout in the following format:
    n1  n2  mash    jaccard containment

* Minimizers
    Given a $(k + w - 1)$-mer, consider the $w$ contained $k$-mers. The (rightmost) $k$-mer with
    minimal hash (for some given hash function) is the minimizer.

* We use minimizers here to sample kmers
    * For proteins, the length is short, so the window size can be set as: `-k 7 -w 1`
    * DNA: `-k 21 -w 5`
    * Increasing the window size speeds up processing

* To get accurate pairwise sequence identities, use clustalo
  https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

Examples:
    1. Calculate distances with default parameters:
       hnsm distance input.fa

    2. Use MurmurHash instead of RapidHash:
       hnsm distance input.fa --hasher murmur

    3. Set custom k-mer size and window size:
       hnsm distance input.fa -k 21 -w 5

    4. Convert distance to similarity:
       hnsm distance input.fa --sim

    5. Use 4 threads for parallel processing:
       hnsm distance input.fa --parallel 4

    6. Compare two FA files:
       hnsm distance file1.fa file2.fa

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
            Arg::new("hasher")
                .long("hasher")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("rapid"),
                    builder::PossibleValue::new("fx"),
                    builder::PossibleValue::new("murmur"),
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
                .help("Convert distance to similarity"),
        )
        .arg(
            Arg::new("zero")
                .long("zero")
                .action(ArgAction::SetTrue)
                .help("Also write results with zero jaccard index"),
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
}

#[derive(Default, Clone)]
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
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let infiles = args
        .get_many::<String>("infiles")
        .unwrap()
        .map(|s| s.as_str())
        .collect::<Vec<_>>();

    //----------------------------
    // Ops
    //----------------------------
    let entries = load_file(infiles.get(0).unwrap(), opt_hasher, opt_kmer, opt_window);
    let others = if infiles.len() == 2 {
        load_file(infiles.get(1).unwrap(), opt_hasher, opt_kmer, opt_window)
    } else {
        entries.clone()
    };

    // Set the number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt_parallel)
        .build_global()?;

    // Use rayon to parallelize the outer loop
    entries.par_iter().for_each(|e1| {
        for e2 in &others {
            let (jaccard, containment, mash) = calc_distances(&e1.set, &e2.set);
            if !is_zero && jaccard == 0. {
                continue;
            }
            let out_string = format!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                e1.name,
                e2.name,
                if is_sim { 1.0 - mash } else { mash },
                jaccard,
                containment
            );
            print!("{}", out_string);
        }
    });

    Ok(())
}

fn load_file(
    infile: &str,
    opt_hasher: &String,
    opt_kmer: usize,
    opt_window: usize,
) -> Vec<MinimizerEntry> {
    let reader = intspan::reader(infile);
    let mut fa_in = fasta::io::Reader::new(reader);

    let mut entries = vec![];

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result.unwrap();

        let name = String::from_utf8(record.name().into()).unwrap();
        let seq = record.sequence();

        let minimizers = match opt_hasher.as_str() {
            "rapid" => hnsm::JumpingMinimizer {
                w: opt_window,
                k: opt_kmer,
                hasher: hnsm::RapidHash,
            }
            .minimizer(&seq[..]),
            "fx" => hnsm::JumpingMinimizer {
                w: opt_window,
                k: opt_kmer,
                hasher: hnsm::FxHash,
            }
            .minimizer(&seq[..]),
            "murmur" => hnsm::JumpingMinimizer {
                w: opt_window,
                k: opt_kmer,
                hasher: hnsm::MurmurHash3,
            }
            .minimizer(&seq[..]),
            _ => unreachable!(),
        };

        let set: rapidhash::RapidHashSet<u64> =
            rapidhash::RapidHashSet::from_iter(minimizers.iter().map(|t| t.1));
        let entry = MinimizerEntry { name, set };
        entries.push(entry);
    }

    entries
}

// Calculate Jaccard, Containment, and Mash distance between two sets
fn calc_distances(
    s1: &rapidhash::RapidHashSet<u64>,
    s2: &rapidhash::RapidHashSet<u64>,
) -> (f64, f64, f64) {
    let inter = s1.intersection(&s2).cloned().count();
    let union = s1.len() + s2.len() - inter;

    let jaccard = inter as f64 / union as f64;
    let containment = inter as f64 / s1.len() as f64;
    // https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation
    let mash = if jaccard == 0.0 {
        1.0
    } else {
        ((-1.0 / 7.0) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()).abs()
    };

    (jaccard, containment, mash)
}
