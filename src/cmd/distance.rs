use clap::*;
use hnsm::Minimizer;
use noodles_fasta as fasta;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("distance")
        .about("Estimate distances between DNA/protein sequences")
        .after_help(
            r###"
* <infile> can be plain text or bgzf but not stdin or gzip
* The outputs:
    n1, n2, mash, jaccard, containment

* Minimizers
    Given a $(k + w - 1)$-mer, consider the $w$ contained $k$-mers. The (rightmost) $k$-mer with
    minimal hash (for some given hash function) is the minimizer.

* We use minimizers here to sample kmers
    * For proteins, the length is short, so the window size can be set as: `-k 7 -w 1`
    * DNA: `-k 21 -w 5`
    * Increasing the window size speeds up processing

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("hasher")
                .long("hasher")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("fx"),
                    builder::PossibleValue::new("murmur"),
                ])
                .default_value("fx")
                .help("Set the hash algorithm"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("Kmer size"),
        )
        .arg(
            Arg::new("window")
                .long("window")
                .short('w')
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Window size"),
        )
        .arg(
            Arg::new("sim")
                .long("sim")
                .action(ArgAction::SetTrue)
                .help("Convert distance to similarity"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("Number of threads"),
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
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = fasta::io::Reader::new(reader);

    let opt_hasher = args.get_one::<String>("hasher").unwrap();
    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_window = *args.get_one::<usize>("window").unwrap();
    let is_sim = args.get_flag("sim");

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    let mut set_of = HashMap::new();
    let mut names = vec![];

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into()).unwrap();
        let seq = record.sequence();

        let minimizers = match opt_hasher.as_str() {
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

        let set: HashSet<u64> = HashSet::from_iter(minimizers.iter().map(|t| t.1));
        names.push(name.clone());
        set_of.insert(name, set);
    }
    // eprintln!("set_of = {:#?}", set_of);

    for n1 in &names {
        for n2 in &names {
            let s1 = set_of.get(n1).unwrap();
            let s2 = set_of.get(n2).unwrap();
            let inter: HashSet<_> = s1.intersection(&s2).collect();
            let union: HashSet<_> = s1.union(&s2).collect();

            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let containment = (inter.len() as f64) / (s1.len() as f64);
            // https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation
            let mash: f64 = if jaccard == 0.0 {
                1.0
            } else {
                ((-1.0 / 7.0f64) * ((2.0 * jaccard) / (1.0f64 + jaccard)).ln()).abs()
            };

            writer.write_fmt(format_args!(
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\n",
                n1,
                n2,
                if is_sim { 1.0 - mash } else { mash },
                jaccard,
                containment
            ))?;
        }
    }

    Ok(())
}
