use clap::*;
use hnsm::Minimizer;
use noodles_fasta as fasta;
use std::collections::{BTreeMap, BTreeSet, HashSet};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("sketch")
        .about("Extract one FA record")
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
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
    let reader = intspan::reader(args.get_one::<String>("infile").unwrap());
    let mut fa_in = fasta::io::Reader::new(reader);

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let mut fac = hnsm::JumpingMinimizer {
        w: 7,
        k: 6,
        hasher: hnsm::FxHash,
    };
    let mut set_of = BTreeMap::new();
    let mut names = vec![];

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into()).unwrap();
        let seq = record.sequence();

        let minimizers = fac.minimizer(&seq[..]);
        let mut set = HashSet::new();
        for (_, hash) in &minimizers {
            set.insert(*hash);
        }
        names.push(name.clone());
        set_of.insert(name, set);
    }
    // eprintln!("set_of = {:#?}", set_of);

    for i in &names {
        for j in &names {
            let set1 = set_of.get(i).unwrap();
            let set2 = set_of.get(j).unwrap();
            let inter: HashSet<_> = set1.intersection(&set2).collect();
            let union: HashSet<_> = set1.union(&set2).collect();

            let dist = 1.0 - ((inter.len() as f64) / (union.len() as f64));

            writer.write_fmt(format_args!("{}\t{}\t{}\n", i, j, dist))?;
        }
    }

    Ok(())
}
