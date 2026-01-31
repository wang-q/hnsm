use clap::*;

use hnsm::libs::synteny::chain::{Anchor, Chain, ChainOpt, DagChainer};
use itertools::Itertools;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dag")
        .about("Synteny chains via DAGchainer algorithm")
        .after_help(
            r###"
Algorithm adopted from `DAGchainer`

# Algorithm & Scoring

1. **The Concept (2D Plane)**
   - Imagine Genome 1 on the X-axis and Genome 2 on the Y-axis.
   - Each gene match is a point (x, y) on this plane.
   - Synteny blocks appear as diagonal lines of points.

2. **The Goal**
   - Find the "best" path connecting these points (Chains).
   - "Best" means maximizing the Total Score.

3. **Scoring Logic**
   - **Reward**: High similarity matches increase the score.
     - Score = Similarity * Max_Score (default 50).
   - **Penalty**: Gaps between points decrease the score.
     - Gap Penalty = Gap_Open + (Distance * Gap_Extend).
     - Distance is calculated in "atomic" units (default 1 unit = 3000 bp).

4. **The Solution (DAGchainer)**
   - Uses Dynamic Programming to find the optimal path.
   - For every point, it calculates: "If I extend a chain from a previous point, what is the best score I can get?"
     - Current_Score = Match_Score + max(Previous_Score - Gap_Penalty).
   - This effectively balances "more matches" vs "tight collinearity".

# Standard format
hnsm synt dag pos.tsv match.tsv
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .index(1)
                .num_args(2)
                .help("Input files: 1. Position file, 2. Match file."),
        )
        .arg(
            Arg::new("go")
                .long("go")
                .num_args(1)
                .default_value("-1.0")
                .value_parser(value_parser!(f32))
                .help("Gap opening penalty"),
        )
        .arg(
            Arg::new("ge")
                .long("ge")
                .num_args(1)
                .default_value("-5.0")
                .value_parser(value_parser!(f32))
                .help("Gap extension penalty"),
        )
        .arg(
            Arg::new("bgs")
                .long("bgs")
                .num_args(1)
                .default_value("10000")
                .value_parser(value_parser!(i32))
                .help("Bp gap size"),
        )
        .arg(
            Arg::new("mms")
                .long("mms")
                .num_args(1)
                .default_value("50.0")
                .value_parser(value_parser!(f32))
                .help("Max match score"),
        )
        .arg(
            Arg::new("mdm")
                .long("mdm")
                .num_args(1)
                .default_value("100000")
                .value_parser(value_parser!(i32))
                .help("Max distance between matches"),
        )
        .arg(
            Arg::new("mas")
                .long("mas")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help("Min alignment score"),
        )
        .arg(
            Arg::new("mna")
                .long("mna")
                .num_args(1)
                .default_value("6")
                .value_parser(value_parser!(i32))
                .help("Min number of aligned pairs"),
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
    let infiles: Vec<_> = args.get_many::<String>("infiles").unwrap().collect();

    let opt_go = *args.get_one::<f32>("go").unwrap();
    let opt_ge = *args.get_one::<f32>("ge").unwrap();
    let opt_bgs = *args.get_one::<i32>("bgs").unwrap();
    let opt_mms = *args.get_one::<f32>("mms").unwrap();
    let opt_mdm = *args.get_one::<i32>("mdm").unwrap();

    let opt_mna = *args.get_one::<i32>("mna").unwrap();
    let opt_mas = if args.contains_id("mas") {
        *args.get_one::<f32>("mas").unwrap()
    } else {
        opt_mna as f32 * 0.5 * opt_mms
    };

    let chain_opt = ChainOpt {
        gap_open_penalty: opt_go,
        gap_extension_penalty: opt_ge,
        bp_gap_size: opt_bgs,
        max_match_score: opt_mms,
        max_dist_between_matches: opt_mdm,
        min_alignment_score: opt_mas,
    };

    //----------------------------
    // Ops
    //----------------------------
    let pos_file = &infiles[0];
    let match_file = &infiles[1];
    let acc_info = read_positions(pos_file)?;
    let (acc_pair_map, mol_pair_map) = parse_match_file(match_file, &acc_info, &chain_opt)?;
    // eprintln!("{:#?}", mol_pair_map);

    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut writer: Box<dyn Write> = if outfile == "stdout" {
        Box::new(io::BufWriter::new(io::stdout()))
    } else {
        Box::new(io::BufWriter::new(std::fs::File::create(outfile)?))
    };

    let mut alignment_count = 0;
    for mol_pair in mol_pair_map.keys() {
        let mut scores = vec![];

        for acc_pair in mol_pair_map.get(mol_pair).unwrap() {
            scores.push(MatchPair {
                pair_key: acc_pair.clone(),
                x: acc_info.get(&acc_pair.0).unwrap().mid as i32,
                y: acc_info.get(&acc_pair.1).unwrap().mid as i32,
                score: *acc_pair_map.get(acc_pair).unwrap(),
            });
        }
        // Remove score entries have the same identities.
        let mut scores = scores
            .iter()
            .unique_by(|e| (e.x, e.y))
            .cloned()
            .collect::<Vec<_>>();

        // The maximum column value of any match. Used to adjust coords for reverse diagonals
        let mut max_y = 0;
        if let Some(value) = scores.iter().map(|e| e.y).max_by(|a, b| a.cmp(b)) {
            max_y = value;
        }
        // eprintln!("scores = {:#?}", scores);

        // Sort scores by x coordinate, then y coordinate
        // This is CRITICAL for the DP algorithm which assumes topological order (or at least X order)
        scores.sort_by(|a, b| {
            if a.x == b.x {
                a.y.cmp(&b.y)
            } else {
                a.x.cmp(&b.x)
            }
        });

        // Convert to Anchors for the library
        let anchors: Vec<Anchor> = scores
            .iter()
            .enumerate()
            .map(|(i, s)| Anchor {
                id: i,
                x: s.x,
                y: s.y,
                score: s.score,
            })
            .collect();

        // Run DagChainer
        let chainer = DagChainer::new(chain_opt.clone());
        let chains = chainer.find_chains(&anchors);

        for chain in chains {
            // Check minimum number of pairs
            if chain.indices.len() < opt_mna as usize {
                continue;
            }

            print_alignment(&mut writer, &scores, &chain, max_y, alignment_count)?;
            alignment_count += 1;
        }
    }
    eprintln!("chain_opt = {:#?}", chain_opt);

    Ok(())
}

#[derive(Debug, Clone)]
struct Feature {
    mol: String,
    acc: String,
    // end5: usize,
    // end3: usize,
    mid: usize,
}

fn store_acc_info(
    mol: &str,
    acc: &str,
    end5: usize,
    end3: usize,
    acc_info: &mut HashMap<String, Feature>,
) {
    if !acc_info.contains_key(acc) {
        // Calculate the midpoint
        let mid_pt = ((end5 + end3) as f64 / 2.0).round() as usize;
        // Create and populate the Feature struct
        let feature = Feature {
            mol: mol.to_string(),
            acc: acc.to_string(),
            // end5,
            // end3,
            mid: mid_pt,
        };

        // Store entry in the provided hashmap
        acc_info.insert(feature.acc.clone(), feature.clone());
    }
}

fn pair_key(
    acc_info: &HashMap<String, Feature>,
    acc_1: &str,
    acc_2: &str,
) -> ((String, String), (String, String)) {
    let acc_pair_key;
    let mol_pair_key;

    let ft_1 = acc_info.get(acc_1).unwrap();
    let ft_2 = acc_info.get(acc_2).unwrap();

    // Determine molecule pair order
    if ft_1.mol < ft_2.mol {
        mol_pair_key = (ft_1.mol.to_string(), ft_2.mol.to_string());
        acc_pair_key = (acc_1.to_string(), acc_2.to_string());
    } else if ft_1.mol == ft_2.mol {
        mol_pair_key = (ft_1.mol.to_string(), ft_2.mol.to_string());
        if ft_1.mid < ft_2.mid {
            acc_pair_key = (acc_1.to_string(), acc_2.to_string());
        } else {
            acc_pair_key = (acc_2.to_string(), acc_1.to_string());
        }
    } else {
        mol_pair_key = (ft_2.mol.to_string(), ft_1.mol.to_string());
        acc_pair_key = (acc_2.to_string(), acc_1.to_string());
    }

    (acc_pair_key, mol_pair_key)
}

fn read_positions(path: &str) -> anyhow::Result<HashMap<String, Feature>> {
    let file = std::fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut acc_info = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.trim().split('\t').collect();

        if parts.len() == 2 {
            // hnsm gff rg format: key \t mol(strand):start-end
            let acc = parts[0];
            let location = parts[1];

            // Parse location: "Human.chr1(+):100-200" or "chr1:100-200"
            // Split by last ':' to separate mol+strand from range
            if let Some((mol_strand, range)) = location.rsplit_once(':') {
                // Parse range "100-200"
                if let Some((start_s, end_s)) = range.split_once('-') {
                    let start: usize = start_s.parse()?;
                    let end: usize = end_s.parse()?;

                    // Parse mol and strand from "Human.chr1(+)"
                    // If (...) is present, remove it.
                    let mol = if let Some(idx) = mol_strand.find('(') {
                        &mol_strand[..idx]
                    } else {
                        mol_strand
                    };

                    store_acc_info(mol, acc, start, end, &mut acc_info);
                }
            }
        } else if parts.len() >= 4 {
            // Standard 4-column format
            let mol = parts[0];
            let acc = parts[1];
            let start: usize = parts[2].parse()?;
            let end: usize = parts[3].parse()?;
            store_acc_info(mol, acc, start, end, &mut acc_info);
        }
    }
    Ok(acc_info)
}

fn parse_match_file(
    file_path: &str,
    acc_info: &HashMap<String, Feature>,
    opt: &ChainOpt,
) -> anyhow::Result<(
    HashMap<(String, String), f32>,
    HashMap<(String, String), Vec<(String, String)>>,
)> {
    let file = std::fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);

    let mut acc_pair_map: HashMap<(String, String), f32> = HashMap::new();
    let mut mol_pair_map: HashMap<(String, String), Vec<(String, String)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let acc_1 = parts[0];
        let acc_2 = parts[1];
        let raw_score: f64 = parts[2].parse().unwrap_or(0.0);

        // Filter self
        if acc_1 == acc_2 {
            continue;
        }

        // Calculate score
        // Only consider Similarity (0.0 - 1.0)
        // Scale to 0 - max_match_score
        let score = (raw_score as f32) * opt.max_match_score;

        if !acc_info.contains_key(acc_1) || !acc_info.contains_key(acc_2) {
            continue;
        }

        let (acc_pair_key, mol_pair_key) = pair_key(acc_info, acc_1, acc_2);

        if acc_pair_map.contains_key(&acc_pair_key) {
            let prev = acc_pair_map.get_mut(&acc_pair_key).unwrap();
            if *prev < score {
                *prev = score;
            }
        } else {
            acc_pair_map.insert(acc_pair_key.clone(), score);
        }

        mol_pair_map
            .entry(mol_pair_key)
            .or_default()
            .push(acc_pair_key);
    }

    for mol_pair in mol_pair_map.keys().cloned().collect::<Vec<_>>() {
        let value = mol_pair_map.get_mut(&mol_pair).unwrap();
        *value = value.iter().unique().cloned().collect();
    }

    Ok((acc_pair_map, mol_pair_map))
}

#[derive(Debug, Clone)]
struct MatchPair {
    pair_key: (String, String),
    x: i32,
    y: i32,
    score: f32,
}

fn print_alignment<W: Write>(
    writer: &mut W,
    scores: &[MatchPair],
    chain: &Chain,
    _max_y: i32,
    alignment_count: usize,
) -> io::Result<()> {
    writeln!(
        writer,
        "> Alignment #{} score = {:.1}",
        alignment_count + 1,
        chain.score
    )?;
    for (i, &index) in chain.indices.iter().enumerate() {
        let print_y = scores[index].y;
        // Note: reverse_order logic removed as it was always false in original code,
        // but if we want to support it, we need to pass a flag.
        // Assuming forward order for now as per original execute function.

        writeln!(
            writer,
            "{}\t{},{}\t{}\t{}\t{:7.1}\t{:7.1}",
            i,
            scores[index].pair_key.0,
            scores[index].pair_key.1,
            scores[index].x,
            print_y,
            scores[index].score,
            chain.path_scores[i],
        )?;
    }
    Ok(())
}
