use clap::*;

use hnsm::libs::synteny::chain::{Anchor, Chain, ChainOpt, DagChainer};
use itertools::Itertools;
use std::collections::HashMap;
use std::io::{self, BufRead};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dag")
        .about("Synteny chains via DAGchainer algorithm")
        .after_help(
            r###"
Algorithm adopted from `DAGchainer`

# Legacy format
cat ~/Scripts/DAGCHAINER/data_sets/Arabidopsis/Arabidopsis.Release5.matchList.filtered |
    tsv-filter --eq 1:1 --eq 5:2 \
    > ath-1-2.tsv

hnsm synt dag ath-1-2.tsv

# Standard format
hnsm synt dag match.tsv --annot annot.tsv
"###,
        )
        .arg(
            Arg::new("annot")
                .long("annot")
                .short('a')
                .num_args(1)
                .help("Annotation file (tsv: mol, acc, start, end). If provided, infile is treated as a match file (tsv: acc1, acc2, score)."),
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
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
    let infile = args.get_one::<String>("infile").unwrap();

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
    let (acc_info, acc_pair_map, mol_pair_map) = if let Some(annot_file) = args.get_one::<String>("annot") {
        let acc_info = read_annotations(annot_file)?;
        let (acc_pair_map, mol_pair_map) = parse_match_file(infile, &acc_info, &chain_opt)?;
        (acc_info, acc_pair_map, mol_pair_map)
    } else {
        parse_legacy_input(infile, &chain_opt)?
    };
    // eprintln!("{:#?}", mol_pair_map);

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

            print_alignment(&scores, &chain, max_y, alignment_count);
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

fn scoring_f(evalue: f64, max_match_score: f32) -> f32 {
    let match_score = -evalue.log10() * 10.0;
    let rounded_score = (match_score + 0.5).floor() / 10.0; // Round to one decimal place
    rounded_score.min(max_match_score as f64) as f32 // Ensure it does not exceed MAX_MATCH_SCORE
}

fn parse_legacy_input(
    file_path: &str,
    opt: &ChainOpt,
) -> anyhow::Result<(
    HashMap<String, Feature>,
    HashMap<(String, String), f32>,
    HashMap<(String, String), Vec<(String, String)>>,
)> {
    let file = std::fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);

    let mut acc_info: HashMap<String, Feature> = HashMap::new();
    let mut acc_pair_map: HashMap<(String, String), f32> = HashMap::new();
    let mut mol_pair_map: HashMap<(String, String), Vec<(String, String)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        // Skip empty lines and those without word characters
        if line.is_empty() || !line.chars().any(|c| c.is_alphanumeric()) {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            continue; // Ensure there are enough fields
        }

        let mol_1 = parts[0];
        let acc_1 = parts[1];
        let end5_1: usize = parts[2].parse().unwrap();
        let end3_1: usize = parts[3].parse().unwrap();
        let mol_2 = parts[4];
        let acc_2 = parts[5];
        let end5_2: usize = parts[6].parse().unwrap();
        let end3_2: usize = parts[7].parse().unwrap();
        let mut score: f64 = parts[8].parse().unwrap();

        // Adjust e_value if it's too low
        if score < 1.0e-250 {
            score = 1.0e-250;
        }

        // Filtering records
        if acc_1 == acc_2 {
            continue; // No self comparisons
        }
        if score > 1.0e-5 {
            continue;
        }

        let score = scoring_f(score, opt.max_match_score);

        // Handle features
        store_acc_info(mol_1, acc_1, end5_1, end3_1, &mut acc_info);
        store_acc_info(mol_2, acc_2, end5_2, end3_2, &mut acc_info);

        let (acc_pair_key, mol_pair_key) = pair_key(&acc_info, acc_1, acc_2);
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

    Ok((acc_info, acc_pair_map, mol_pair_map))
}

fn read_annotations(path: &str) -> anyhow::Result<HashMap<String, Feature>> {
    let file = std::fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut acc_info = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() < 4 {
            continue;
        }

        let mol = parts[0];
        let acc = parts[1];
        let start: usize = parts[2].parse()?;
        let end: usize = parts[3].parse()?;

        store_acc_info(mol, acc, start, end, &mut acc_info);
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
        let mut score: f64 = parts[2].parse().unwrap_or(1.0); // Default score if parsing fails? Or expect valid score.

        // Adjust e_value if it's too low (assuming input is e-value if very small)
        // If input is bitscore (e.g. > 10), we should probably handle it differently.
        // Original logic assumes E-value.
        // Let's assume input is E-value for consistency with legacy, OR pre-calculated score.
        // If score is > 1.0 (bitscore?), we might want to take it as is.
        // But `scoring_f` expects E-value.
        // Let's assume input is E-value.
        if score < 1.0e-250 {
            score = 1.0e-250;
        }
        
        // Filter self
        if acc_1 == acc_2 {
            continue;
        }

        // Calculate score
        // If input is already a score (large positive), `scoring_f` might produce weird results if it expects small E-values.
        // `scoring_f`: -log10(evalue) * 10.
        // If input is 50.0 (bitscore), -log10(50) is negative.
        // We should probably check if score looks like an E-value or a Score.
        // For now, let's strictly follow legacy behavior: Input is E-value.
        let score = scoring_f(score, opt.max_match_score);

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

fn print_alignment(
    scores: &[MatchPair],
    chain: &Chain,
    _max_y: i32,
    alignment_count: usize,
) {
    println!(
        "> Alignment #{} score = {:.1}",
        alignment_count + 1,
        chain.score
    );
    for (i, &index) in chain.indices.iter().enumerate() {
        let print_y = scores[index].y; 
        // Note: reverse_order logic removed as it was always false in original code, 
        // but if we want to support it, we need to pass a flag. 
        // Assuming forward order for now as per original execute function.
        
        println!(
            "{}\t{},{}\t{}\t{}\t{:7.1}\t{:7.1}",
            i,
            scores[index].pair_key.0,
            scores[index].pair_key.1,
            scores[index].x,
            print_y,
            scores[index].score,
            chain.path_scores[i],
        );
    }
}
