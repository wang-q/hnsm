use clap::*;

use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{self, BufRead};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("Chains of syntenic genes")
        .after_help(
            r###"
Algorithm adopted from `DAGchainer`

cat ~/Scripts/DAGCHAINER/data_sets/Arabidopsis/Arabidopsis.Release5.matchList.filtered |
    tsv-filter --eq 1:1 --eq 5:2 \
    > ath-1-2.tsv

cargo run --bin hnsm chain ath-1-2.tsv


"###,
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

#[derive(Debug)]
pub struct ChainOpt {
    gap_open_penalty: f32,
    gap_extension_penalty: f32,
    bp_gap_size: i32,
    max_match_score: f32,
    max_dist_between_matches: i32,
    min_alignment_score: f32,
    reverse_order: bool,
    max_y: i32,
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

    let mut chain_opt = ChainOpt {
        gap_open_penalty: opt_go,
        gap_extension_penalty: opt_ge,
        bp_gap_size: opt_bgs,
        max_match_score: opt_mms,
        max_dist_between_matches: opt_mdm,
        min_alignment_score: opt_mas,
        reverse_order: false,
        max_y: 0,
    };

    //----------------------------
    // Ops
    //----------------------------
    let (acc_info, acc_pair_map, mol_pair_map) = parse_input_file(infile, &chain_opt)?;
    // eprintln!("{:#?}", mol_pair_map);

    for mol_pair in mol_pair_map.keys() {
        let mut scores = vec![];

        for acc_pair in mol_pair_map.get(mol_pair).unwrap() {
            scores.push(Score {
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
        if let Some(value) = scores.iter().map(|e| e.y).max_by(|a, b| a.cmp(b)) {
            chain_opt.max_y = value;
        }
        // eprintln!("scores = {:#?}", scores);
        print_chains(&mut scores, &chain_opt);
    }
    eprintln!("chain_opt = {:#?}", chain_opt);

    Ok(())
}

#[derive(Debug, Clone)]
struct Feature {
    mol: String,
    acc: String,
    end5: usize,
    end3: usize,
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
            end5,
            end3,
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
    let mut acc_pair_key = (String::new(), String::new());
    let mut mol_pair_key = (String::new(), String::new());

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

fn parse_input_file(
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
        // if mol_1 == mol_2 & &!include_self {
        //     continue;
        // }
        // if mol_1 != mol_2 & &tandem_only {
        //     continue;
        // }
        let score = scoring_f(score, opt.max_match_score);

        // Handle features
        store_acc_info(mol_1, acc_1, end5_1, end3_1, &mut acc_info);
        store_acc_info(mol_2, acc_2, end5_2, end3_2, &mut acc_info);

        let (acc_pair_key, mol_pair_key) = pair_key(&acc_info, acc_1, acc_2);
        if acc_pair_map.contains_key(&acc_pair_key) {
            let prev = acc_pair_map.get_mut(&acc_pair_key).unwrap();
            // take the lowest e_value for the match pair
            // Important when there are multiple HSPs reported between two accessions
            if *prev < score {
                *prev = score;
            }
        } else {
            acc_pair_map.insert(acc_pair_key.clone(), score);
        }

        mol_pair_map
            .entry(mol_pair_key)
            .or_insert(Vec::new())
            .push(acc_pair_key);
    }

    for mol_pair in mol_pair_map.keys().cloned().collect::<Vec<_>>() {
        let value = mol_pair_map.get_mut(&mol_pair).unwrap();
        *value = value.iter().unique().cloned().collect();
    }

    Ok((acc_info, acc_pair_map, mol_pair_map))
}

#[derive(Debug, Clone)]
struct Score {
    pair_key: (String, String),
    x: i32,
    y: i32,
    score: f32,
}

#[derive(Debug, Default)]
struct Path {
    score: f32,
    rc: i32,
    sub: usize,
}

// reverse complement the second coordinate set.
fn adjust_scores(mut scores: Vec<Score>) -> anyhow::Result<Vec<Score>> {
    // if unsafe { REVERSE_ORDER } {
    //     for score in scores.iter_mut() {
    //         score.y = unsafe { MAX_Y - score.y + 1 };
    //     }
    // }
    scores.sort_by(|a, b| {
        if a.x == b.x {
            a.y.cmp(&b.y)
        } else {
            a.x.cmp(&b.x)
        }
    });
    Ok(scores)
}

//  Find and output highest scoring chains in scores treating it as a DAG
fn print_chains(scores: &mut Vec<Score>, options: &ChainOpt) {
    loop {
        let mut updated = false;

        // Initialize path scores and 'from' indices
        let n = scores.len();
        let mut path_scores = vec![0.0; n];
        let mut from_indices = vec![-1; n];
        for i in 0..n {
            path_scores[i] = scores[i].score;
            from_indices[i] = -1_i32;
        }

        for j in 1..n {
            for i in (0..j).rev() {
                let del_x = scores[j].x - scores[i].x - 1;
                let del_y = scores[j].y - scores[i].y - 1;

                if del_x < 0 || del_y < 0 {
                    continue;
                }

                // Check maximum distances
                if del_x > options.max_dist_between_matches
                    && del_y > options.max_dist_between_matches
                {
                    break;
                }
                if del_x > options.max_dist_between_matches
                    || del_y > options.max_dist_between_matches
                {
                    continue;
                }

                let num_gaps = ((del_x + del_y + (del_x - del_y).abs()) as f32
                    / (2 * options.bp_gap_size) as f32
                    + 0.5) as i32;
                let mut new_score = path_scores[i] + scores[j].score;

                if num_gaps > 0 {
                    new_score += options.gap_open_penalty
                        + (num_gaps as f32 * options.gap_extension_penalty);
                }

                if new_score > path_scores[j] {
                    path_scores[j] = new_score;
                    from_indices[j] = i as i32;
                    updated = true;
                }
            }
        }

        let high_scores: Vec<Path> = path_scores
            .iter()
            .enumerate()
            .filter(|&(_, &score)| score >= options.min_alignment_score)
            .map(|(sub, &score)| Path {
                score,
                sub,
                rc: scores[sub].x + scores[sub].y,
            })
            .collect();

        let mut high: Vec<Path> = high_scores;
        high.sort_by(|a, b| {
            if a.score != b.score {
                a.score
                    .partial_cmp(&b.score)
                    .unwrap_or(Ordering::Equal)
                    .reverse()
            } else {
                a.rc.cmp(&b.rc)
            }
        });

        let mut ali_ct = 0;
        for entry in high {
            if from_indices[entry.sub] != -1 {
                let alignment_path = build_alignment_path(&from_indices, entry.sub);
                print_alignment(&scores, &path_scores, alignment_path, options, ali_ct);
                ali_ct += 1;
            }
        }

        if !updated {
            break;
        }

        // Retain only updated scores
        let mut index = 0;
        scores.retain(|_| {
            index += 1;
            from_indices[index - 1] != -1
        });
    }
}

fn build_alignment_path(from_indices: &Vec<i32>, start_index: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut current = start_index;

    while from_indices[current] >= 0 {
        path.push(current);
        current = from_indices[current] as usize;
    }
    path.push(current); // Include the start path.reverse();
    path
}

fn print_alignment(
    scores: &Vec<Score>,
    path_scores: &Vec<f32>,
    path: Vec<usize>,
    options: &ChainOpt,
    alignment_count: usize,
) {
    println!(
        "> Alignment #{} score = {:.1}",
        alignment_count + 1,
        path_scores[*path.first().unwrap()]
    );
    for &index in &path {
        let print_y = if options.reverse_order {
            options.max_y - scores[index].y + 1
        } else {
            scores[index].y
        };
        println!(
            "{}\t{},{}\t{}\t{}\t{:7.1}\t{:7.1}",
            path.iter().position(|&x| x == index).unwrap(),
            scores[index].pair_key.0,
            scores[index].pair_key.1,
            scores[index].x,
            print_y,
            scores[index].score,
            path_scores[index],
        );
    }
}
