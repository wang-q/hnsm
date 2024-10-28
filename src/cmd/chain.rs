use clap::*;

use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead};
use std::path;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("Chains of syntenic genes")
        .after_help(
            r###"
Algorithm adopted from `DAGchainer`

cat ~/Scripts/DAGCHAINER/data_sets/Arabidopsis/Arabidopsis.Release5.matchList.filtered |
    tsv-filter --eq 1:1 --eq 5:1 \
    > ath-1-1.tsv

cargo run --bin hnsm chain ath-1-1.tsv


"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("eps")
                .long("eps")
                .num_args(1)
                .default_value("0.05")
                .value_parser(value_parser!(f32))
                .help("The maximum distance between two points"),
        )
        .arg(
            Arg::new("min_points")
                .long("min_points")
                .num_args(1)
                .default_value("1")
                .value_parser(value_parser!(usize))
                .help("core point"),
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

const MAX_MATCH_SCORE: f64 = 50.0;

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();

    let (acc_info, acc_pair, mol_pair) = parse_input_file(infile)?;
    eprintln!("acc_pair = {:#?}", mol_pair);

    Ok(())
}

#[derive(Debug)]
pub struct ChainOpt {
    indel_score: f32,
    // gap extension penalty
    gap_open_penalty: f32,
    // gap open penalty
    bp_gap_size: i32,
    // length for a single gap in base pairs
    min_alignment_score: i32,
    max_dist_between_matches: i32,
    //  If set true by -r option, then use 2nd coordinates in reverse order
    reverse_order: bool,
    max_y: i32, // The maximum column value of any match. Used to adjust coords for reverse diagonals
}

impl Default for ChainOpt {
    fn default() -> ChainOpt {
        ChainOpt {
            indel_score: 1.0,
            gap_open_penalty: 1.0,
            bp_gap_size: -1,
            min_alignment_score: -1,
            max_dist_between_matches: 100_000,
            reverse_order: false,
            max_y: 0,
        }
    }
}

// //
// static mut INDEL_SCORE: f32 = 1.0;
//
// static mut GAP_OPEN_PENALTY: f32 = 1.0;
// static mut BP_GAP_SIZE: i32 = -1;
// static mut MIN_ALIGNMENT_SCORE: i32 = -1;
// static mut MAX_DIST_BETWEEN_MATCHES: i32 = 100_000;
//
// static mut REVERSE_ORDER: bool = false;
// static mut MAX_Y: i32 = 0;

// reading in coordinate file, add each match to a list.
fn read_scores(filename: &str) -> anyhow::Result<Vec<Score>> {
    let path = path::Path::new(filename);
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut scores = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        let pair_id = parts.next().unwrap().parse::<i32>().unwrap();
        let x = parts.next().unwrap().parse::<i32>().unwrap();
        let y = parts.next().unwrap().parse::<i32>().unwrap();
        let score = parts.next().unwrap().parse::<f32>().unwrap();

        // if y > unsafe { MAX_Y } {
        //     unsafe {
        //         MAX_Y = y;
        //     }
        // }
        scores.push(Score {
            pair_id,
            x,
            y,
            score,
        });
    }

    Ok(scores)
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

fn pair_key(acc_info: &HashMap<String, Feature>, acc_1: &str, acc_2: &str) -> (String, String) {
    let mut acc_pair_key = String::new();
    let mut mol_pair_key = String::new();

    let ft_1 = acc_info.get(acc_1).unwrap();
    let ft_2 = acc_info.get(acc_2).unwrap();

    // Determine molecule pair order
    if ft_1.mol < ft_2.mol {
        mol_pair_key = format!("{},{}", ft_1.mol, ft_2.mol);
        acc_pair_key = format!("{},{}", acc_1, acc_2);
    } else if ft_1.mol == ft_2.mol {
        mol_pair_key = format!("{},{}", ft_1.mol, ft_2.mol);
        if ft_1.mid < ft_2.mid {
            acc_pair_key = format!("{},{}", acc_1, acc_2);
        } else {
            acc_pair_key = format!("{},{}", acc_2, acc_1);
        }
    } else {
        mol_pair_key = format!("{},{}", ft_2.mol, ft_1.mol);
        acc_pair_key = format!("{},{}", acc_2, acc_1);
    }

    (acc_pair_key, mol_pair_key)
}

fn scoring_f(evalue: f64) -> f64 {
    let match_score = -evalue.log10() * 10.0;
    let rounded_score = (match_score + 0.5).floor() / 10.0; // Round to one decimal place
    rounded_score.min(MAX_MATCH_SCORE) // Ensure it does not exceed MAX_MATCH_SCORE
}

fn parse_input_file(
    file_path: &str,
) -> anyhow::Result<(
    HashMap<String, Feature>,
    HashMap<String, f64>,
    HashMap<String, Vec<String>>,
)> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    let mut acc_info: HashMap<String, Feature> = HashMap::new();
    let mut acc_pair: HashMap<String, f64> = HashMap::new();
    let mut mol_pair: HashMap<String, Vec<String>> = HashMap::new();

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
        let score = scoring_f(score);

        // Handle features
        store_acc_info(mol_1, acc_1, end5_1, end3_1, &mut acc_info);
        store_acc_info(mol_2, acc_2, end5_2, end3_2, &mut acc_info);

        let (acc_pair_key, mol_pair_key) = pair_key(&acc_info, acc_1, acc_2);
        if acc_pair.contains_key(&acc_pair_key) {
            let prev = acc_pair.get_mut(&acc_pair_key).unwrap();
            // take the lowest e_value for the match pair
            // Important when there are multiple HSPs reported between two accessions
            if *prev < score {
                *prev = score;
            }
        } else {
            acc_pair.insert(acc_pair_key.clone(), score);
        }

        mol_pair
            .entry(mol_pair_key)
            .or_insert(Vec::new())
            .push(acc_pair_key);
    }

    Ok((acc_info, acc_pair, mol_pair))
}


#[derive(Debug, Clone)]
struct Score {
    pair_id: i32,
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

// // reverse complement the second coordinate set.
// fn adjust_scores(mut scores: Vec<Score>) -> anyhow::Result<Vec<Score>> {
//     if unsafe { REVERSE_ORDER } {
//         for score in scores.iter_mut() {
//             score.y = unsafe { MAX_Y - score.y + 1 };
//         }
//     }
//     scores.sort_by(|a, b| {
//         if a.x == b.x {
//             a.y.cmp(&b.y)
//         } else {
//             a.x.cmp(&b.x)
//         }
//     });
//     Ok(scores)
// }
//
// // consecutive score entries have the same identities.
// fn remove_duplicates(scores: Vec<Score>) -> Vec<Score> {
//     let mut unique_scores = Vec::new();
//     let mut previous_score: Option<Score> = None;
//
//     for score in scores.iter() {
//         match previous_score {
//             Some(ref prev) if prev.x == score.x && prev.y == score.y => {
//                 eprintln!("Duplicate score entries: {:?} and {:?}", prev, score);
//                 eprintln!("Discarding latter");
//                 continue;
//             }
//             _ => {
//                 unique_scores.push(score.clone());
//             }
//         }
//         previous_score = Some(score.clone());
//     }
//
//     unique_scores
// }
//
// //  Find and output highest scoring chains in  score  treating it as a DAG
// fn print_chains(score: &mut Vec<Score>) {
//     let mut path_score = vec![0.0; score.len()];
//     let mut from = vec![-1; score.len()];
//     let mut high: Vec<Path> = Vec::new();
//     let mut done: bool;
//     let mut ali_ct = 0;
//
//     loop {
//         done = true;
//         let n = score.len();
//         for i in 0..n {
//             path_score[i] = score[i].score;
//             from[i] = -1;
//         }
//
//         for j in 1..n {
//             for i in (0..j).rev() {
//                 let del_x = score[j].x - score[i].x - 1;
//                 let del_y = score[j].y - score[i].y - 1;
//
//                 if del_x >= 0 && del_y >= 0 {
//                     if del_x > unsafe { MAX_DIST_BETWEEN_MATCHES }
//                         && del_y > unsafe { MAX_DIST_BETWEEN_MATCHES }
//                     {
//                         break;
//                     }
//                     if del_x > unsafe { MAX_DIST_BETWEEN_MATCHES }
//                         || del_y > unsafe { MAX_DIST_BETWEEN_MATCHES }
//                     {
//                         continue;
//                     }
//
//                     let num_gaps = ((del_x + del_y + (del_x - del_y).abs())
//                         / (2 * unsafe { BP_GAP_SIZE }) as f32
//                         + 0.5) as i32;
//                     let mut x = path_score[i] + score[j].score;
//
//                     if num_gaps > 0 {
//                         x += unsafe { GAP_OPEN_PENALTY }
//                             + (num_gaps as f32 * unsafe { INDEL_SCORE });
//                     }
//
//                     if x > path_score[j] {
//                         path_score[j] = x;
//                         from[j] = i;
//                     }
//                 }
//             }
//         }
//
//         high.clear();
//         for i in 0..n {
//             if path_score[i] >= unsafe { MIN_ALIGNMENT_SCORE } as f32 {
//                 high.push(Path {
//                     score: path_score[i],
//                     sub: i,
//                     rc: score[i].x + score[i].y,
//                 });
//             }
//         }
//
//         high.sort_by(|a, b| {
//             if a.score > b.score {
//                 Ordering::Less
//             } else if a.score < b.score {
//                 Ordering::Greater
//             } else {
//                 a.rc.cmp(&b.rc)
//             }
//         });
//
//         let m = high.len();
//         for i in 0..m {
//             if from[high[i].sub] != -2 {
//                 let mut ans = Vec::new();
//                 let mut j = high[i].sub;
//                 while from[j] >= 0 {
//                     ans.push(j);
//                     j = from[j];
//                 }
//                 ans.push(j);
//                 if from[j] == -2 {
//                     done = false;
//                     break;
//                 } else {
//                     ans.reverse();
//                     let s = ans.len();
//                     println!(
//                         ">Alignment #{} score = {:.1}",
//                         ali_ct + 1,
//                         path_score[high[i].sub]
//                     );
//                     for j in 0..s {
//                         let print_y = if unsafe { REVERSE_ORDER } {
//                             unsafe { MAX_Y - score[ans[j]].y + 1 }
//                         } else {
//                             score[ans[j]].y
//                         };
//                         // Mark as printed
//                         println!(
//                             "{:3}: {} {:6} {:6} {:7.1} {:7.1}",
//                             j,
//                             score[ans[j]].pair_id,
//                             score[ans[j]].x,
//                             print_y,
//                             score[ans[j]].score,
//                             path_score[ans[j]]
//                         );
//                     }
//                 }
//             }
//         }
//
//         if !done {
//             let mut j = 0;
//             for i in 0..n {
//                 if from[i] != -2 {
//                     if i != j {
//                         score[j] = score[i].clone(); // Use clone to properly copy the value }
//                         j += 1;
//                     }
//                 }
//                 score.truncate(j);
//             }
//         }
//     }
// }
