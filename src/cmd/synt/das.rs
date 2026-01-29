use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("das")
        .about("Domain architecture similarity")
        .after_help(
            r###"
Domain architecture similarity via dynamic programming


cargo run --bin hnsm synt das 1 --sep ""

cargo run --bin hnsm synt das 1 --sep "" --ma 2 --mm=-1.0 --gp=-1.0


"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("ma")
                .long("ma")
                .num_args(1)
                .default_value("1.0")
                .value_parser(value_parser!(f32))
                .help("Match score"),
        )
        .arg(
            Arg::new("mm")
                .long("mm")
                .num_args(1)
                .default_value("-0.2")
                .value_parser(value_parser!(f32))
                .help("Mismatch score"),
        )
        .arg(
            Arg::new("gp")
                .long("gp")
                .num_args(1)
                .default_value("-0.01")
                .value_parser(value_parser!(f32))
                .help("Gap penalty"),
        )
        .arg(
            Arg::new("sep")
                .long("sep")
                .num_args(1)
                .default_value("\t")
                .help("Gap penalty"),
        )
        .arg(
            Arg::new("header")
                .long("header")
                .short('H')
                .action(ArgAction::SetTrue)
                .help("Inputs have header lines"),
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

#[derive(Debug, Clone)]
struct DasOpt {
    ma: f32,
    mm: f32,
    gp: f32,
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    // let infile = args.get_one::<String>("infile").unwrap();

    let das_opt = DasOpt {
        ma: *args.get_one::<f32>("ma").unwrap(),
        mm: *args.get_one::<f32>("mm").unwrap(),
        gp: *args.get_one::<f32>("gp").unwrap(),
    };

    let opt_sep = args.get_one::<String>("sep").unwrap();
    // let is_header = args.get_flag("header");
    //
    // let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Ops
    //----------------------------
    let a = "CATGT";
    let b = "ACGCTG";

    let s = a.split(opt_sep).map(str::to_string).collect::<Vec<_>>();
    let t = b.split(opt_sep).map(str::to_string).collect::<Vec<_>>();

    let mat = sim_mat(&s, &t, &das_opt);
    // eprintln!("mat = {:#?}", mat);
    let (sa, ta) = align(&mat, &s, &t, &das_opt);
    println!("{}", sa.join(opt_sep));
    println!("{}", ta.join(opt_sep));

    // let reader = intspan::reader(infile);
    // let mut lines = reader.lines().filter_map(Result::ok);
    // while let Some(line) = lines.next() {
    //     // Read sequence
    //     let (sh, th, s, t) = if is_header {
    //         let header1 = line;
    //         let seq1 = lines
    //             .next()
    //             .unwrap()
    //             .split(opt_sep)
    //             .map(str::to_string)
    //             .collect::<Vec<_>>();
    //         let header2 = lines.next().unwrap();
    //         let seq2 = lines
    //             .next()
    //             .unwrap()
    //             .split(opt_sep)
    //             .map(str::to_string)
    //             .collect::<Vec<_>>();
    //         (header1, header2, seq1, seq2)
    //     } else {
    //         let seq1 = lines
    //             .next()
    //             .unwrap()
    //             .split(opt_sep)
    //             .map(str::to_string)
    //             .collect::<Vec<_>>();
    //         let seq2 = lines
    //             .next()
    //             .unwrap()
    //             .split(opt_sep)
    //             .map(str::to_string)
    //             .collect::<Vec<_>>();
    //         (String::new(), String::new(), seq1, seq2)
    //     };
    //
    //     let mat = similarity(&s, &t, &das_opt);
    //     let (sa, ta) = align(&mat, &s, &t, &das_opt);
    //
    //     // if output_domain_distance {
    //     //     let dd = sa.iter().filter(|&&x| x == "-").count();
    //     //     if has_header {
    //     //         println!("{}\t{}\t{}", sh, th, dd);
    //     //     } else {
    //     //         println!("{}", dd);
    //     //     }
    //     // } else {
    //     if is_header {
    //         println!("{}", sh);
    //         println!("{}", sa.join(opt_sep));
    //         println!("{}", th);
    //         println!("{}", ta.join(opt_sep));
    //     } else {
    //         println!("{}", sa.join(opt_sep));
    //         println!("{}", ta.join(opt_sep));
    //     }
    //     // }
    // }

    Ok(())
}

fn compare(c1: &str, c2: &str, das_opt: &DasOpt) -> f32 {
    if c1 == c2 {
        das_opt.ma
    } else {
        das_opt.mm
    }
}

fn sim_mat(s: &[String], t: &[String], das_opt: &DasOpt) -> Vec<Vec<f32>> {
    let m = s.len();
    let n = t.len();
    let mut mat = vec![vec![0.0; n + 1]; m + 1];

    for i in 0..=m {
        mat[i][0] = das_opt.gp * i as f32;
    }
    for j in 0..=n {
        mat[0][j] = das_opt.gp * j as f32;
    }

    for i in 1..=m {
        for j in 1..=n {
            let p = compare(&s[i - 1], &t[j - 1], das_opt);
            mat[i][j] = 0.0_f32
                .max(mat[i - 1][j] + das_opt.gp)
                .max(mat[i][j - 1] + das_opt.gp)
                .max(mat[i - 1][j - 1] + p);
        }
    }

    mat
}

fn align(
    mat: &[Vec<f32>],
    s: &[String],
    t: &[String],
    das_opt: &DasOpt,
) -> (Vec<String>, Vec<String>) {
    let (mut sa, mut ta) = (Vec::new(), Vec::new());
    let (mut i, mut j) = (s.len(), t.len());

    while i != 0 || j != 0 {
        if i == 0 {
            // Case 3: last element of the 2nd array is paired with a gap
            sa.push("-".to_string());
            ta.push(t[j - 1].to_string());
            j -= 1;
        } else if j == 0 {
            // Case 2: last element of the 1st array is paired with a gap
            sa.push(s[i - 1].to_string());
            ta.push("-".to_string());
            i -= 1;
        } else {
            let p = compare(&s[i - 1], &t[j - 1], das_opt);

            if mat[i][j] == mat[i - 1][j - 1] + p {
                sa.push(s[i - 1].to_string());
                ta.push(t[j - 1].to_string());
                i -= 1;
                j -= 1;
            } else if mat[i - 1][j] > mat[i][j - 1] {
                sa.push(s[i - 1].to_string());
                ta.push("-".to_string());
                i -= 1;
            } else {
                sa.push("-".to_string());
                ta.push(t[j - 1].to_string());
                j -= 1;
            }
        }
    }

    let mut ret_s = vec![];
    let mut ret_t = vec![];
    while let Some(s) = sa.pop() {
        ret_s.push(s);
    }
    while let Some(t) = ta.pop() {
        ret_t.push(t);
    }
    (ret_s, ret_t)
}
