use clap::*;
use noodles_fasta as fasta;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("n50")
        .about("Count total bases in FA file(s)")
        .after_help(
            r#"
* N50 is the default output value, set a single `-N 0` to skip this
* To calculate both N50 and N90, enter `-N 50 -N 90`
* Turn on other options to compute more statitics
* E-size is defined as the expected contig length at which a random position locates

"#,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("noheader")
                .long("noheader")
                .short('H')
                .action(ArgAction::SetTrue)
                .help("Do not display headers"),
        )
        .arg(
            Arg::new("nx")
                .long("nx")
                .short('N')
                .num_args(1)
                .default_value("50")
                .action(ArgAction::Append)
                .value_parser(value_parser!(usize))
                .help("Compute Nx statistic"),
        )
        .arg(
            Arg::new("sum")
                .long("sum")
                .short('S')
                .action(ArgAction::SetTrue)
                .help("Compute the sum of the sizes of all records"),
        )
        .arg(
            Arg::new("average")
                .long("average")
                .short('A')
                .action(ArgAction::SetTrue)
                .help("Compute the average length of all records"),
        )
        .arg(
            Arg::new("esize")
                .long("esize")
                .short('E')
                .action(ArgAction::SetTrue)
                .help("Compute the E-size (from GAGE)"),
        )
        .arg(
            Arg::new("count")
                .long("count")
                .short('C')
                .action(ArgAction::SetTrue)
                .help("Count records"),
        )
        .arg(
            Arg::new("genome")
                .long("genome")
                .short('g')
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Size of the genome, not the total size of the files"),
        )
        .arg(
            Arg::new("transpose")
                .long("transpose")
                .short('t')
                .action(ArgAction::SetTrue)
                .help("Transpose the outputs"),
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
    let is_noheader = args.get_flag("noheader");
    let is_sum = args.get_flag("sum");
    let is_average = args.get_flag("average");
    let is_esize = args.get_flag("esize");
    let is_count = args.get_flag("count");
    let is_transpose = args.get_flag("transpose");

    let opt_nx: Vec<_> = args
        .get_many::<usize>("nx")
        .unwrap()
        .map(|el| *el)
        .collect();

    let opt_genome = if args.contains_id("genome") {
        *args.get_one::<usize>("genome").unwrap()
    } else {
        usize::MAX
    };

    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Operating
    //----------------------------
    let mut lens = vec![];
    let mut record_cnt = 0;
    let mut total_size = 0;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);
        let mut fa_in = fasta::io::Reader::new(reader);

        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result?;

            let len = record.sequence().len();

            lens.push(len);
            record_cnt += 1;
            total_size += len;
        }
    }
    lens.sort_unstable_by(|a, b| b.cmp(a));
    // eprintln!("lens = {:#?}", lens);

    // reach n_given% of total_size or genome_size
    let mut goals = vec![];
    for el in opt_nx.iter() {
        let goal = if opt_genome != usize::MAX {
            (*el as f64) * (opt_genome as f64) / 100.0
        } else {
            (*el as f64) * (total_size as f64) / 100.0
        } as usize;
        goals.push(goal);
    }

    let mut cumul_size = 0; // the cumulative size
    let mut e_size = 0.0;
    let mut nx_sizes = vec![0; goals.len()];

    for cur_size in lens {
        let prev_cumul_size = cumul_size;
        cumul_size += cur_size;

        e_size = (prev_cumul_size as f64) / (cumul_size as f64) * e_size
            + (cur_size as f64 * cur_size as f64) / cumul_size as f64;

        for (i, goal) in goals.iter().enumerate() {
            if nx_sizes[i] == 0 && cumul_size > *goal {
                nx_sizes[i] = cur_size;
            }
        }
    }

    //----------------------------
    // Output
    //----------------------------
    let mut outputs = vec![];

    // set N == 0 to skip this
    if !(opt_nx.len() == 1 && opt_nx[0] == 0) {
        for (i, nx) in opt_nx.iter().enumerate() {
            let mut row = vec![];
            if !is_noheader {
                row.push(format!("N{}", nx));
            }
            row.push(format!("{}", nx_sizes[i]));
            outputs.push(row);
        }
    }

    if is_sum {
        let mut row = vec![];
        if !is_noheader {
            row.push(format!("S"));
        }
        row.push(format!("{}", total_size));
        outputs.push(row);
    }

    if is_average {
        let mut row = vec![];
        if !is_noheader {
            row.push(format!("A"));
        }
        row.push(format!("{:.2}", total_size as f64 / record_cnt as f64));
        outputs.push(row);
    }

    if is_esize {
        let mut row = vec![];
        if !is_noheader {
            row.push(format!("E"));
        }
        row.push(format!("{:.2}", e_size));
        outputs.push(row);
    }

    if is_count {
        let mut row = vec![];
        if !is_noheader {
            row.push(format!("C"));
        }
        row.push(format!("{:.2}", record_cnt));
        outputs.push(row);
    }

    if is_transpose {
        outputs = transpose(outputs);
    }

    for row in outputs {
        writer.write_fmt(format_args!("{}\n", row.join("\t")))?;
    }

    Ok(())
}

// https://stackoverflow.com/questions/64498617/how-to-transpose-a-vector-of-vectors-in-rust
fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}
