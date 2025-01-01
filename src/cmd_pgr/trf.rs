use clap::*;
use cmd_lib::*;
use std::io::BufRead;
use std::{env, vec};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("trf")
        .about("Identify tandem repeats in a genome")
        .after_help(
            r###"
This command identifies tandem repeats in a genome via `trf`.

* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * trf
    * spanr

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input file to process."),
        )
        .arg(
            Arg::new("match")
                .long("match")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("Matching weight"),
        )
        .arg(
            Arg::new("mismatch")
                .long("mismatch")
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("Mismatching penalty"),
        )
        .arg(
            Arg::new("delta")
                .long("delta")
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("Indel penalty"),
        )
        .arg(
            Arg::new("pm")
                .long("pm")
                .num_args(1)
                .default_value("80")
                .value_parser(value_parser!(usize))
                .help("Match probability"),
        )
        .arg(
            Arg::new("pi")
                .long("pi")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Indel probability"),
        )
        .arg(
            Arg::new("minscore")
                .long("minscore")
                .num_args(1)
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("Minimum alignment score to report"),
        )
        .arg(
            Arg::new("maxperiod")
                .long("maxperiod")
                .num_args(1)
                .default_value("2000")
                .value_parser(value_parser!(usize))
                .help("Maximum period size to report"),
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
    let outfile = args.get_one::<String>("outfile").unwrap();

    let opt_match = *args.get_one::<usize>("match").unwrap();
    let opt_mismatch = *args.get_one::<usize>("mismatch").unwrap();
    let opt_delta = *args.get_one::<usize>("delta").unwrap();
    let opt_pm = *args.get_one::<usize>("pm").unwrap();
    let opt_pi = *args.get_one::<usize>("pi").unwrap();
    let opt_minscore = *args.get_one::<usize>("minscore").unwrap();
    let opt_maxperiod = *args.get_one::<usize>("maxperiod").unwrap();

    //----------------------------
    // Paths
    //----------------------------
    let curdir = env::current_dir()?;
    let pgr = env::current_exe()?.display().to_string();
    let tempdir = tempfile::Builder::new().prefix("pgr_trf_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"pgr\"     = ${pgr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_infile = intspan::absolute_path(args.get_one::<String>("infile").unwrap())?
        .display()
        .to_string();
    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        intspan::absolute_path(outfile)?.display().to_string()
    };

    //----------------------------
    // Ops
    //----------------------------
    run_cmd!(info "==> Switch to outdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Split by names")?;
    run_cmd!(
        hnsm split name ${abs_infile} -o .
    )?;

    run_cmd!(info "==> Process each chromosomes")?;
    run_cmd!(
        hnsm size ${abs_infile} -o chr.sizes
    )?;
    let mut chrs: Vec<String> = vec![];
    for line in intspan::read_lines("chr.sizes") {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            chrs.push(fields[0].to_string());
        }
    }

    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        run_cmd!(
            trf ${chr}.fa ${opt_match} ${opt_mismatch} ${opt_delta} ${opt_pm} ${opt_pi} ${opt_minscore} ${opt_maxperiod} -d -h -ngs > trf.${i}.dat
        )?;

        // 198 229 12 2.7 12 90 0 50 34 46 3 15 1.62 CATTACCACCAC CATTAGCACCACCATTACCACCACCATCACCA ATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACG TTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAA
        // start end
        // period copy_number consensus_pattern_size
        // perc_matches perc_indels
        // alignment_score
        // perc_a perc_c perc_g perc_t
        // entropy
        // consensus_pattern
        // repeat_seq
        // 15 fields
        // The last 2 fields were introduced by -ngs
        // Matched with `hnsm range mg1655.fa NC_000913:198-229`

        let reader = intspan::reader(&format!("trf.{}.dat", i));

        let rg_file = format!("trf.{}.rg", i);
        let mut writer = intspan::writer(&rg_file);
        for line in reader.lines().map_while(Result::ok) {
            let fields: Vec<&str> = line.split_ascii_whitespace().collect();
            if fields.len() < 15 {
                continue;
            }

            let start = fields[0].parse::<usize>()?;
            let end = fields[1].parse::<usize>()?;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }
        rg_files.push(rg_file);
    }

    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        spanr cover $[rg_files] -o ${abs_outfile}
    )?;

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}

// use std::io::{Read, Write};
// fn pause() {
//     let mut stdin = std::io::stdin();
//     let mut stdout = std::io::stdout();
//
//     // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
//     write!(stdout, "Press any key to continue...").unwrap();
//     stdout.flush().unwrap();
//
//     // Read a single byte and discard
//     let _ = stdin.read(&mut [0u8]).unwrap();
// }
