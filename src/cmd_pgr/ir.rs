use clap::*;
use cmd_lib::*;
use std::io::BufRead;
use std::{env, vec};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("ir")
        .about("Identify interspersed repeats in a genome")
        .after_help(
            r###"
This command identifies interspersed repeats in a genome, mimicking the functionality of `RepeatMasker`.

* <repeat> is path to the fasta file containing repeats from Dfam, RepBase, or TnCentral。
* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * FastK / Profex / Fastrm
    * spanr

"###,
        )
        .arg(
            Arg::new("repeat")
                .required(true)
                .num_args(1)
                .index(1)
                .help("The repeats database"),
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Input file to process"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Size of the k-mer"),
        )
        .arg(
            Arg::new("fk")
                .long("fk")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("Fill holes between repetitive k-mers"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .num_args(1)
                .default_value("300")
                .value_parser(value_parser!(usize))
                .help("Minimum length of repetitive fragments"),
        )
        .arg(
            Arg::new("ff")
                .long("ff")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Fill holes between repetitive fragments"),
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

    let opt_kmer = *args.get_one::<usize>("kmer").unwrap();
    let opt_fk = *args.get_one::<usize>("fk").unwrap();
    let opt_min = *args.get_one::<usize>("min").unwrap();
    let opt_ff = *args.get_one::<usize>("ff").unwrap();

    //----------------------------
    // Paths
    //----------------------------
    let curdir = env::current_dir()?;
    let pgr = env::current_exe()?.display().to_string();
    let tempdir = tempfile::Builder::new().prefix("pgr_rm_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"pgr\"     = ${pgr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_repeat = intspan::absolute_path(args.get_one::<String>("repeat").unwrap())?
        .display()
        .to_string();
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
    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> FastK on repeat")?;
    run_cmd!(
        FastK -t -k${opt_kmer} -Nrepeat ${abs_repeat}
    )?;

    run_cmd!(info "==> FastK on genome")?;
    run_cmd!(
        FastK -p:repeat -k${opt_kmer} -Ngenome ${abs_infile}
    )?;

    run_cmd!(info "==> Process each chromosome")?;
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

    let re_prof: regex::Regex = regex::Regex::new(
        r"(?xi)
            (?<start>\d+)       # start
            \s*-\s*             # spacer
            (?<end>\d+)         # end
            ",
    )?;

    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        let sn = i + 1;
        run_cmd!(
            Profex -z genome ${sn} > prof.${sn}.txt
        )?;

        let reader = intspan::reader(&format!("prof.{}.txt", sn));

        let rg_file = format!("prof.{}.rg", sn);
        let mut writer = intspan::writer(&rg_file);

        for line in reader.lines().map_while(Result::ok) {
            let Some(caps) = re_prof.captures(&line) else {
                continue;
            };

            let start = caps["start"].parse::<usize>()? + 1;
            let end = caps["end"].parse::<usize>()? + 1;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }
        rg_files.push(rg_file);
    }

    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        spanr cover $[rg_files] |
            spanr span --op fill -n ${opt_fk} stdin |
            spanr span --op excise -n ${opt_min} stdin |
            spanr span --op fill -n ${opt_ff} stdin -o ${abs_outfile}
    )?;

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}
