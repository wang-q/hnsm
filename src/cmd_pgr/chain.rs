use clap::*;
use cmd_lib::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("chain")
        .about("UCSC chain/net pipeline")
        .after_help(
            r###"
This command implements the UCSC pipeline for pairwise genome alignments, psl-chain-net-axt-maf.

* <target> and <query> are fasta files
* <psl> can be a .psl file or a directory containing multiple .psl files
* Default names of target and query in the output .maf are derived from the basename of <target> and <query>

* `--lineargap` and `--minscore`:
    * Human18vsChimp2 use `loose` and 1000
    * Human19vsChimp3 use `medium` and 5000
    * `loose` corresponds to chicken/human linear gap costs
    * `medium` corresponds to mouse/human linear gap costs

* The following binaries from the kent-tools are required and should be found in $PATH:
    * axtChain
    * chainAntiRepeat
    * chainMergeSort
    * chainPreNet
    * chainNet
    * netSyntenic
    * netChainSubset
    * chainStitchId
    * netSplit
    * netToAxt
    * axtSort
    * axtToMaf
    * netFilter
    * chainSplit

Definitions:

* The *target* is the reference genome sequence
* The *query* is some other genome sequence

* A *chain* is a sequence of non-overlapping gapless blocks, with single- or double-sided gaps between blocks.
  Within a chain, target and query coords are monotonically non-decreasing.
* A *net* is a hierarchical collection of chains.

References:

* [Chains Nets](https://genomewiki.ucsc.edu/index.php/Chains_Nets)
* [Prebuild binaries](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)

"###,
        )
        .arg(
            Arg::new("target")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Path to the target genome FA file"),
        )
        .arg(
            Arg::new("query")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Path to the query genome FA file"),
        )
        .arg(
            Arg::new("psl")
                .required(true)
                .num_args(1)
                .index(3)
                .help("Path to the PSL file or directory containing PSL files"),
        )
        .arg(
            Arg::new("lineargap")
                .long("lineargap")
                .action(ArgAction::Set)
                .value_parser([
                    builder::PossibleValue::new("loose"),
                    builder::PossibleValue::new("medium"),
                ])
                .default_value("loose")
                .help("Linear gap cost setting for axtChain"),
        )
        .arg(
            Arg::new("minscore")
                .long("minscore")
                .num_args(1)
                .default_value("1000")
                .value_parser(value_parser!(usize))
                .help("Minimum alignment score for axtChain"),
        )
        .arg(
            Arg::new("tname")
                .long("tname")
                .short('t')
                .num_args(1)
                .help("Custom name for the target genome"),
        )
        .arg(
            Arg::new("qname")
                .long("qname")
                .short('q')
                .num_args(1)
                .help("Custom name for the query genome"),
        )
        .arg(
            Arg::new("syn")
                .long("syn")
                .action(ArgAction::SetTrue)
                .help("Generate syntenic alignments"),
        )
        .arg(
            Arg::new("outdir")
                .short('o')
                .long("outdir")
                .num_args(1)
                .default_value("stdout")
                .help("Output location. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_lineargap = args.get_one::<String>("lineargap").unwrap();
    let opt_minscore = *args.get_one::<usize>("minscore").unwrap();

    let is_syn = args.get_flag("syn");

    //----------------------------
    // Paths
    //----------------------------
    let curdir = std::env::current_dir()?;
    let pgr = std::env::current_exe()?.display().to_string();
    let tempdir = tempfile::Builder::new().prefix("pgr_chain_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"pgr\"     = ${pgr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    run_cmd!(info "==> Absolute paths")?;
    let abs_target = intspan::absolute_path(args.get_one::<String>("target").unwrap())?
        .display()
        .to_string();
    let abs_query = intspan::absolute_path(args.get_one::<String>("query").unwrap())?
        .display()
        .to_string();

    let opt_tname = if let Some(tname) = args.get_one::<String>("tname") {
        if tname.is_empty() {
            "".to_string()
        } else {
            format!("{}.", tname)
        }
    } else {
        format!("{}.", get_basename(&abs_target).unwrap())
    };
    let opt_qname = if let Some(qname) = args.get_one::<String>("qname") {
        if qname.is_empty() {
            "".to_string()
        } else {
            format!("{}.", qname)
        }
    } else {
        format!("{}.", get_basename(&abs_query).unwrap())
    };

    let abs_psl = intspan::absolute_path(args.get_one::<String>("psl").unwrap())?
        .display()
        .to_string();
    let infiles = if std::path::Path::new(&abs_psl).is_dir() {
        list_files_ext(&abs_psl, "psl")?
    } else {
        vec![abs_psl]
    };

    let abs_outdir = if outdir == "stdout" {
        outdir.to_string()
    } else {
        intspan::absolute_path(outdir)?.display().to_string()
    };

    //----------------------------
    // Ops
    //----------------------------
    run_cmd!(info "==> Switch to tempdir")?;
    std::env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Target .sizes anc .2bit")?;
    run_cmd!(
        hnsm size ${abs_target} -o target.chr.sizes;
        faToTwoBit ${abs_target} target.chr.2bit;
    )?;
    run_cmd!(info "==> Query .sizes anc .2bit")?;
    run_cmd!(
        hnsm size ${abs_query} -o query.chr.sizes;
        faToTwoBit ${abs_query} query.chr.2bit;
    )?;

    run_cmd!(info "==> axtChain")?;
    // axtChain - Chain together axt alignments.
    // usage:
    //   axtChain -linearGap=loose in.axt tNibDir qNibDir out.chain
    // Where tNibDir/qNibDir are either directories full of nib files, or the
    // name of a .2bit file
    //
    // chainAntiRepeat - Get rid of chains that are primarily the results of
    // repeats and degenerate DNA
    // usage:
    //    chainAntiRepeat tNibDir qNibDir inChain outChain
    // options:
    //    -minScore=N - minimum score (after repeat stuff) to pass
    //    -noCheckScore=N - score that will pass without checks (speed tweak)
    std::fs::create_dir_all("pslChain")?;
    for infile in infiles {
        let stem = get_basename(&infile).unwrap();
        run_cmd!(
            axtChain -minScore=${opt_minscore} -linearGap=${opt_lineargap} -psl ${infile} target.chr.2bit query.chr.2bit stdout |
                chainAntiRepeat target.chr.2bit query.chr.2bit stdin pslChain/${stem}.chain
        )?;
    }

    run_cmd!(info "==> chainMergeSort and chainPreNet")?;
    {
        // This step would open all .chain files and reach system's maxfile limit.
        // So merge 100 files a time.
        //
        // chainMergeSort - Combine sorted files into larger sorted file
        // usage:
        //    chainMergeSort file(s)
        // Output goes to standard output
        // options:
        //    -saveId - keep the existing chain ids.
        //    -inputList=somefile - somefile contains list of input chain files.
        //    -tempDir=somedir/ - somedir has space for temporary sorting data, default ./
        let mut files = list_files_ext("pslChain", "chain")?;
        let mut sn = 1;
        let mut merge_files = vec![];
        while !files.is_empty() {
            let batching: Vec<_> = files.drain(0..100.min(files.len())).collect();

            intspan::write_lines(
                "chainList.tmp",
                &batching.iter().map(AsRef::as_ref).collect(),
            )?;
            run_cmd!(
                chainMergeSort -inputList=chainList.tmp > all.${sn}.chain.tmp
            )?;
            merge_files.push(format!("all.{}.chain.tmp", sn));

            sn += 1;
        }

        run_cmd!(
            chainMergeSort $[merge_files] > all.chain
        )?;

        // chainPreNet - Remove chains that don't have a chance of being netted
        // usage:
        //   chainPreNet in.chain target.sizes query.sizes out.chain
        run_cmd!(
            chainPreNet all.chain target.chr.sizes query.chr.sizes all.pre.chain
        )?;
    }

    run_cmd!(info "==> chain-net")?;
    {
        // chainNet - Make alignment nets out of chains
        // usage:
        //   chainNet in.chain target.sizes query.sizes target.net query.net
        //
        // netSyntenic - Add synteny info to net.
        // usage:
        //   netSyntenic in.net out.net
        run_cmd!(
            chainNet -minSpace=1 all.pre.chain target.chr.sizes query.chr.sizes stdout query.chainnet |
                netSyntenic stdin noClass.net
        )?;

        // netChainSubset - Create chain file with subset of chains that appear in
        // the net
        // usage:
        //    netChainSubset in.net in.chain out.chain
        // options:
        //    -gapOut=gap.tab - Output gap sizes to file
        //    -type=XXX - Restrict output to particular type in net file
        //    -splitOnInsert - Split chain when get an insertion of another chain
        //    -wholeChains - Write entire chain references by net, don't split
        //     when a high-level net is encoundered.  This is useful when nets
        //     have been filtered.
        //    -skipMissing - skip chains that are not found instead of generating
        //     an error.  Useful if chains have been filtered.
        //
        // chainStitchId - Join chain fragments with the same chain ID into a single
        //    chain per ID.  Chain fragments must be from same original chain but
        //    must not overlap.  Chain fragment scores are summed.
        // usage:
        //    chainStitchId in.chain out.chain
        run_cmd!(
            netChainSubset -verbose=0 noClass.net all.chain stdout |
                chainStitchId stdin over.chain
        )?;

        // netSplit - Split a genome net file into chromosome net files
        // usage:
        //   netSplit in.net outDir
        std::fs::create_dir_all("net")?;
        run_cmd!(
            netSplit noClass.net net > /dev/null
        )?;
    }

    run_cmd!(info "==> netToAxt")?;
    {
        std::fs::create_dir_all("axtNet")?;

        let files = list_files_ext("net", "net")?;

        // netToAxt - Convert net (and chain) to axt.
        // usage:
        //   netToAxt in.net in.chain target.2bit query.2bit out.axt
        // note:
        // directories full of .nib files (an older format)
        // may also be used in place of target.2bit and query.2bit.
        //
        // axtSort - Sort axt files
        // usage:
        //   axtSort in.axt out.axt
        for file in files {
            let stem = get_basename(&file).unwrap();
            run_cmd!(
                netToAxt ${file} all.pre.chain target.chr.2bit query.chr.2bit stdout |
                    axtSort stdin axtNet/${stem}.axt
            )?;
        }
    }

    run_cmd!(info "==> axt-maf")?;
    if !is_syn {
        run_cmd!(info "==> axtToMaf")?;

        let files = list_files_ext("axtNet", "axt")?;
        for file in files {
            let stem = get_basename(&file).unwrap();
            if abs_outdir == "stdout" {
                if opt_tname.is_empty() {
                    run_cmd!(
                        axtToMaf ${file} target.chr.sizes query.chr.sizes stdout
                    )?;
                } else {
                    run_cmd!(
                        axtToMaf -tPrefix=${opt_tname} -qPrefix=${opt_qname} ${file} target.chr.sizes query.chr.sizes stdout
                    )?;
                }
            } else {
                if opt_tname.is_empty() {
                    run_cmd!(
                        axtToMaf ${file} target.chr.sizes query.chr.sizes ${abs_outdir}/${stem}.maf
                    )?;
                } else {
                    run_cmd!(
                        axtToMaf -tPrefix=${opt_tname} -qPrefix=${opt_qname} ${file} target.chr.sizes query.chr.sizes ${abs_outdir}/${stem}.maf
                    )?;
                }
            }
        }
    } else {
        std::fs::create_dir_all("synNet")?;
        std::fs::create_dir_all("chain")?;

        run_cmd!(info "==> synNet.maf")?;

        // netFilter - Filter out parts of net.  What passes
        // filter goes to standard output.  Note a net is a
        // recursive data structure.  If a parent fails to pass
        // the filter, the children are not even considered.
        // usage:
        //    netFilter in.net(s)
        run_cmd!(
            netFilter -syn noClass.net |
                netSplit stdin synNet > /dev/null
        )?;

        // chainSplit - Split chains up by target or query sequence
        // usage:
        //    chainSplit outDir inChain(s)
        // options:
        //    -q  - Split on query (default is on target)
        //    -lump=N  Lump together so have only N split files.
        run_cmd!(
            chainSplit synNet all.chain
        )?;

        let files = list_files_ext("synNet", "net")?;
        for file in files {
            let stem = get_basename(&file).unwrap();
            let chain_file = format!("{}.chain", file.strip_suffix(".net").unwrap());
            if abs_outdir == "stdout" {
                if opt_tname.is_empty() {
                    run_cmd!(
                        netToAxt ${file} ${chain_file} target.chr.2bit query.chr.2bit stdout |
                            axtSort stdin stdout |
                            axtToMaf stdin target.chr.sizes query.chr.sizes stdout
                    )?;
                } else {
                    run_cmd!(
                        netToAxt ${file} ${chain_file} target.chr.2bit query.chr.2bit stdout |
                            axtSort stdin stdout |
                            axtToMaf -tPrefix=${opt_tname} -qPrefix=${opt_qname} stdin target.chr.sizes query.chr.sizes stdout
                    )?;
                }
            } else {
                if opt_tname.is_empty() {
                    run_cmd!(
                        netToAxt ${file} ${chain_file} target.chr.2bit query.chr.2bit stdout |
                            axtSort stdin stdout |
                            axtToMaf stdin target.chr.sizes query.chr.sizes ${abs_outdir}/${stem}.maf
                    )?;
                } else {
                    run_cmd!(
                        netToAxt ${file} ${chain_file} target.chr.2bit query.chr.2bit stdout |
                            axtSort stdin stdout |
                            axtToMaf -tPrefix=${opt_tname} -qPrefix=${opt_qname} stdin target.chr.sizes query.chr.sizes ${abs_outdir}/${stem}.maf
                    )?;
                }
            }
        }
    }

    //----------------------------
    // Done
    //----------------------------
    std::env::set_current_dir(&curdir)?;

    Ok(())
}

fn list_files_ext(dir: &str, extension: &str) -> Result<Vec<String>, std::io::Error> {
    let mut files = Vec::new();
    let dir_path = std::path::Path::new(dir);

    if dir_path.is_dir() {
        for entry in std::fs::read_dir(dir_path)? {
            let entry = entry?;
            let path = entry.path();

            if path.is_file() {
                if let Some(ext) = path.extension() {
                    if ext == extension {
                        files.push(path.to_string_lossy().into_owned());
                    }
                }
            }
        }
    }

    Ok(files)
}

fn get_basename(file_path: &str) -> Option<String> {
    let path = std::path::Path::new(file_path);

    if let Some(file_name) = path.file_stem().and_then(std::ffi::OsStr::to_str) {
        return Some(file_name.split('.').next().unwrap().to_string());
    }

    None
}
