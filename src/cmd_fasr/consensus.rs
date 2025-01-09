use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("consensus")
        .about("Generate consensus sequences using POA")
        .after_help(
            r###"
This subcommand generates consensus sequences from block FA files using the POA (Partial Order Alignment) algorithm.

Input files can be gzipped. If the input file is 'stdin', data is read from standard input.

Note:
- Requires `spoa` to be installed and available in $PATH.
    * The original `poa` was unstable and sometimes crashed
- Supports parallel processing for improved performance.
    * Running in parallel mode with 1 reader, 1 writer and the corresponding number of workers
    * The order of output may be different from the original
- If outgroups are present, they are handled appropriately.

Examples:
1. Generate consensus sequences from a block FA file:
   fasr consensus tests/fasr/example.fas

2. Generate consensus sequences with outgroups:
   fasr consensus tests/fasr/example.fas --outgroup

3. Run in parallel with 4 threads:
   fasr consensus tests/fasr/example.fas --parallel 4

4. Output results to a file:
   fasr consensus tests/fasr/example.fas -o output.fas


"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
        .arg(
            Arg::new("cname")
                .long("cname")
                .num_args(1)
                .default_value("consensus")
                .help("Name of the consensus"),
        )
        .arg(
            Arg::new("has_outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("Indicates the presence of outgroups at the end of each block"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1")
                .help("Number of threads for parallel processing"),
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
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    //----------------------------
    // Operating
    //----------------------------
    if opt_parallel == 1 {
        // Single-threaded mode
        let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

        for infile in args.get_many::<String>("infiles").unwrap() {
            let mut reader = intspan::reader(infile);
            while let Ok(block) = hnsm::next_fas_block(&mut reader) {
                let out_string = proc_block(&block, args)?;
                writer.write_all(out_string.as_ref())?;
            }
        }
    } else {
        // Parallel mode
        proc_block_p(args)?;
    }

    Ok(())
}

fn proc_block(block: &hnsm::FasBlock, args: &ArgMatches) -> anyhow::Result<String> {
    //----------------------------
    // Args
    //----------------------------
    let cname = args.get_one::<String>("cname").unwrap();
    let has_outgroup = args.get_flag("has_outgroup");

    //----------------------------
    // Ops
    //----------------------------
    let mut seqs = vec![];

    let outgroup = if has_outgroup {
        Some(block.entries.iter().last().unwrap())
    } else {
        None
    };

    for entry in &block.entries {
        seqs.push(entry.seq().as_ref());
    }
    if outgroup.is_some() {
        seqs.pop().unwrap(); // Remove the outgroup sequence
    }

    // Generate consensus sequence
    let mut cons = hnsm::get_consensus_poa(&seqs).unwrap();
    cons = cons.replace('-', "");

    let mut range = block.entries.first().unwrap().range().clone();

    //----------------------------
    // Output
    //----------------------------
    let mut out_string = "".to_string();
    if range.is_valid() {
        *range.name_mut() = cname.to_string();
        out_string += format!(">{}\n{}\n", range, cons).as_ref();
    } else {
        out_string += format!(">{}\n{}\n", cname, cons).as_ref();
    }
    if outgroup.is_some() {
        out_string += outgroup.unwrap().to_string().as_ref();
    }

    // end of a block
    out_string += "\n";

    Ok(out_string)
}

// Adopt from https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html#create-a-parallel-pipeline
fn proc_block_p(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    // Channel 1 - Read files to blocks
    let (snd1, rcv1) = crossbeam::channel::bounded::<hnsm::FasBlock>(10);
    // Channel 2 - Results
    let (snd2, rcv2) = crossbeam::channel::bounded(10);

    crossbeam::scope(|s| {
        //----------------------------
        // Reader thread
        //----------------------------
        s.spawn(|_| {
            for infile in args.get_many::<String>("infiles").unwrap() {
                let mut reader = intspan::reader(infile);
                while let Ok(block) = hnsm::next_fas_block(&mut reader) {
                    snd1.send(block).unwrap();
                }
            }
            // Close the channel - this is necessary to exit the for-loop in the worker
            drop(snd1);
        });

        //----------------------------
        // Worker threads
        //----------------------------
        for _ in 0..parallel {
            // Send to sink, receive from source
            let (sendr, recvr) = (snd2.clone(), rcv1.clone());
            // Spawn workers in separate threads
            s.spawn(move |_| {
                // Receive until channel closes
                for block in recvr.iter() {
                    let out_string = proc_block(&block, args).unwrap();
                    sendr.send(out_string).unwrap();
                }
            });
        }
        // Close the channel, otherwise sink will never exit the for-loop
        drop(snd2);

        //----------------------------
        // Writer thread
        //----------------------------
        for out_string in rcv2.iter() {
            writer.write_all(out_string.as_ref()).unwrap();
        }
    })
    .unwrap();

    Ok(())
}
