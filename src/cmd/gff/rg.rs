use clap::*;
use noodles_gff as gff;
use std::io::Write;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("rg")
        .about("Extract ranges from GFF files")
         .after_help(
            r###"
Extract ranges from GFF files.

Examples:
1. Extract default 'gene' features:
   hnsm gff rg tests/gff_rg/test.gff

2. Extract 'mRNA' features with assembly name 'Human':
   hnsm gff rg tests/gff_rg/test.gff --tag mRNA --asm Human -o output.tsv

3. Simplify sequence IDs:
   hnsm gff rg tests/gff_rg/test.gff --simplify

"###,
        )
       .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("tag")
                .long("tag")
                .num_args(1)
                .default_value("gene")
                .help("Features to retain"),
        )
        .arg(
            Arg::new("asm")
                .long("asm")
                .num_args(1)
                .help("Assembly name. Default: inferred from filename"),
        )
        .arg(
            Arg::new("simplify")
                .long("simplify")
                .short('s')
                .action(ArgAction::SetTrue)
                .help("Simplify sequence names"),
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
    let outfile = args.get_one::<String>("outfile").unwrap();
    let tag = args.get_one::<String>("tag").unwrap();
    let is_simplify = args.get_flag("simplify");

    let asm = if let Some(g) = args.get_one::<String>("asm") {
        g.clone()
    } else {
        let path = std::path::Path::new(infile);
        let file_name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        file_name.split('.').next().unwrap_or("unknown").to_string()
    };

    eprintln!("gff rg: infile={}, outfile={}, tag={}, asm={}, simplify={}", infile, outfile, tag, asm, is_simplify);

    let reader = intspan::reader(infile);
    let mut reader = gff::io::Reader::new(reader);
    let mut writer = intspan::writer(outfile);

    for result in reader.record_bufs() {
        let record = result?;
        if record.ty().to_ascii_lowercase() != tag.to_ascii_lowercase().as_bytes() {
            continue;
        }

        // ID
        let mut id = match record.attributes().get(b"ID") {
            Some(gff::feature::record_buf::attributes::field::Value::String(s)) => s.to_string(),
            _ => "NA".to_string(),
        };
        if is_simplify {
            if let Some((_, s)) = id.split_once(':') {
                id = s.to_string();
            }
            if let Some(i) = id.find(&[' ', '.', ',', '-'][..]) {
                id = id[..i].to_string();
            }
        }
        
        // Range
        let seq_name = record.reference_sequence_name();
        let strand = match record.strand() {
            gff::feature::record::Strand::Reverse => "-",
            _ => "+",
        };
        let start = record.start();
        let end = record.end();

        writeln!(
            writer,
            "{}\t{}.{}({}):{}-{}",
            id, asm, seq_name, strand, start, end
        )?;
    }

    Ok(())
}
