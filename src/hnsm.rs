extern crate clap;
use clap::*;

mod cmd;

fn main() -> anyhow::Result<()> {
    let app = Command::new("hnsm")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Homogeneous Nucleic acid Smart Matching")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        .subcommand(cmd::count::make_subcommand())
        .subcommand(cmd::filter::make_subcommand())
        .subcommand(cmd::masked::make_subcommand())
        .subcommand(cmd::n50::make_subcommand())
        .subcommand(cmd::one::make_subcommand())
        .subcommand(cmd::order::make_subcommand())
        .subcommand(cmd::rc::make_subcommand())
        .subcommand(cmd::replace::make_subcommand())
        .subcommand(cmd::sixframe::make_subcommand())
        .subcommand(cmd::size::make_subcommand())
        .subcommand(cmd::some::make_subcommand())
        .subcommand(cmd::split::make_subcommand())
        .after_help(
            r###"
* <infiles> are paths to fasta files, .fa.gz is supported
    * infile == stdin means reading from STDIN

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        // fasta
        Some(("count", sub_matches)) => cmd::count::execute(sub_matches),
        Some(("filter", sub_matches)) => cmd::filter::execute(sub_matches),
        Some(("masked", sub_matches)) => cmd::masked::execute(sub_matches),
        Some(("n50", sub_matches)) => cmd::n50::execute(sub_matches),
        Some(("one", sub_matches)) => cmd::one::execute(sub_matches),
        Some(("order", sub_matches)) => cmd::order::execute(sub_matches),
        Some(("rc", sub_matches)) => cmd::rc::execute(sub_matches),
        Some(("replace", sub_matches)) => cmd::replace::execute(sub_matches),
        Some(("sixframe", sub_matches)) => cmd::sixframe::execute(sub_matches),
        Some(("size", sub_matches)) => cmd::size::execute(sub_matches),
        Some(("some", sub_matches)) => cmd::some::execute(sub_matches),
        Some(("split", sub_matches)) => cmd::split::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}
