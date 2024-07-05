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
        .subcommand(cmd::sixframe::make_subcommand())
        .subcommand(cmd::size::make_subcommand());

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("sixframe", sub_matches)) => cmd::sixframe::execute(sub_matches),
        Some(("size", sub_matches)) => cmd::size::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}
