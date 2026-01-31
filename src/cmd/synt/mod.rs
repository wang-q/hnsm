use clap::Command;

pub mod dag;
pub mod das;
pub mod dna;
pub mod merge;
pub mod ribbon;
pub mod circle;

pub fn make_subcommand() -> Command {
    Command::new("synt")
        .about("Synteny analysis commands")
        .subcommand_required(true)
        .subcommand(dag::make_subcommand())
        .subcommand(das::make_subcommand())
        .subcommand(dna::make_subcommand())
        .subcommand(merge::make_subcommand())
        .subcommand(ribbon::make_subcommand())
        .subcommand(circle::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("dag", sub_matches)) => dag::execute(sub_matches),
        Some(("das", sub_matches)) => das::execute(sub_matches),
        Some(("dna", sub_matches)) => dna::execute(sub_matches),
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        Some(("ribbon", sub_matches)) => ribbon::execute(sub_matches),
        Some(("circle", sub_matches)) => circle::execute(sub_matches),
        _ => unreachable!(),
    }
}
