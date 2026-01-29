use clap::Command;

pub mod das;
pub mod dna;
pub mod merge;
pub mod view;

pub fn make_subcommand() -> Command {
    Command::new("synt")
        .about("Synteny analysis commands")
        .subcommand_required(true)
        .subcommand(das::make_subcommand())
        .subcommand(dna::make_subcommand())
        .subcommand(merge::make_subcommand())
        .subcommand(view::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("das", sub_matches)) => das::execute(sub_matches),
        Some(("dna", sub_matches)) => dna::execute(sub_matches),
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        Some(("view", sub_matches)) => view::execute(sub_matches),
        _ => unreachable!(),
    }
}
