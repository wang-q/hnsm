use clap::Command;

pub mod dna;
pub mod merge;

pub fn make_subcommand() -> Command {
    Command::new("synt")
        .about("Synteny analysis commands")
        .subcommand_required(true)
        .subcommand(dna::make_subcommand())
        .subcommand(merge::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("dna", sub_matches)) => dna::execute(sub_matches),
        Some(("merge", sub_matches)) => merge::execute(sub_matches),
        _ => unreachable!(),
    }
}
