use clap::Command;

pub mod phylip;

pub fn make_subcommand() -> Command {
    Command::new("mat")
        .about("Matrix operations")
        .subcommand_required(true)
        .subcommand(phylip::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("phylip", sub_matches)) => phylip::execute(sub_matches),
        _ => unreachable!(),
    }
}
