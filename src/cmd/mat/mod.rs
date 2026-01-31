use clap::Command;

pub mod compare;
pub mod format;
pub mod pair;
pub mod phylip;
pub mod subset;

pub fn make_subcommand() -> Command {
    Command::new("mat")
        .about("Matrix operations")
        .subcommand_required(true)
        .subcommand(compare::make_subcommand())
        .subcommand(format::make_subcommand())
        .subcommand(pair::make_subcommand())
        .subcommand(phylip::make_subcommand())
        .subcommand(subset::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("compare", sub_matches)) => compare::execute(sub_matches),
        Some(("format", sub_matches)) => format::execute(sub_matches),
        Some(("pair", sub_matches)) => pair::execute(sub_matches),
        Some(("phylip", sub_matches)) => phylip::execute(sub_matches),
        Some(("subset", sub_matches)) => subset::execute(sub_matches),
        _ => unreachable!(),
    }
}
