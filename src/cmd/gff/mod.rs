use clap::Command;

pub mod rg;

pub fn make_subcommand() -> Command {
    Command::new("gff")
        .about("GFF file operations")
        .subcommand_required(true)
        .subcommand(rg::make_subcommand())
}

pub fn execute(matches: &clap::ArgMatches) -> anyhow::Result<()> {
    match matches.subcommand() {
        Some(("rg", sub_matches)) => rg::execute(sub_matches),
        _ => unreachable!(),
    }
}
