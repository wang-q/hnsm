use clap::*;

pub mod cc;
pub mod dbscan;

/// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("clust")
        .about("Clustering commands")
        .subcommand_required(true)
        .subcommand(cc::make_subcommand())
        .subcommand(dbscan::make_subcommand())
}

/// Execute pkg command
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("cc", sub_args)) => cc::execute(sub_args),
        Some(("dbscan", sub_args)) => dbscan::execute(sub_args),
        _ => unreachable!("Exhausted list of subcommands and subcommand_required prevents `None`"),
    }
}
