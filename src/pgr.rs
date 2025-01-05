extern crate clap;
use clap::*;

mod cmd_pgr;

fn main() -> anyhow::Result<()> {
    let app = Command::new("pgr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("`pgr` - Population Genomes Refiner")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        .subcommand(cmd_pgr::chain::make_subcommand())
        .subcommand(cmd_pgr::ir::make_subcommand())
        .subcommand(cmd_pgr::rept::make_subcommand())
        .subcommand(cmd_pgr::trf::make_subcommand())
        .after_help(
            r###"
`pgr` can also refer to `Populations defined by Gene transfeR`,
as a reimplementation of PopCOGenT (populations as clusters of gene transfer)

Subcommand groups:

* Genome alignments:
    * chain

* Repeats:
    * ir / rept / trf

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("chain", sub_matches)) => cmd_pgr::chain::execute(sub_matches),
        Some(("ir", sub_matches)) => cmd_pgr::ir::execute(sub_matches),
        Some(("rept", sub_matches)) => cmd_pgr::rept::execute(sub_matches),
        Some(("trf", sub_matches)) => cmd_pgr::trf::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: paralog
