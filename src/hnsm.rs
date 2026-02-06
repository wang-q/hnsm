extern crate clap;
use clap::*;

mod cmd;

fn main() -> anyhow::Result<()> {
    let app = Command::new("hnsm")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Homogeneous Nucleic/amino acids Smart Matching")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        // info
        // records
        // transform
        // index
        // distance
        .subcommand(cmd::manifold::make_subcommand())
        // synteny
        .subcommand(cmd::synt::make_subcommand())
        // gff
        .subcommand(cmd::gff::make_subcommand())
        .after_help(
            r###"Subcommand groups:
* Distance
    * manifold
* Synteny
    * synt dag / synt das / synt mmg / synt merge / synt ribbon / synt circle
* GFF
    * gff rg
"###,
        );

    // Check which subcommand the user ran...
    match app.get_matches().subcommand() {
        Some(("synt", sub_matches)) => cmd::synt::execute(sub_matches),
        Some(("gff", sub_matches)) => cmd::gff::execute(sub_matches),
        Some(("manifold", sub_matches)) => cmd::manifold::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

