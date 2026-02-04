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
        .subcommand(cmd::prefilter::make_subcommand())
        // distance
        .subcommand(cmd::dist::make_subcommand())
        .subcommand(cmd::manifold::make_subcommand())
        // clustering
        .subcommand(cmd::clust::make_subcommand())
        // synteny
        .subcommand(cmd::synt::make_subcommand())
        // gff
        .subcommand(cmd::gff::make_subcommand())
        // mat
        .subcommand(cmd::mat::make_subcommand())
        .after_help(
            r###"Subcommand groups:
* Indexing
    * prefilter
* Distance
    * dist hv / dist seq / dist vector
    * manifold
* Clustering
    * clust cc / clust dbscan / clust km / clust mcl
* Matrix
    * mat compare / mat format / mat pair / mat phylip / mat subset
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
        Some(("mat", sub_matches)) => cmd::mat::execute(sub_matches),
        Some(("clust", sub_matches)) => cmd::clust::execute(sub_matches),
        Some(("manifold", sub_matches)) => cmd::manifold::execute(sub_matches),
        Some(("prefilter", sub_matches)) => cmd::prefilter::execute(sub_matches),
        Some(("dist", sub_matches)) => cmd::dist::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: Remove fully contained sequences
