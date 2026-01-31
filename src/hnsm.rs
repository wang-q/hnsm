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
        .subcommand(cmd::size::make_subcommand())
        .subcommand(cmd::count::make_subcommand())
        .subcommand(cmd::masked::make_subcommand())
        .subcommand(cmd::n50::make_subcommand())
        // records
        .subcommand(cmd::one::make_subcommand())
        .subcommand(cmd::some::make_subcommand())
        .subcommand(cmd::order::make_subcommand())
        .subcommand(cmd::split::make_subcommand())
        // transform
        .subcommand(cmd::replace::make_subcommand())
        .subcommand(cmd::rc::make_subcommand())
        .subcommand(cmd::filter::make_subcommand())
        .subcommand(cmd::dedup::make_subcommand())
        .subcommand(cmd::mask::make_subcommand())
        .subcommand(cmd::sixframe::make_subcommand())
        // index
        .subcommand(cmd::gz::make_subcommand())
        .subcommand(cmd::range::make_subcommand())
        .subcommand(cmd::prefilter::make_subcommand())
        // fastq
        .subcommand(cmd::interleave::make_subcommand())
        .subcommand(cmd::fq2fa::make_subcommand())
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

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter / dedup / mask / sixframe
    * indexing: gz / range / prefilter
* Fastq files: interleave / fq2fa
* Distance
    * DNA/protein: dist hv / dist seq
    * vectors: dist vector
    * manifold
* Clustering
    * clust cc / clust dbscan / clust km / clust mcl
* Matrix
    * mat compare / mat format / mat pair / mat phylip / mat subset
* Synteny
    * synt chain / synt das / synt dna / synt merge / synt view

"###,
        );

    // Check which subcommand the user ran...
    match app.get_matches().subcommand() {
        Some(("synt", sub_matches)) => cmd::synt::execute(sub_matches),
        Some(("gff", sub_matches)) => cmd::gff::execute(sub_matches),
        Some(("mat", sub_matches)) => cmd::mat::execute(sub_matches),
        Some(("clust", sub_matches)) => cmd::clust::execute(sub_matches),
        Some(("count", sub_matches)) => cmd::count::execute(sub_matches),
        Some(("dedup", sub_matches)) => cmd::dedup::execute(sub_matches),
        Some(("n50", sub_matches)) => cmd::n50::execute(sub_matches),
        Some(("one", sub_matches)) => cmd::one::execute(sub_matches),
        Some(("some", sub_matches)) => cmd::some::execute(sub_matches),
        Some(("order", sub_matches)) => cmd::order::execute(sub_matches),
        Some(("split", sub_matches)) => cmd::split::execute(sub_matches),
        Some(("replace", sub_matches)) => cmd::replace::execute(sub_matches),
        Some(("rc", sub_matches)) => cmd::rc::execute(sub_matches),
        Some(("filter", sub_matches)) => cmd::filter::execute(sub_matches),
        Some(("manifold", sub_matches)) => cmd::manifold::execute(sub_matches),
        Some(("mask", sub_matches)) => cmd::mask::execute(sub_matches),
        Some(("masked", sub_matches)) => cmd::masked::execute(sub_matches),
        Some(("sixframe", sub_matches)) => cmd::sixframe::execute(sub_matches),
        Some(("size", sub_matches)) => cmd::size::execute(sub_matches),
        Some(("gz", sub_matches)) => cmd::gz::execute(sub_matches),
        Some(("range", sub_matches)) => cmd::range::execute(sub_matches),
        Some(("prefilter", sub_matches)) => cmd::prefilter::execute(sub_matches),
        Some(("interleave", sub_matches)) => cmd::interleave::execute(sub_matches),
        Some(("fq2fa", sub_matches)) => cmd::fq2fa::execute(sub_matches),
        Some(("dist", sub_matches)) => cmd::dist::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: Remove fully contained sequences
