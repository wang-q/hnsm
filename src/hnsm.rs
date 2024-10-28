#![feature(array_chunks)]
#![feature(slice_as_chunks)]
// Add these imports to use the stdsimd library
#![feature(portable_simd)]

extern crate clap;
use clap::*;

mod cmd;

fn main() -> anyhow::Result<()> {
    let app = Command::new("hnsm")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Homogeneous Nucleic acids/amino acids Smart Matching")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        .subcommand(cmd::chain::make_subcommand())
        .subcommand(cmd::cluster::make_subcommand())
        .subcommand(cmd::count::make_subcommand())
        .subcommand(cmd::das::make_subcommand())
        .subcommand(cmd::distance::make_subcommand())
        .subcommand(cmd::filter::make_subcommand())
        .subcommand(cmd::gz::make_subcommand())
        .subcommand(cmd::masked::make_subcommand())
        .subcommand(cmd::n50::make_subcommand())
        .subcommand(cmd::one::make_subcommand())
        .subcommand(cmd::order::make_subcommand())
        .subcommand(cmd::range::make_subcommand())
        .subcommand(cmd::rc::make_subcommand())
        .subcommand(cmd::replace::make_subcommand())
        .subcommand(cmd::manifold::make_subcommand())
        .subcommand(cmd::similarity::make_subcommand())
        .subcommand(cmd::sixframe::make_subcommand())
        .subcommand(cmd::size::make_subcommand())
        .subcommand(cmd::some::make_subcommand())
        .subcommand(cmd::split::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter
    * indexing: gz / range

* Clustering
    * vectors: similarity
    * DNA/protein: distance / identity
    * cluster
    * reduction

* Synteny
    * das
    * chain

* <infiles> are paths to fasta files, .fa.gz is supported
    * infile == stdin means reading from STDIN
    * `hnsm gz` writes out the BGZF format and `hnsm range` reads it

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        // info
        Some(("size", sub_matches)) => cmd::size::execute(sub_matches),
        Some(("count", sub_matches)) => cmd::count::execute(sub_matches),
        Some(("masked", sub_matches)) => cmd::masked::execute(sub_matches),
        Some(("n50", sub_matches)) => cmd::n50::execute(sub_matches),
        // records
        Some(("one", sub_matches)) => cmd::one::execute(sub_matches),
        Some(("some", sub_matches)) => cmd::some::execute(sub_matches),
        Some(("order", sub_matches)) => cmd::order::execute(sub_matches),
        Some(("split", sub_matches)) => cmd::split::execute(sub_matches),
        // transform
        Some(("replace", sub_matches)) => cmd::replace::execute(sub_matches),
        Some(("rc", sub_matches)) => cmd::rc::execute(sub_matches),
        Some(("filter", sub_matches)) => cmd::filter::execute(sub_matches),
        //
        Some(("sixframe", sub_matches)) => cmd::sixframe::execute(sub_matches),
        // index
        Some(("gz", sub_matches)) => cmd::gz::execute(sub_matches),
        Some(("range", sub_matches)) => cmd::range::execute(sub_matches),
        // clustering
        Some(("distance", sub_matches)) => cmd::distance::execute(sub_matches),
        Some(("similarity", sub_matches)) => cmd::similarity::execute(sub_matches),
        Some(("cluster", sub_matches)) => cmd::cluster::execute(sub_matches),
        Some(("manifold", sub_matches)) => cmd::manifold::execute(sub_matches),
        // Synteny
        Some(("das", sub_matches)) => cmd::das::execute(sub_matches),
        Some(("chain", sub_matches)) => cmd::chain::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}

// TODO:
//  interleave
//  sort
//  identity: accurate pairwise sequence identity
//    https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
