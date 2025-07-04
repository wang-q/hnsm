# Change Log

## Unreleased - ReleaseDate

## 0.3.4 - 2025-06-05

* Enhance `hnsm fq2fa` command
  * Add support for processing multiple input files

## 0.3.3 - 2025-06-05

* Add `hnsm fq2fa` command for converting FASTQ to FASTA format
* Refactor clustering commands
    * Merge `cc` and `dbscan` modes into a unified `clust` command
    * Improve component sorting in connected components clustering
    * Support multiple output formats
* Move matrix-related functionality to `intspan` crate
* Update dependencies
    * intspan 0.8.4
    * noodles-bgzf 0.36.0

## 0.3.2 - 2025-04-02

* Improve documentation for all commands
* Update dependencies
    * rapidhash 1.4.0
    * rand 0.9.0
* Improve CI/CD
    * Use cargo-zigbuild for Linux builds
    * Update Rust toolchain to nightly-2025-01-16

## 0.3.1 - 2025-02-08

* Move `pgr` and `fasr` out

## 0.3.0 - 2025-01-20

* Add `hnsm mask`
* Add `hnsm sixframe`
* Add `hnsm hv`
* Add `hnsm prefilter`
* Add `--merge`, `--list`, and `--parallel` to `hnsm distance`

* Improve help texts
* `libs/loc.rs`
    * Use IndexMap
    * Add `read_offset()`
    * Add `records_offset()`
* `libs/hash.rs`
    * Add `seq_mins()`

## 0.2.0 - 2024-12-20

* Add `hnsm convert`
* Add `hnsm interleave`

## 0.1.11 - 2024-11-04

* Add `hnsm dedup`

* Unfinished `hnsm das`
* Unfinished `hnsm chain`
* Unfinished `hnsm manifold`

* Add `--mode cc` to `hnsm cluster`

## 0.1.10 - 2024-10-07

* Add `hnsm similarity`
* Add `hnsm distance`
* Add `hnsm cluster`

## 0.1.8 - 2024-07-19

* LRU cache for `hnsm range`
* Tweaks publish.yml

## 0.1.1 - 2024-07-13

* Index Subcommands
    * gz
    * range

## 0.1.0 - 2024-07-11

* Skeletons, need to be filled

* Fasta Subcommands
    * size
    * one
    * some
    * order
    * masked
    * rc
    * count
    * replace
    * filter
    * split name
    * split about
    * n50
