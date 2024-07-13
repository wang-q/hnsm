# hnsm

[![Build](https://github.com/wang-q/hnsm/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![codecov](https://codecov.io/gh/wang-q/hnsm/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/hnsm)
[![license](https://img.shields.io/github/license/wang-q/hnsm)](https://github.com//wang-q/hnsm)

`hnsm` - **H**omogeneous **N**ucleic acid **S**mart **M**atching

## Install

Current release: 0.1.2

```shell
cargo install --path . --force --offline

# test
cargo test -- --test-threads=1

# build under WSL 2
export CARGO_TARGET_DIR=/tmp
cargo build

```

## Synopsis

### `hnsm help`

```text
$ hnsm help
Homogeneous Nucleic acid Smart Matching

Usage: hnsm [COMMAND]

Commands:
  count     Count base statistics in FA file(s)
  filter    Filter records in FA file(s)
  gz        Compressing a file using the blocked gzip format (BGZF)
  masked    Masked regions in FA file(s)
  n50       Count total bases in FA file(s)
  one       Extract one FA record
  order     Extract some FA records by the given order
  range     Extract sequences defined by the range(s)
  rc        Reverse complement a FA file
  replace   Replace headers of a FA file
  sixframe  Six-Frame Translation
  size      Count total bases in FA file(s)
  some      Extract some FA records
  split     Split FA file(s) into several files
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version


* <infiles> are paths to fasta files, .fa.gz is supported
    * infile == stdin means reading from STDIN
    * `hnsm gz` writes out the BGZF format and `hnsm range` reads it

```

## Examples

### Fasta files

```shell
hnsm size tests/fasta/ufasta.fa
hnsm count tests/fasta/ufasta.fa.gz
hnsm masked tests/fasta/ufasta.fa
hnsm n50 tests/fasta/ufasta.fa -N 90 -N 50 -S -t

hnsm one tests/fasta/ufasta.fa read12
hnsm some tests/fasta/ufasta.fa tests/fasta/list.txt
hnsm order tests/fasta/ufasta.fa tests/fasta/list.txt

hnsm replace tests/fasta/ufasta.fa tests/fasta/replace.tsv
hnsm rc tests/fasta/ufasta.fa

hnsm filter -a 10 -z 50 -U tests/fasta/ufasta.fa
hnsm filter -a 1 -u tests/fasta/ufasta.fa tests/fasta/ufasta.fa.gz
hnsm filter --iupac --upper tests/fasta/filter.fa

hnsm filter -a 400 tests/fasta/ufasta.fa |
    hnsm split name stdin -o tmp
hnsm split about -c 2000 tests/fasta/ufasta.fa -o tmp

cargo run --bin hnsm sixframe

cargo run --bin hnsm sort

```

### Index

`samtools faidx` is designed for fast randomized extraction of sequences from reference sequences,
and requires that the sequence file be "well-formatted", i.e., all sequence lines must be the same
length, which is to facilitate random access to disk files. For a mammal reference genome, this
requirement is reasonable; loading a 100M chromosome into memory would take up more resources and
reduce speed.

However, for bacterial genome or metagenome sequences, loading a complete sequence has no impact,
and `hnsm range` will use the LRU cache to store the recently used sequences to reduce disk accesses
and thus speed up the process. In addition, plain text files use the same indexing format as BGZF.

```shell
# gz
bgzip -c tests/index/final.contigs.fa > tests/index/final.contigs.fa.gz;
bgzip -r tests/index/final.contigs.fa.gz

hnsm gz tests/index/final.contigs.fa -o tmp

# range
samtools faidx tests/index/final.contigs.fa
samtools faidx tests/index/final.contigs.fa \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"
samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg

hnsm range tests/index/final.contigs.fa.gz \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170(-):1-20" "k81_158:70001-70020"
hnsm range tests/index/final.contigs.fa.gz -r tests/index/sample.rg

```

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024.
