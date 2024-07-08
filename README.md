# hnsm

[![Build](https://github.com/wang-q/hnsm/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![codecov](https://codecov.io/gh/wang-q/hnsm/branch/main/graph/badge.svg?token=LtxYK5Fff0)](https://codecov.io/gh/wang-q/hnsm)
[![license](https://img.shields.io/github/license/wang-q/hnsm)](https://github.com//wang-q/hnsm)
[![Lines of code](https://tokei.rs/b1/github/wang-q/hnsm?category=code)](https://github.com//wang-q/hnsm)

`hnsm` - **H**omogeneous **N**ucleic acid **S**mart **M**atching

## Install

Current release: 0.0.1

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
  masked    Masked regions in FA file(s)
  one       Extract one FA record
  order     Extract some FA records by the given order
  rc        Reverse complement a FA file
  sixframe  Six-Frame Translation
  size      Count total bases in FA file(s)
  some      Extract some FA records
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version

```

## Examples

```shell
cargo run --bin hnsm help

hnsm size tests/fasta/ufasta.fa
hnsm count tests/fasta/ufasta.fa

hnsm one tests/fasta/ufasta.fa read12

hnsm some tests/fasta/ufasta.fa tests/fasta/list.txt

cargo run --bin hnsm order tests/fasta/ufasta.fa tests/fasta/list.txt

hnsm masked tests/fasta/ufasta.fa

hnsm rc tests/fasta/ufasta.fa


cargo run --bin hnsm sixframe

```

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024.
