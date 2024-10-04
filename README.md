# hnsm

[![Build](https://github.com/wang-q/hnsm/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![codecov](https://codecov.io/gh/wang-q/hnsm/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/hnsm)
[![license](https://img.shields.io/github/license/wang-q/hnsm)](https://github.com//wang-q/hnsm)

`hnsm` - **H**omogeneous **N**ucleic acids **S**mart **M**atching

## Install

Current release: 0.1.8

```shell
rustup update -- nightly
cargo +nightly install --path . --force #--offline

# local repo
# cargo clean
# rm Cargo.lock

# test
cargo +nightly test -- --test-threads=1

# build under WSL 2
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo +nightly build

# build for CentOS 7
# rustup target add x86_64-unknown-linux-gnu
# pip3 install cargo-zigbuild
cargo +nightly zigbuild --target x86_64-unknown-linux-gnu.2.17 --release
ll $CARGO_TARGET_DIR/x86_64-unknown-linux-gnu/release/

```

## Synopsis

### `hnsm help`

```text
Homogeneous Nucleic acids Smart Matching

Usage: hnsm [COMMAND]

Commands:
  count       Count base statistics in FA file(s)
  distance    Estimate distances between DNA/protein sequences
  filter      Filter records in FA file(s)
  gz          Compressing a file using the blocked gzip format (BGZF)
  masked      Masked regions in FA file(s)
  n50         Count total bases in FA file(s)
  one         Extract one FA record
  order       Extract some FA records by the given order
  range       Extract sequences defined by the range(s)
  rc          Reverse complement a FA file
  replace     Replace headers of a FA file
  similarity  Similarity of vectors
  sixframe    Six-Frame Translation
  size        Count total bases in FA file(s)
  some        Extract some FA records
  split       Split FA file(s) into several files
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version


Subcommand groups:

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter
    * indexing: gz / range

* Clustering
    * vectors: similarity
    * DNA/protein: distance / identity

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

hnsm filter -a 10 -z 50 -U tests/fasta/ufasta.fa
hnsm filter -a 1 -u tests/fasta/ufasta.fa tests/fasta/ufasta.fa.gz
hnsm filter --iupac --upper tests/fasta/filter.fa

hnsm replace tests/fasta/ufasta.fa tests/fasta/replace.tsv
hnsm rc tests/fasta/ufasta.fa

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

### Clustering

#### Similarity and dissimilarity (distance) of vectors

```shell
# rustup update -- nightly
cargo +nightly bench --bench simd

cargo run --bin nwr similarity tests/assembly/domain.tsv --mode euclid --bin

cargo run --bin nwr similarity tests/assembly/domain.tsv --mode cosine --bin

cargo run --bin nwr similarity tests/assembly/domain.tsv --mode jaccard --bin

hyperfine --warmup 1 \
    -n p1 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 1 --mode jaccard --bin > /dev/null' \
    -n p2 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 2 --mode jaccard --bin > /dev/null' \
    -n p3 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 3 --mode jaccard --bin > /dev/null' \
    -n p4 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 4 --mode jaccard --bin > /dev/null' \
    -n p6 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 6 --mode jaccard --bin > /dev/null' \
    -n p8 \
    'nwr similarity data/Domian_content_1000.tsv --parallel 8 --mode jaccard --bin > /dev/null' \
    --export-markdown sim.md.tmp

```

| Command |       Mean [s] | Min [s] | Max [s] |    Relative |
|:--------|---------------:|--------:|--------:|------------:|
| `p1`    | 17.364 ± 1.244 |  16.065 |  20.367 | 2.11 ± 0.18 |
| `p2`    | 10.467 ± 1.405 |   9.421 |  14.045 | 1.27 ± 0.18 |
| `p3`    |  8.226 ± 0.380 |   7.615 |   8.722 |        1.00 |
| `p4`    |  8.430 ± 0.842 |   7.777 |  10.614 | 1.02 ± 0.11 |
| `p6`    |  8.330 ± 0.827 |   7.648 |  10.371 | 1.01 ± 0.11 |
| `p8`    | 10.268 ± 2.407 |   8.415 |  15.486 | 1.25 ± 0.30 |

| Command |       Mean [s] | Min [s] | Max [s] |    Relative |
|:--------|---------------:|--------:|--------:|------------:|
| `p1`    | 27.631 ± 0.533 |  27.123 |  29.022 | 6.92 ± 0.17 |
| `p2`    | 14.673 ± 0.241 |  14.417 |  15.138 | 3.67 ± 0.08 |
| `p3`    |  9.778 ± 0.075 |   9.635 |   9.885 | 2.45 ± 0.04 |
| `p4`    |  7.472 ± 0.160 |   7.260 |   7.686 | 1.87 ± 0.05 |
| `p6`    |  5.151 ± 0.123 |   4.975 |   5.397 | 1.29 ± 0.04 |
| `p8`    |  3.995 ± 0.062 |   3.915 |   4.105 |        1.00 |

#### Pairwise distances computed by MEGA

```shell
cargo run --bin hnsm distance tests/fasta/IBPA.fa -k 7 -w 1

# distance matrix
brew install csvtk
brew install wang-q/tap/tsv-utils
cargo install affinityprop

cargo run --bin hnsm distk tests/fasta/IBPA.fa -k 7 -w 1 --sim |
    tsv-select -f 1-3 |
    csvtk spread -H -t -k 2 -v 3 |
    sed '1d' \
    > tests/fasta/IBPA.fa.sim

affinityprop -s 3 --damping 0.1 --input tests/fasta/IBPA.fa.sim

```

```text
[ 1] #IBPA_ECOLI
[ 2] #IBPA_ECOLI_GA
[ 3] #IBPA_ECOLI_GA_LV
[ 4] #IBPA_ECOLI_GA_LV_ST
[ 5] #IBPA_ECOLI_GA_LV_RK
[ 6] #IBPA_ESCF3
[ 7] #A0A192CFC5_ECO25
[ 8] #Q2QJL7_ACEAC
[ 9] #A0A010SUI8_PSEFL
[10] #K1J4J6_9GAMM

[         1      2      3      4      5      6      7      8      9     10 ]
[ 1]
[ 2]  0.0602
[ 3]  0.1750 0.1078
[ 4]  0.2195 0.1493 0.0372
[ 5]  0.3249 0.2472 0.1242 0.0837
[ 6]  0.0000 0.0602 0.1750 0.2195 0.3249
[ 7]  0.0000 0.0602 0.1750 0.2195 0.3249 0.0000
[ 8]  0.8522 0.9614 1.0840 1.0625 1.1991 0.8522 0.8522
[ 9]  0.7768 0.8595 1.0080 1.0282 1.1133 0.7768 0.7841 0.9763
[10]  0.4583 0.5306 0.6785 0.7230 0.8351 0.4583 0.4583 0.7230 0.7153

```

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024.
