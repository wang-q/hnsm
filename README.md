# hnsm

[![Build](https://github.com/wang-q/hnsm/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![codecov](https://codecov.io/gh/wang-q/hnsm/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/hnsm)
[![license](https://img.shields.io/github/license/wang-q/hnsm)](https://github.com//wang-q/hnsm)

`hnsm` - **H**omogeneous **N**ucleic acids/amino acids **S**mart **M**atching

<!-- TOC -->
* [hnsm](#hnsm)
  * [Install](#install)
  * [Synopsis](#synopsis)
    * [`hnsm help`](#hnsm-help)
    * [`fasr help`](#fasr-help)
  * [Examples](#examples)
    * [Fasta files](#fasta-files)
    * [Index](#index)
    * [Fastq](#fastq)
    * [Clustering](#clustering)
      * [Similarity and dissimilarity (distance) of vectors](#similarity-and-dissimilarity-distance-of-vectors)
      * [Pairwise distances by Minimizer](#pairwise-distances-by-minimizer)
      * [Matrix conversion](#matrix-conversion)
      * [DBSCAN](#dbscan)
      * [PCoA](#pcoa)
    * [Block Fasta files](#block-fasta-files)
  * [Author](#author)
  * [License](#license)
<!-- TOC -->

## Install

Current release: 0.1.11

```shell
cargo install --path . --force #--offline

# local repo
# cargo clean
# rm Cargo.lock

# test
cargo test -- --test-threads=1

# build under WSL 2
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build

```

## Synopsis

### `hnsm help`

```text
Homogeneous Nucleic acids Smart Matching

Usage: hnsm [COMMAND]

Commands:
  cluster     Clustering based on pairwise distances
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
    * cluster

* <infiles> are paths to fasta files, .fa.gz is supported
    * infile == stdin means reading from STDIN
    * `hnsm gz` writes out the BGZF format and `hnsm range` reads it

```

### `fasr help`

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

hnsm dedup tests/fasta/dedup.fa
hnsm dedup tests/fasta/dedup.fa -s -b -f stdout

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

### Fastq

```shell
hnsm interleave tests/fasta/ufasta.fa.gz tests/fasta/ufasta.fa

hnsm interleave tests/fasta/ufasta.fa

hnsm interleave --fq tests/fasta/ufasta.fa

hnsm interleave --fq tests/fastq/R1.fq.gz tests/fastq/R2.fq.gz

hnsm interleave --fq tests/fastq/R1.fq.gz

```

### Clustering

#### Similarity and dissimilarity (distance) of vectors

```shell
cargo bench --bench simd

cargo run --bin hnsm similarity tests/clust/domain.tsv --mode euclid --bin

cargo run --bin hnsm similarity tests/clust/domain.tsv --mode cosine --bin

cargo run --bin hnsm similarity tests/clust/domain.tsv --mode jaccard --bin

hyperfine --warmup 1 \
    -n p1 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 1 --mode jaccard --bin > /dev/null' \
    -n p2 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 2 --mode jaccard --bin > /dev/null' \
    -n p3 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 3 --mode jaccard --bin > /dev/null' \
    -n p4 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 4 --mode jaccard --bin > /dev/null' \
    -n p6 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 6 --mode jaccard --bin > /dev/null' \
    -n p8 \
    'hnsm similarity data/Domian_content_1000.tsv --parallel 8 --mode jaccard --bin > /dev/null' \
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

#### Pairwise distances by Minimizer

```text
$ hnsm distance tests/clust/IBPA.fa -k 7 -w 1 |
    hnsm cluster stdin --mode matrix
IBPA_ECOLI      0       0.0669  0.2014  0.2106  0.4405  0       0.0011  0.7124  0.5454  0.3675
IBPA_ECOLI_GA   0.0669  0       0.114   0.1464  0.2917  0.0669  0.068   1       0.7023  0.5974
IBPA_ECOLI_GA_LV        0.2014  0.114   0       0.0344  0.132   0.2014  0.2024  1       0.7023  1
IBPA_ECOLI_GA_LV_ST     0.2106  0.1464  0.0344  0       0.0895  0.2106  0.2117  1       0.7023  1
IBPA_ECOLI_GA_LV_RK     0.4405  0.2917  0.132   0.0895  0       0.4405  0.4416  1       0.7023  1
IBPA_ESCF3      0       0.0669  0.2014  0.2106  0.4405  0       0.0011  0.7124  0.5454  0.3675
A0A192CFC5_ECO25        0.0011  0.068   0.2024  0.2117  0.4416  0.0011  0       0.7134  0.5464  0.3686
Q2QJL7_ACEAC    0.7124  1       1       1       1       0.7124  0.7134  0       1       0.6134
A0A010SUI8_PSEFL        0.5454  0.7023  0.7023  0.7023  0.7023  0.5454  0.5464  1       0       0.7023
K1J4J6_9GAMM    0.3675  0.5974  1       1       1       0.3675  0.3686  0.6134  0.7023  0

```

```text
$ clustalo -i tests/clust/IBPA.fa --auto --full --distmat-out=tests/clust/IBPA.mat
10
IBPA_ECOLI          0.000000 0.058394 0.160584 0.197080 0.277372 0.000000 0.000000 0.583942 0.540146 0.372263
IBPA_ECOLI_GA       0.058394 0.000000 0.102190 0.138686 0.218978 0.058394 0.058394 0.627737 0.576642 0.416058
IBPA_ECOLI_GA_LV    0.160584 0.102190 0.000000 0.036496 0.116788 0.160584 0.160584 0.671533 0.635036 0.496350
IBPA_ECOLI_GA_LV_ST 0.197080 0.138686 0.036496 0.000000 0.080292 0.197080 0.197080 0.656934 0.642336 0.518248
IBPA_ECOLI_GA_LV_RK 0.277372 0.218978 0.116788 0.080292 0.000000 0.277372 0.277372 0.700730 0.671533 0.569343
IBPA_ESCF3          0.000000 0.058394 0.160584 0.197080 0.277372 0.000000 0.000000 0.583942 0.540146 0.372263
A0A192CFC5_ECO25    0.000000 0.058394 0.160584 0.197080 0.277372 0.000000 0.000000 0.589928 0.546763 0.372263
Q2QJL7_ACEAC        0.583942 0.627737 0.671533 0.656934 0.700730 0.583942 0.589928 0.000000 0.628378 0.518248
A0A010SUI8_PSEFL    0.540146 0.576642 0.635036 0.642336 0.671533 0.540146 0.546763 0.628378 0.000000 0.496350
K1J4J6_9GAMM        0.372263 0.416058 0.496350 0.518248 0.569343 0.372263 0.372263 0.518248 0.496350 0.000000

```

#### Matrix conversion

```shell
hnsm convert tests/clust/IBPA.fa.tsv --mode matrix

hnsm convert tests/clust/IBPA.fa.tsv --mode lower

hnsm convert tests/clust/IBPA.mat --mode pair

```

#### DBSCAN

```shell
hnsm cluster tests/clust/IBPA.fa.tsv --mode dbscan --eps 0.05 --min_points 2

cat tests/clust/IBPA.fa.tsv |
    tsv-filter --le 3:0.05 |
    hnsm cluster stdin --mode cc

```

#### PCoA

```shell
cargo run --bin hnsm manifold tests/clust/IBPA.fa.tsv --mode pcoa --dim 2

```

### Block Fasta files

```shell
fasr maf2fas tests/fasr/example.maf

fasr axt2fas tests/fasr/RM11_1a.chr.sizes tests/fasr/example.axt --qname RM11_1a

fasr filter tests/fasr/example.fas --ge 10

fasr name tests/fasr/example.fas --count

fasr cover tests/fasr/example.fas

fasr cover tests/fasr/example.fas --name S288c --trim 10

fasr concat tests/fasr/name.lst tests/fasr/example.fas

fasr subset tests/fasr/name.lst tests/fasr/example.fas
fasr subset tests/fasr/name.lst tests/fasr/refine.fas --required

fasr link tests/fasr/example.fas --pair
fasr link tests/fasr/example.fas --best

fasr replace tests/fasr/replace.tsv tests/fasr/example.fas
fasr replace tests/fasr/replace.fail.tsv tests/fasr/example.fas

samtools faidx tests/fasr/NC_000932.fa NC_000932:1-10

fasr check tests/fasr/NC_000932.fa tests/fasr/A_tha.pair.fas

fasr create tests/fasr/genome.fa tests/fasr/I.connect.tsv --name S288c

# Create a fasta file containing multiple genomes
cat tests/fasr/genome.fa | sed 's/^>/>S288c./' > tests/fasr/genomes.fa
samtools faidx tests/fasr/genomes.fa S288c.I:1-100

cargo run --bin fasr create tests/fasr/genomes.fa tests/fasr/I.name.tsv --multi

fasr separate tests/fasr/example.fas -o . --suffix .tmp

spoa tests/fasr/refine.fasta -r 1

cargo run --bin fasr consensus tests/fasr/example.fas
cargo run --bin fasr consensus tests/fasr/refine.fas
cargo run --bin fasr consensus tests/fasr/refine.fas --outgroup -p 2

cargo run --bin fasr refine tests/fasr/example.fas
cargo run --bin fasr refine tests/fasr/example.fas --msa none --chop 10
cargo run --bin fasr refine tests/fasr/refine2.fas --msa clustalw --outgroup
cargo run --bin fasr refine tests/fasr/example.fas --quick

cargo run --bin fasr split tests/fasr/example.fas --simple
cargo run --bin fasr split tests/fasr/example.fas -o . --chr --suffix .tmp

cargo run --bin fasr slice tests/fasr/slice.json tests/fasr/slice.fas --name S288c

cargo run --bin fasr join tests/fasr/S288cvsYJM789.slice.fas --name YJM789
cargo run --bin fasr join \
    tests/fasr/S288cvsRM11_1a.slice.fas \
    tests/fasr/S288cvsYJM789.slice.fas \
    tests/fasr/S288cvsSpar.slice.fas

cargo run --bin fasr stat tests/fasr/example.fas --outgroup

cargo run --bin fasr variation tests/fasr/example.fas
cargo run --bin fasr variation tests/fasr/example.fas --outgroup

cargo run --bin fasr xlsx tests/fasr/example.fas
cargo run --bin fasr xlsx tests/fasr/example.fas --outgroup

cargo run --bin fasr pl-p2m tests/fasr/S288cvsRM11_1a.slice.fas tests/fasr/S288cvsSpar.slice.fas

```

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024.
