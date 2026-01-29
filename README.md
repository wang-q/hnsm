# hnsm - **H**omogeneous **N**ucleic acids/amino acids **S**mart **M**atching

[![Publish](https://github.com/wang-q/hnsm/actions/workflows/publish.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![Build](https://github.com/wang-q/hnsm/actions/workflows/build.yml/badge.svg)](https://github.com/wang-q/hnsm/actions)
[![codecov](https://codecov.io/gh/wang-q/hnsm/branch/master/graph/badge.svg?token=8toyNHCsVU)](https://codecov.io/gh/wang-q/hnsm)
[![Lines of code](https://www.aschey.tech/tokei/github/wang-q/hnsm)](https://github.com//wang-q/hnsm)
[![license](https://img.shields.io/github/license/wang-q/hnsm)](https://github.com//wang-q/hnsm)

`hnsm` is a high-performance bioinformatics toolkit for sequence analysis and comparison. It
provides:

- 🚀 Fast sequence comparison
- 🧠 Smart matching algorithms for sequence clustering
- 🧬 Support for both DNA and protein sequences
- 📊 Comprehensive analysis tools for genomic data
- ⚡️ Parallel processing for large datasets

<!-- TOC -->
* [hnsm - **H**omogeneous **N**ucleic acids/amino acids **S**mart **M**atching](#hnsm---homogeneous-nucleic-acidsamino-acids-smart-matching)
  * [Install](#install)
  * [Synopsis](#synopsis)
    * [`hnsm help`](#hnsm-help)
  * [Examples](#examples)
    * [FA files](#fa-files)
    * [Index](#index)
    * [FQ](#fq)
    * [Distances and Clustering](#distances-and-clustering)
      * [Similarity and dissimilarity (distance) of vectors](#similarity-and-dissimilarity-distance-of-vectors)
      * [Pairwise distances by Minimizer](#pairwise-distances-by-minimizer)
      * [Clustering](#clustering)
      * [PCoA](#pcoa)
    * [Genomes](#genomes)
    * [Proteomes](#proteomes)
    * [Assemblies](#assemblies)
  * [Author](#author)
  * [License](#license)
<!-- TOC -->

## Install

Current release: 0.3.4

```bash
cargo install --path . --force #--offline

# local repo
# cargo clean
# rm Cargo.lock

# test
cargo test -- --test-threads=1

# bench
cargo bench --bench simd

# local docs
# Linux
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps --document-private-items --open
# Windows Powershell
$env:RUSTDOCFLAGS="--html-in-header ./katex-header.html"
cargo doc --no-deps --document-private-items --open

# build under WSL 2
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build

```

## Synopsis

### `hnsm help`

```console
$ hnsm help
Homogeneous Nucleic acids/amino acids Smart Matching

Usage: hnsm [COMMAND]

Commands:
  size        Count total bases in FA file(s)
  count       Count base statistics in FA file(s)
  masked      Identify masked regions in FA file(s)
  n50         Calculate N50 and other assembly statistics
  one         Extract one FA record by name
  some        Extract some FA records based on a list of names
  order       Extract some FA records in the order specified by a list
  split       Split FA file(s) into several files
  replace     Replace headers of a FA file based on a TSV mapping
  rc          Reverse complement sequences in FA file(s)
  filter      Filter and format sequences in FA file(s)
  dedup       Deduplicate records in FA file(s)
  mask        Mask regions in FA file(s)
  sixframe    Translate DNA sequences in six frames
  gz          Compressing a file using the BGZF format
  range       Extract sequence regions by coordinates
  prefilter   Prefilter genome/metagenome assembly by amino acid minimizers
  interleave  Interleave paired-end sequences
  distance    Estimate sequence distances using minimizers
  hv          Estimate distances between DNA/protein files using hypervectors
  similarity  Calculate similarity between vectors
  manifold    Manifold learning based on pairwise distances
  clust       Clustering commands
  mat         Matrix commands
  das         Domain architecture similarity
  chain       Chains of syntenic genes
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version


Subcommand groups:

* Fasta files
    * info: size / count / masked / n50
    * records: one / some / order / split
    * transform: replace / rc / filter / dedup / mask / sixframe
    * indexing: gz / range / prefilter
* Fastq files: interleave

* Distance
    * DNA/protein: distance / hv
    * vectors: similarity
    * manifold
* Clustering: clust cc / clust dbscan
* Matrix: mat pair / mat phylip

* Synteny
    * das
    * chain

```

## Examples

### FA files

```bash
hnsm size tests/fasta/ufasta.fa
hnsm count tests/fasta/ufasta.fa.gz
hnsm masked tests/fasta/ufasta.fa
hnsm n50 tests/fasta/ufasta.fa -N 90 -N 50 -S -t

hnsm size tests/genome/mg1655.fa.gz -o tests/genome/mg1655.size.tsv

hnsm one tests/fasta/ufasta.fa read12
hnsm some tests/fasta/ufasta.fa tests/fasta/list.txt
hnsm order tests/fasta/ufasta.fa tests/fasta/list.txt

hnsm filter -a 10 -z 50 -U tests/fasta/ufasta.fa
hnsm filter -a 1 -u tests/fasta/ufasta.fa tests/fasta/ufasta.fa.gz
hnsm filter --iupac --upper tests/fasta/filter.fa

hnsm dedup tests/fasta/dedup.fa
hnsm dedup tests/fasta/dedup.fa -s -b -f stdout

hnsm mask --hard tests/fasta/ufasta.fa tests/fasta/mask.json

hnsm replace tests/fasta/ufasta.fa tests/fasta/replace.tsv
hnsm rc tests/fasta/ufasta.fa

hnsm filter -a 400 tests/fasta/ufasta.fa |
    hnsm split name stdin -o tmp
hnsm split about -c 2000 tests/fasta/ufasta.fa -o tmp

hnsm sixframe tests/fasta/trans.fa
hnsm sixframe tests/fasta/trans.fa --len 3 --start --end

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

```bash
# gz
bgzip -c tests/index/final.contigs.fa > tests/index/final.contigs.fa.gz;
bgzip -r tests/index/final.contigs.fa.gz

hnsm gz tests/index/final.contigs.fa -o tmp

# range
samtools faidx tests/index/final.contigs.fa
samtools faidx tests/index/final.contigs.fa \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170:1-20" "k81_158:70001-70020"
samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg

hnsm range tests/index/final.contigs.fa.gz \
    "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_170(-):1-20" "k81_158:70001-70020"
hnsm range tests/index/final.contigs.fa.gz -r tests/index/sample.rg

```

### FQ

```bash
hnsm interleave tests/fasta/ufasta.fa.gz tests/fasta/ufasta.fa

hnsm interleave tests/fasta/ufasta.fa

hnsm interleave --fq tests/fasta/ufasta.fa

hnsm interleave --fq tests/fastq/R1.fq.gz tests/fastq/R2.fq.gz

hnsm interleave --fq tests/fastq/R1.fq.gz

cargo run --bin hnsm fq2fa tests/fastq/R1.fq.gz

```

### Distances and Clustering

#### Similarity and dissimilarity (distance) of vectors

```bash
hnsm similarity tests/clust/domain.tsv --mode euclid --bin

hnsm similarity tests/clust/domain.tsv --mode cosine --bin

hnsm similarity tests/clust/domain.tsv --mode jaccard --bin

hyperfine --warmup 1 \
    -n p1 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 1 > /dev/null' \
    -n p2 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 2 > /dev/null' \
    -n p3 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 3 > /dev/null' \
    -n p4 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 4 > /dev/null' \
    -n p6 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 6 > /dev/null' \
    -n p8 \
    'hnsm similarity data/Domian_content_1000.tsv --mode jaccard --bin -p 8 > /dev/null' \
    --export-markdown sim.md.tmp

```

* AMD Ryzen 7 8745HS

| Command |       Mean [s] | Min [s] | Max [s] |    Relative |
|:--------|---------------:|--------:|--------:|------------:|
| `p1`    | 12.340 ± 0.270 |  11.762 |  12.697 | 3.79 ± 0.15 |
| `p2`    |  6.952 ± 0.219 |   6.553 |   7.268 | 2.14 ± 0.10 |
| `p3`    |  5.383 ± 0.344 |   5.171 |   6.287 | 1.66 ± 0.12 |
| `p4`    |  4.375 ± 0.162 |   4.046 |   4.585 | 1.35 ± 0.07 |
| `p6`    |  3.602 ± 0.069 |   3.445 |   3.712 | 1.11 ± 0.04 |
| `p8`    |  3.252 ± 0.106 |   2.984 |   3.390 |        1.00 |

#### Pairwise distances by Minimizer

```console
$ hnsm distance tests/clust/IBPA.fa -k 7 -w 1 |
    hnsm mat phylip stdin
  10
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

```console
$ clustalo -i tests/clust/IBPA.fa --auto --full --distmat-out=tests/clust/IBPA.phy
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

#### Clustering

```bash
hnsm clust dbscan tests/clust/IBPA.fa.tsv --eps 0.05 --min_points 2
hnsm clust dbscan tests/clust/IBPA.fa.tsv --eps 0.05 --min_points 2 --format pair

cat tests/clust/IBPA.fa.tsv |
    tsv-filter --le 3:0.05 |
    tsv-select -f 1-2 \
    > tests/clust/IBPA.fa.05.tsv

hnsm clust cc tests/clust/IBPA.fa.05.tsv

```

#### PCoA

```bash
cargo run --bin hnsm manifold tests/clust/IBPA.fa.tsv --mode pcoa --dim 2

```

### Genomes

* genomes

```bash
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz |
    gzip -dc |
    hnsm filter stdin -s |
    hnsm gz stdin -o tests/genome/mg1655.fa

curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz |
    gzip -dc |
    hnsm filter stdin -s |
    hnsm gz stdin -o tests/genome/sakai.fa

```

* mash

```bash
hnsm distance tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz --hasher mod -k 21 -w 1
#NC_002695       NC_000913       0.0221  0.4580  0.5881
#NC_002127       NC_000913       0.6640  0.0000  0.0006
#NC_002128       NC_000913       0.4031  0.0001  0.0053

hnsm rc tests/genome/mg1655.fa.gz |
    hnsm distance tests/genome/sakai.fa.gz stdin --hasher mod -k 21 -w 1
#NC_002695       RC_NC_000913    0.0221  0.4580  0.5881
#NC_002127       RC_NC_000913    0.6640  0.0000  0.0006
#NC_002128       RC_NC_000913    0.4031  0.0001  0.0053

hnsm rc tests/genome/mg1655.fa.gz |
    hnsm distance tests/genome/mg1655.fa.gz stdin --hasher mod -k 21 -w 1
#NC_000913       RC_NC_000913    0.0000  1.0000  1.0000
hnsm rc tests/genome/mg1655.fa.gz |
    hnsm distance tests/genome/mg1655.fa.gz stdin --hasher rapid -k 21 -w 1
#NC_000913       RC_NC_000913    0.2289  0.0041  0.0082

hnsm distance tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz --merge --hasher mod -k 21 -w 1
#tests/genome/sakai.fa.gz   tests/genome/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779

hnsm distance tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz --merge --hasher rapid -k 21 -w 1
#tests/genome/sakai.fa.gz   tests/genome/mg1655.fa.gz  5394043 4562542 3071076 6885509 0.0230  0.4460  0.5693

echo -e "tests/genome/sakai.fa.gz\ntests/genome/mg1655.fa.gz" |
    hnsm distance stdin --merge --list --hasher mod -k 21 -w 1
#tests/genome/sakai.fa.gz   tests/genome/sakai.fa.gz   5302382 5302382 5302382 5302382 0.0000  1.0000  1.0000
#tests/genome/sakai.fa.gz   tests/genome/mg1655.fa.gz  5302382 4543891 3064483 6781790 0.0226  0.4519  0.5779
#tests/genome/mg1655.fa.gz  tests/genome/sakai.fa.gz   4543891 5302382 3064483 6781790 0.0226  0.4519  0.6744
#tests/genome/mg1655.fa.gz  tests/genome/mg1655.fa.gz  4543891 4543891 4543891 4543891 0.0000  1.0000  1.0000

```

* synteny


```bash
hnsm synt dna tests/genome/small_1.fa tests/genome/small_2.fa -k 21

hnsm synt dna tests/genome/small_1.fa tests/genome/small_2.fa -k 21 -o tests/synt/small_1_2.tsv

hnsm synt dna tests/genome/small_1.fa tests/genome/small_2.fa -k 21 --rounds 500,10 --chain-gap 500 -o tests/synt/small_1_2_500.tsv

# Merge fragmented blocks
hnsm synt merge tests/synt/small_1_2.tsv -o tests/synt/small_1_2.merged.tsv

# Merge with divergence-based parameters (e.g., 5% divergence -> chain-gap 100000)
hnsm synt merge tests/synt/small_1_2_500.tsv -d 5.0 -o tests/synt/small_1_2.d5.tsv

cargo run --bin hnsm synt view tests/synt/small_1_2.tsv -o tmp.svg

# Read synteny blocks
hnsm synt dna tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz -k 21 --chain-gap 100 --min-weight 2 --max-freq 100 --rounds 1000,100,10 -v -o tests/synt/mg1655_sakai.tsv

hnsm synt merge tests/synt/mg1655_sakai.tsv --chain-gap 1000 -o tests/synt/mg1655_sakai.1000.tsv

hnsm synt view tests/synt/mg1655_sakai.tsv -o tests/synt/mg1655_sakai.svg

```

### Proteomes

```bash
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz \
    > tests/clust/mg1655.pro.fa.gz

curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_protein.faa.gz \
    > tests/clust/sakai.pro.fa.gz

curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_protein.faa.gz \
    > tests/clust/pao1.pro.fa.gz

hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 |
    rgr filter stdin --ne 3:1

hyperfine --warmup 1 \
    -n p1 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 1 > /dev/null' \
    -n p2 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 2 > /dev/null' \
    -n p3 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 3 > /dev/null' \
    -n p4 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 4 > /dev/null' \
    -n p6 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 6 > /dev/null' \
    -n p8 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 8 > /dev/null' \
    --export-markdown dis.md.tmp

hyperfine --warmup 1 \
    -n p1 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 1 > /dev/null' \
    -n p2 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 2 > /dev/null' \
    -n p3 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 3 > /dev/null' \
    -n p4 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 4 > /dev/null' \
    -n p6 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 6 > /dev/null' \
    -n p8 \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 --zero -p 8 > /dev/null' \
    --export-markdown dis.md.tmp

```

* AMD Ryzen 7 8745HS

| Command |       Mean [s] | Min [s] | Max [s] |    Relative |
|:--------|---------------:|--------:|--------:|------------:|
| `p1`    | 24.840 ± 0.377 |  24.161 |  25.319 | 8.36 ± 0.45 |
| `p2`    | 11.150 ± 0.153 |  10.958 |  11.481 | 3.75 ± 0.20 |
| `p3`    |  6.877 ± 0.083 |   6.743 |   7.022 | 2.31 ± 0.12 |
| `p4`    |  5.178 ± 0.107 |   5.033 |   5.317 | 1.74 ± 0.10 |
| `p6`    |  3.745 ± 0.205 |   3.508 |   4.080 | 1.26 ± 0.10 |
| `p8`    |  2.972 ± 0.155 |   2.842 |   3.347 |        1.00 |

| Command |       Mean [s] | Min [s] | Max [s] |    Relative |
|:--------|---------------:|--------:|--------:|------------:|
| `p1`    | 32.044 ± 0.776 |  30.726 |  32.970 | 8.78 ± 0.25 |
| `p2`    | 15.014 ± 0.110 |  14.842 |  15.212 | 4.12 ± 0.07 |
| `p3`    |  9.182 ± 0.080 |   8.987 |   9.246 | 2.52 ± 0.04 |
| `p4`    |  6.885 ± 0.087 |   6.735 |   6.974 | 1.89 ± 0.04 |
| `p6`    |  4.569 ± 0.110 |   4.492 |   4.831 | 1.25 ± 0.04 |
| `p8`    |  3.648 ± 0.055 |   3.587 |   3.774 |        1.00 |

* Hypervector

```bash
hnsm hv tests/clust/IBPA.fa
#tests/clust/IBPA.fa     tests/clust/IBPA.fa     776     776     776     776     0.0000  1.0000  1.0000
hnsm distance tests/clust/IBPA.fa --merge
#tests/clust/IBPA.fa     tests/clust/IBPA.fa     763     763     763     763     0.0000  1.0000  1.0000

hnsm hv tests/clust/mg1655.pro.fa.gz
#tests/clust/mg1655.pro.fa.gz    tests/clust/mg1655.pro.fa.gz    1240734 1240734 1240734 1240734 0.0000  1.0000  1.0000
hnsm distance tests/clust/mg1655.pro.fa.gz --merge
#tests/clust/mg1655.pro.fa.gz    tests/clust/mg1655.pro.fa.gz    1267403 1267403 1267403 1267403 0.0000  1.0000  1.0000

hnsm hv tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 1
#tests/clust/mg1655.pro.fa.gz    tests/clust/pao1.pro.fa.gz      1240734 1733273 81195   2892811 0.4154  0.0281  0.0654
hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 1 --merge
#tests/clust/mg1655.pro.fa.gz    tests/clust/pao1.pro.fa.gz      1267403 1770832 60605   2977630 0.4602  0.0204  0.0478

hyperfine --warmup 1 \
    -n distance \
    'hnsm distance tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 1 --merge > /dev/null' \
    -n hv \
    'hnsm hv tests/clust/mg1655.pro.fa.gz tests/clust/pao1.pro.fa.gz -k 7 -w 2 -p 1 > /dev/null' \
    --export-markdown dis.md.tmp

```

| Command    |    Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:-----------|-------------:|---------:|---------:|------------:|
| `distance` |  161.5 ± 3.3 |    154.6 |    167.4 |        1.00 |
| `hv`       | 1571.1 ± 4.4 |   1564.7 |   1577.3 | 9.73 ± 0.20 |

* Six-frame

```bash
hnsm sixframe tests/genome/sakai.fa.gz --len 35 |
    hnsm distance stdin tests/clust/mg1655.pro.fa.gz -k 7 -w 2 -p 4 |
    wc -l
#21124

hnsm range tests/genome/sakai.fa.gz "NC_002695(+):4468532-4468696|frame=1"

```

### Assemblies

```bash
cargo run --bin hnsm prefilter tests/index/final.contigs.fa tests/clust/IBPA.fa

# SRR6323163 - APH(3')-IIIa
# 3300030246 - acrB
hnsm prefilter tests/clust/SRR6323163.fa.gz "tests/clust/APH(3')-IIIa.fa"
hnsm prefilter tests/clust/SRR6323163.fa.gz "tests/clust/acrB.fa"

hnsm range tests/clust/SRR6323163.fa.gz "k141_4576(-):285-455|frame=2"

hnsm prefilter 3300035148.fna.gz "tests/clust/APH(3')-IIIa.fa" -c 1000000 -p 8

```

## Author

Qiang Wang <wang-q@outlook.com>

## License

MIT.

Copyright by Qiang Wang.

Written by Qiang Wang <wang-q@outlook.com>, 2024-
