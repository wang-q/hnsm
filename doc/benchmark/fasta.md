# `hnsm` and `faops`

Benchmarks between C and Rust implementations

```shell
brew install faops

brew install hyperfine
brew install neofetch

```

## System info

* Ryzen 7 5800

```text
wangq@R7
--------
OS: Ubuntu 20.04.5 LTS on Windows 10 x86_64
Kernel: 5.15.153.1-microsoft-standard-WSL2
Uptime: 16 days, 6 hours, 30 mins
Packages: 1408 (dpkg), 191 (brew), 5 (snap)
Shell: bash 5.0.17
Theme: Adwaita [GTK3]
Icons: Adwaita [GTK3]
Terminal: Windows Terminal
CPU: AMD Ryzen 7 5800 (16) @ 3.393GHz
GPU: fb45:00:00.0 Microsoft Corporation Device 008e
Memory: 562MiB / 32030MiB

```

## `hnsm size`

```shell
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "hnsm  size .fa" \
    'cat tests/fasta/ufasta.fa | hnsm  size stdin > /dev/null' \
    -n "faops size .fa" \
    'cat tests/fasta/ufasta.fa | faops size stdin > /dev/null' \
    -n "hnsm  size .fa.gz" \
    'hnsm  size tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops size .fa.gz" \
    'faops size tests/fasta/ufasta.fa.gz > /dev/null'

cat size.md.tmp

```

| Command             | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------|----------:|---------:|---------:|------------:|
| `hnsm  size .fa`    | 2.7 ± 0.3 |      2.1 |      4.2 | 1.01 ± 0.13 |
| `faops size .fa`    | 2.6 ± 0.2 |      2.0 |      4.2 |        1.00 |
| `hnsm  size .fa.gz` | 2.7 ± 0.2 |      2.3 |      3.6 | 1.03 ± 0.12 |
| `faops size .fa.gz` | 2.8 ± 0.2 |      2.4 |      4.5 | 1.07 ± 0.13 |

## `hnsm some`

```shell
hyperfine --warmup 10 --export-markdown some.md.tmp \
    -n "hnsm some" \
    'hnsm  some tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some" \
    'faops some tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null' \
    -n "hnsm some -i" \
    'hnsm  some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some -i" \
    'faops some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null'

cat some.md.tmp

```

| Command         | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:----------------|----------:|---------:|---------:|------------:|
| `hnsm some`     | 3.9 ± 0.3 |      3.3 |      6.0 | 1.01 ± 0.12 |
| `faops some`    | 4.1 ± 0.3 |      3.6 |      5.9 | 1.06 ± 0.12 |
| `hnsm some -i`  | 3.9 ± 0.3 |      3.4 |      5.8 |        1.00 |
| `faops some -i` | 4.1 ± 0.4 |      3.6 |      6.3 | 1.07 ± 0.13 |

## `hnsm n50`

```shell
hyperfine --warmup 10 --export-markdown n50.md.tmp \
    -n "hnsm n50 .gz" \
    'hnsm  n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops n50 .gz" \
    'faops n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "hnsm n50 -E -S -A" \
    'hnsm  n50 -E -S -A tests/fasta/ufasta.fa > /dev/null' \
    -n "faops n50 -E -S -A" \
    'faops n50 -E -S -A tests/fasta/ufasta.fa > /dev/null'

cat n50.md.tmp

```

| Command              | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:---------------------|----------:|---------:|---------:|------------:|
| `hnsm n50 .gz`       | 2.7 ± 0.2 |      2.3 |      3.7 | 1.00 ± 0.09 |
| `faops n50 .gz`      | 2.7 ± 0.2 |      2.3 |      3.6 |        1.00 |
| `hnsm n50 -E -S -A`  | 3.2 ± 0.2 |      2.9 |      4.2 | 1.22 ± 0.10 |
| `faops n50 -E -S -A` | 2.8 ± 0.2 |      2.4 |      3.8 | 1.04 ± 0.09 |
