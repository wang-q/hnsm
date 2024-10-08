# `hnsm` and `faops`

```shell
brew install samtools

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

## `bgzip`

```shell
hyperfine --warmup 5 --export-markdown gz.md.tmp \
    -n "bgzip" \
    'rm -f tests/index/final.contigs.fa.gz*;
     bgzip -c tests/index/final.contigs.fa > tests/index/final.contigs.fa.gz;
     bgzip -r tests/index/final.contigs.fa.gz' \
    -n "bgzip --threads 4" \
    'rm -f tests/index/final.contigs.fa.gz*;
     bgzip -c --threads 4 tests/index/final.contigs.fa > tests/index/final.contigs.fa.gz;
     bgzip -r tests/index/final.contigs.fa.gz' \
    -n "hnsm gz" \
    'rm -f tests/index/final.contigs.fa.gz*;
     hnsm gz tests/index/final.contigs.fa' \
    -n "hnsm gz -p 4" \
    'rm -f tests/index/final.contigs.fa.gz*;
     hnsm gz -p 4 tests/index/final.contigs.fa'

cat gz.md.tmp

```

| Command             |   Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------|------------:|---------:|---------:|------------:|
| `bgzip`             | 134.1 ± 1.6 |    132.0 |    138.2 | 2.70 ± 0.07 |
| `bgzip --threads 4` |  67.1 ± 0.9 |     65.8 |     68.9 | 1.35 ± 0.04 |
| `hnsm gz`           | 95.7 ± 16.1 |     90.6 |    183.4 | 1.93 ± 0.33 |
| `hnsm gz -p 4`      |  49.7 ± 1.2 |     47.5 |     52.8 |        1.00 |

## `.loc` and `.fai`

```shell
hyperfine --warmup 5 --export-markdown faidx.md.tmp \
    -n "samtools faidx .fa" \
    'rm -f tests/index/final.contigs.fa.fai;
     samtools faidx tests/index/final.contigs.fa' \
    -n "samtools faidx .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.fai;
     samtools faidx tests/index/final.contigs.fa.gz' \
    -n "hnsm range .fa" \
    'rm -f tests/index/final.contigs.fa.loc;
     hnsm range tests/index/final.contigs.fa' \
    -n "hnsm range .fa.gz" \
    'rm -f tests/index/final.contigs.fa.gz.loc;
     hnsm range tests/index/final.contigs.fa.gz'

cat faidx.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 22.9 ± 1.5 |     20.1 |     27.8 |        1.00 |
| `samtools faidx .fa.gz` | 30.3 ± 1.5 |     27.6 |     36.1 | 1.32 ± 0.11 |
| `hnsm range .fa`        | 25.7 ± 1.4 |     23.2 |     33.2 | 1.12 ± 0.10 |
| `hnsm range .fa.gz`     | 25.2 ± 1.6 |     23.0 |     35.9 | 1.10 ± 0.10 |

## `hnsm range`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "hnsm range .fa" \
    'hnsm range tests/index/final.contigs.fa \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"' \
    -n "hnsm range .fa.gz" \
    'hnsm range tests/index/final.contigs.fa.gz \
        "k81_130" "k81_130:11-20" "k81_170:304-323" "k81_158:70001-70020"'

cat rg.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 12.1 ± 0.8 |     10.3 |     14.9 | 1.21 ± 0.11 |
| `samtools faidx .fa.gz` | 14.3 ± 0.8 |     12.6 |     17.1 | 1.42 ± 0.12 |
| `hnsm range .fa`        | 10.0 ± 0.7 |      8.7 |     12.8 |        1.00 |
| `hnsm range .fa.gz`     | 14.1 ± 0.7 |     12.0 |     15.8 | 1.40 ± 0.12 |

## `hnsm range -r`

```shell
hyperfine --warmup 5 --export-markdown rg.md.tmp \
    -n "samtools faidx .fa" \
    'samtools faidx tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "samtools faidx .fa.gz" \
    'samtools faidx tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null' \
    -n "hnsm range .fa" \
    'hnsm range tests/index/final.contigs.fa -r tests/index/sample.rg > /dev/null' \
    -n "hnsm range .fa.gz" \
    'hnsm range tests/index/final.contigs.fa.gz -r tests/index/sample.rg > /dev/null'

cat rg.md.tmp

```

| Command                 |  Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:------------------------|-----------:|---------:|---------:|------------:|
| `samtools faidx .fa`    | 13.7 ± 0.7 |     11.6 |     16.8 | 1.16 ± 0.09 |
| `samtools faidx .fa.gz` | 16.2 ± 0.9 |     14.4 |     19.1 | 1.38 ± 0.11 |
| `hnsm range .fa`        | 11.8 ± 0.7 |     10.3 |     13.9 |        1.00 |
| `hnsm range .fa.gz`     | 15.7 ± 0.8 |     14.2 |     18.8 | 1.33 ± 0.10 |
