# `hnsm` and `faops`

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

```bash
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "hnsm  size ufasta.fa" \
    'cat tests/fasta/ufasta.fa | hnsm  size stdin > /dev/null' \
    -n "faops size ufasta.fa" \
    'cat tests/fasta/ufasta.fa | faops size stdin > /dev/null' \
    -n "hnsm  size ufasta.fa.gz" \
    'hnsm  size tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops size ufasta.fa.gz" \
    'faops size tests/fasta/ufasta.fa.gz > /dev/null'

```

| Command                   | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------------|----------:|---------:|---------:|------------:|
| `hnsm  size ufasta.fa`    | 2.7 ± 0.3 |      2.3 |      3.9 | 1.02 ± 0.14 |
| `faops size ufasta.fa`    | 2.7 ± 0.2 |      2.1 |      3.7 |        1.00 |
| `hnsm  size ufasta.fa.gz` | 2.8 ± 0.3 |      2.4 |      3.9 | 1.05 ± 0.14 |
| `faops size ufasta.fa.gz` | 2.9 ± 0.3 |      2.4 |      3.9 | 1.07 ± 0.14 |
