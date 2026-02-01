# `hnsm synt` - 共线性分析命令组使用文档

## 1. 概述 (Overview)

`hnsm synt` 子命令组用于统一管理共线性分析相关的功能。

命令结构如下：

*   `hnsm synt mmg`: 基于 Minimizer 图的 DNA 共线性分析。
*   `hnsm synt das`: 基于动态规划的结构域架构相似性分析。
*   `hnsm synt merge`: 通用的共线性块合并/缝合工具。
*   `hnsm synt view`: 共线性结果可视化工具 (SVG)。

## 2. `hnsm synt mmg` 设计详情

`hnsm synt mmg` 是基于 Minimizer 图的高效 DNA 共线性分析工具。

请参考命令行帮助获取详细的算法说明和参数列表：

```bash
hnsm synt mmg --help
```

### 2.1. 结果解读示例 (Output Interpretation Example)

由于采用了分层迭代策略，不同轮次（Round）发现的共线性块可能会在坐标上相互重叠。
*   大窗口（如 1000）发现的块通常代表宏观的共线性骨架。
*   小窗口（如 10）发现的块则填充了细节，或者延伸到了大窗口未能覆盖的边缘区域。
*   **Masking 机制**：每一轮发现的区域会被“掩盖”（Masking），在后续小窗口轮次中，只会对未被掩盖的区域（即之前未能确定共线性的复杂或边缘区域）进行重新扫描和建图。这保证了算法能专注于那些难以处理的区域，而不是重复计算已经明确的核心骨架。
*   用户后续可以使用 `hnsm synt merge` 对这些重叠或碎片化的块进行合并与精简。

样例：

* `hnsm synt dna tests/genome/small_1.fa tests/genome/small_2.fa -k 21`

```tsv
# Block_ID	Range	Count	Round
# 第一轮 (Round 1001) 发现的核心区域
0	small_1.seq1(+):500-1089	589	1001
0	small_2.seq1(+):500-1089	589	1001
# 第二轮 (Round 101) 扩展到了更大的范围
1	small_1.seq1(+):50-1539	899	101
1	small_2.seq1(+):50-1539	899	101
# 第三轮 (Round 11) 进一步延伸至边缘
# 注意：Count (89) 较小是因为大部分锚点已被前两轮 Masking，这里仅包含新发现的边缘或填补空隙的锚点
# 为何看起来是连续的？
# 尽管中间区域已被 Mask（无锚点），但只要两侧剩余的 Minimizer 距离小于 `--chain-gap`，
# 它们仍会被连接在一起，从而形成跨越 Mask 区域的连续 Block。
2	small_1.seq1(+):5-1584	89	11
2	small_2.seq1(+):5-1584	89	11
```

* `hnsm synt dna tests/genome/small_1.fa tests/genome/small_2.fa -k 21 --rounds 500,10 --chain-gap 500`

```tsv
# Block_ID      Range   Count   Round
0       small_1.seq1(+):250-1339        1089    501
0       small_2.seq1(+):250-1339        1089    501
1       small_1.seq1(+):5-249   244     11
1       small_2.seq1(+):5-249   244     11
2       small_1.seq1(+):1340-1584       244     11
2       small_2.seq1(+):1340-1584       244     11
```

## 3. `hnsm synt merge` 设计详情

`hnsm synt merge` 旨在作为一个通用的、独立的后处理工具，用于将初步的共线性块（Synteny Blocks）“缝合”成更长、更完整的共线性链（Chains）。

### 3.1. 设计动机 (Motivation)

#### 3.1.1. 解决边缘碎片化问题
在使用分层迭代策略（如 `hnsm synt dna` 的 rounds 机制）时，真实的共线性区域往往因为不同分辨率的窗口采样而被切割成多个碎片。
*   **现象**：一个大块的主体在 Round 1 被发现，但边缘部分因为 minimizer 竞争而在 Round 2 或 3 被独立发现。
*   **需求**：需要一个后处理步骤，基于物理距离和线性关系，将这些碎片重新合并。

#### 3.1.2. 参数调优灵活性
共线性构建通常需要尝试不同的距离阈值（Gap Size）。将其独立后，用户可以在不重新运行昂贵的图构建过程的情况下，快速调整参数并查看结果。

### 3.2. 核心功能 (Core Features)

#### 3.2.1. 输入格式 (Input)
支持 `hnsm synt dna` 输出的标准 TSV 格式：
*   Block_ID
*   Range (Genome.Chr(Strand):Start-End)
*   Count
*   Round

#### 3.2.2. 算法逻辑 (Algorithm)
1.  **预处理**：读取输入块，解析 Range 信息，按 Block ID 聚合。
2.  **简单合并 (Simple Merge)**：
    *   逻辑：如果两个 Block 属于同一对序列、同一方向，且它们在所有基因组上的距离都小于 `chain-gap`，则进行合并。
    *   **动态默认值**：支持通过 `--divergence` 自动设置合理的 `chain-gap`。
        *   `< 1.0%`: 10000 bp
        *   `1.0 - 10.0%`: 100000 bp
        *   `> 10.0%`: 1000000 bp
    *   适用：解决同一线性路径上的碎片化问题。

#### 3.2.3. 输出 (Output)
*   格式与输入一致（TSV）。
*   合并后的块，其 Count 为子块 Count 之和。
*   Range 更新为合并后的最大范围。

### 3.3. 接口设计 (CLI Design)

```bash
hnsm synt merge [OPTIONS] <infile>
```

#### 参数选项
*   `--divergence, -d <FLOAT>`: 预估序列差异度（%）。用于自动设置默认的 `chain-gap`。
*   `--chain-gap <INT>`: 允许合并的最大间隙（bp）。如果指定，覆盖 `--divergence` 推导的默认值。
*   `--outfile, -o <FILE>`: 输出文件路径。

### 3.4. 工作流示例 (Workflows)

```bash
# 1. 生成初步块（可能包含碎片）
hnsm synt dna genome1.fa genome2.fa -o raw_blocks.tsv

# 2. 缝合碎片
# 自动根据 5% 差异度设置参数
hnsm synt merge raw_blocks.tsv -d 5 -o merged.tsv

# 或手动指定 Gap
hnsm synt merge raw_blocks.tsv --chain-gap 50000 -o merged_manual.tsv
```

## 5. `hnsm synt das` 设计详情

`hnsm synt das` 用于计算结构域架构相似性（Domain Architecture Similarity），采用动态规划算法对序列特征进行比对。

### 5.1. 参数说明

*   `infile` (Required): 输入文件。
*   `--ma <FLOAT>`: 匹配得分 (Match score)，默认 1.0。
*   `--mm <FLOAT>`: 错配得分 (Mismatch score)，默认 -0.2。
*   `--gp <FLOAT>`: 空隙惩罚 (Gap penalty)，默认 -0.01。
*   `--sep <STR>`: 分隔符，默认 "\t"。
*   `--header, -H`: 标记输入文件是否包含表头。
*   `--outfile, -o <FILE>`: 输出文件路径 (默认: stdout)。

## FastGA + UCSC chainNet 共线性检测

```bash
FastGA -v -psl tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz > tmp.psl

pgr chain -t="" tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz tmp.psl > tmp.syn.maf


```

## 基于蛋白相似度的共线性检测

### E. coli mg1655 基因组内部

```bash
# 生成基因位置
hnsm gff rg tests/genome/mg1655.gff.gz --tag cds --key protein_id --asm mg1655 -s --ss \
    -o tests/genome/mg1655.rg.tsv

# 计算蛋白相似度矩阵
hnsm dist seq tests/genome/mg1655.pro.fa.gz -k 7 -w 2 --sim -p 8 |
    rgr filter stdin --ff-ne 1:2 \
    > tests/genome/mg1655.sim.tsv

# dag chain
# 接收两个文件：1. 基因位置文件 2. 蛋白相似度矩阵文件
hnsm synt dag mg1655.rg.tsv mg1655.sim.tsv -o tests/genome/mg1655.dag.tsv

# 显示共线性块
cargo run --bin hnsm synt ribbon tests/genome/mg1655.dag.tsv tests/genome/mg1655.size.tsv -o mg1655.ribbon.svg

cargo run --bin hnsm synt circle tests/genome/mg1655.dag.tsv tests/genome/mg1655.size.tsv -o mg1655.circle.svg


```

## 6. 背景与 ntSynt 对比 (Background & Comparison)

`hnsm synt dna` 的设计深受 **ntSynt** (Coombe et al., 2025) 的启发。我们的目标是在 Rust 生态中重现并优化其基于 Minimizer 图的宏观共线性检测算法，提供更高效、更易用的单一二进制工具。

### 6.1. 核心相似点 (Similarities)

*   **算法逻辑**: 两者都采用 "Sketching -> Graph Building -> Path Finding" 的流程。
*   **分层迭代**: 都使用从大窗口到小窗口的迭代策略（Iterative Refinement），并结合 Global Masking 技术来逐步解析复杂区域。
*   **内存优化**: 都利用 Bloom Filter 来过滤 Singleton k-mers，显著降低大基因组分析的内存消耗。
*   **参数模型**: 默认参数（如 block size, chain gap）都基于序列差异度（Divergence）进行动态调整。

### 6.2. 主要差异与改进 (Differences & Improvements)

| 特性 | ntSynt (Original) | hnsm synt dna (Rust Implementation) |
| :--- | :--- | :--- |
| **语言与架构** | Python 脚本 + C++ 二进制 (Hybrid Pipeline) | 纯 Rust 实现 (Single Binary)，无外部依赖 |
| **图数据结构** | `igraph` (Python/C binding) | `petgraph` (Rust)，配合定制的内存优化邻接表 |
| **哈希算法** | `ntHash` | `RapidHash` / `MurmurHash3` |
| **软掩膜支持** | 需预处理 | 内置 `--soft-mask`，动态过滤小写碱基 k-mer |
| **后处理** | 脚本处理 | 独立的 `hnsm synt merge` 命令，支持灵活缝合 |
| **可视化** | 依赖 Python/R (ntSynt-viz) | 内置 `hnsm synt view` 命令，直接生成 SVG |
| **输入/输出** | 依赖文件系统中间文件 | 流式处理，支持标准输入/输出 |
| **安装与部署** | 需配置 Python 环境和编译 C++ | `cargo install` 或单文件分发 |

### 6.3. 演进过程

`hnsm` 逐步移植并重构了 ntSynt 的关键组件：
1.  **哈希与索引**: 使用 `minimizer_iter` 和 `rapidhash` 替代了 `indexlr`。
2.  **图算法**: 使用 Rust 的 `petgraph` 库重写了构图和路径查找逻辑，并针对共线性图的稀疏特性进行了 O(E) 遍历优化。
3.  **工程化**: 将原本分散的预处理、过滤、比对步骤整合为单一命令 `hnsm synt dna`，并统一了 I/O 格式。
