# `hnsm synt` - 共线性分析命令组设计方案

## 1. 概述 (Overview)

`hnsm synt` 子命令组用于统一管理共线性分析相关的功能。

命令结构如下：

*   `hnsm synt dna`: 基于 Minimizer 图的 DNA 共线性分析。
*   `hnsm synt das`: 基于动态规划的结构域架构相似性分析。
*   `hnsm synt merge`: 通用的共线性块合并/缝合工具。

## 2. `hnsm synt dna` 设计详情

`hnsm synt dna` 是基于 Minimizer 图的高效 DNA 共线性分析工具，专为处理全基因组序列设计。

### 2.1. 核心算法 (Algorithm)

该工具采用分层迭代（Hierarchical Iterative）策略，结合 Minimizer 图与线性路径查找算法，从粗到细逐步构建共线性块。

1.  **分层迭代 (Iterative Refinement)**:
    *   采用多轮（Rounds）策略，窗口大小（Window Size, `w`）从大到小变化（例如 1000 -> 100 -> 10）。
    *   **全局掩膜 (Global Masking)**: 每一轮发现的共线性区域会被记录在全局掩膜（Coverage Mask）中。后续轮次在提取 Minimizer 时，会跳过已覆盖的区域，从而专注于未解决的复杂区域或更精细的结构。

2.  **Minimizer 提取与过滤 (Pass 1 & 2)**:
    *   **Pass 1 (Counting)**: 遍历所有序列，统计 Minimizer 频率。使用 Bloom Filter 过滤掉仅出现一次的 Minimizer（Singletons），以减少内存占用。
    *   **Pass 2 (Graph Building)**: 再次遍历序列，构建 Minimizer 图。
        *   **Repeat Filtering**: 忽略频率超过 `max-freq` 的高频 Minimizer（通常是转座子或简单重复序列）。
        *   **Mask Filtering**: 忽略落在全局掩膜区域内的 Minimizer。

3.  **图构建与简化 (Graph Processing)**:
    *   **节点 (Nodes)**: 代表唯一的 Minimizer Hash。
    *   **边 (Edges)**: 代表两个 Minimizer 在基因组上的邻接关系。
    *   **权重过滤 (Weight Pruning)**: 仅保留权重（支持的基因组数量）大于 `min-weight` 的边，去除偶然的噪音连接。
    *   **传递归约 (Transitive Reduction)**: 移除冗余的传递边（如 A->B->C 存在时，移除 A->C），简化图结构。

4.  **线性路径查找 (Linear Path Finding)**:
    *   使用 O(E) 的全图扫描算法。
    *   识别图中互为唯一邻居（Reciprocal Best Hits in Graph）的节点链。
    *   这些链代表了多基因组间高度保守的线性共线性骨架。

5.  **块构建 (Block Construction)**:
    *   将抽象的 Hash 路径映射回具体的基因组坐标。
    *   使用二分查找（Binary Search）高效定位每个 Hash 在各基因组中的具体位置（Range）。
    *   输出包含所有输入基因组对应坐标的共线性块。

### 2.2. 参数说明 (Parameters)

*   `--divergence, -d <FLOAT>`: 预估的序列差异度（%）。此参数会自动调整 rounds, block-size, merge 等默认值。
*   `--kmer, -k <INT>`: Minimizer k-mer 大小（默认 24）。
*   `--max-freq <INT>`: 过滤高频 k-mer 的阈值（默认 1000）。
*   `--rounds, -r <STR>`: 自定义分层迭代的窗口大小序列（如 "10000,1000,100"）。
*   `--block-size, -b <INT>`: 最小共线性块大小（bp）。
*   `--min-weight <INT>`: 最小边权重（支持的基因组数量，默认 2）。

### 2.3. 输出格式 (Output)

输出为 Block TSV 格式。每一行代表一个基因组片段（Range），共享相同 ID 的行属于同一个共线性块（Block）。

字段说明：
1.  **Block_ID**: 共线性块的唯一标识符（整数）。
2.  **Range**: 格式为 `Genome.Chr(Strand):Start-End`。
    *   `Genome.Chr`: 序列名称，通常由文件名和染色体名拼接而成（如 `S288c.I`）。
    *   `Strand`: `+` 或 `-`。注意：第一个基因组的链方向会被标准化为 `+`。
    *   `Start-End`: 1-based 起止坐标。
3.  **Count**: 块中包含的 Minimizer 数量。
4.  **Round**: （可选）发现该块时的窗口大小（迭代轮次）。

**注意**:
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

`hnsm synt merge` 旨在作为一个通用的、独立的后处理工具，用于将碎片化的同源锚点（Anchors）或初步共线性块（Synteny Blocks）“缝合”成更长、更完整的共线性链（Chains）。

该设计遵循 "Do one thing and do it well" 的 Unix 哲学，将**特征发现**（Feature Detection）与**共线性构建**（Synteny Chaining）解耦。

### 3.1. 设计动机 (Motivation)

#### 3.1.1. 解决边缘碎片化问题
在使用分层迭代策略（如 `hnsm synt dna` 的 rounds 机制）时，真实的共线性区域往往因为不同分辨率的窗口采样而被切割成多个碎片。
*   **现象**：一个大块的主体在 Round 1 被发现，但边缘部分因为 minimizer 竞争而在 Round 2 或 3 被独立发现。
*   **需求**：需要一个后处理步骤，基于物理距离和线性关系，将这些碎片重新合并。

#### 3.1.2. 支持多种数据源 (Universality)
共线性分析不应局限于 DNA minimizer。通过标准化输入接口，`merge` 命令可以支持：
*   **DNA Synteny**: `hnsm synt dna` 的输出。
*   **Protein Synteny**: 基于 BLASTP, Diamond, MMseqs2 等比对结果。
*   **Gene Orthology**: 基于基因位置的同源性分析（如 MCScanX 流程）。
*   **Generic Anchors**: 任何具有坐标信息的同源片段（GFF, BED）。

#### 3.1.3. 参数调优灵活性
共线性构建通常需要尝试不同的距离阈值（Gap Size）和惩罚参数。将其独立后，用户可以在不重新运行昂贵的序列比对/图构建过程的情况下，快速调整参数并查看结果。

### 3.2. 核心功能 (Core Features)

#### 3.2.1. 输入格式 (Input)
支持标准的 TSV 格式（类似 BLAST m8 或 BEDPE），至少包含以下字段：
*   Query ID (`q_name`)
*   Query Start (`q_start`)
*   Query End (`q_end`)
*   Target ID (`t_name`)
*   Target Start (`t_start`)
*   Target End (`t_end`)
*   Strand (`strand`: +/-)
*   Score (`score`: 可选，用于加权链)

#### 3.2.2. 算法逻辑 (Algorithm)
1.  **预处理**：读取输入锚点，按 Query/Target ID 分组并排序。
2.  **简单合并 (Simple Merge)**：
    *   针对 `hnsm synt dna` 的主要场景。
    *   逻辑：如果两个 Block 属于同一对序列、同一方向，且它们在 Query 和 Target 上的距离都小于 `max_gap`，则直接合并。
    *   适用：解决同一线性路径上的碎片化问题。
3.  **动态规划链 (DAG Chaining)** (进阶功能)：
    *   类似 DAGChainer 或 MCScanX。
    *   构建有向无环图，寻找得分最高的路径。
    *   适用：处理复杂的重排、过滤噪音、从散乱的蛋白比对中构建共线性。

#### 3.2.3. 输出 (Output)
*   **Chained Blocks**: 合并后的超级块，格式与输入类似，但坐标范围更大。
*   **Statistics**: 报告合并前后的块数量、覆盖度变化。

### 3.3. 接口设计 (CLI Design)

```bash
hnsm synt merge [OPTIONS] <INPUT>
```

#### 参数选项
*   `--max-gap <INT>`: 允许合并的最大间隙（bp），默认如 100000。
*   `--input-fmt <STR>`: 输入格式预设 (auto, synteny, blast, gff)。
*   `--algo <STR>`: 算法选择
    *   `simple`: 仅基于距离的贪婪合并（默认，速度快）。
    *   `dag`: 基于动态规划的全局最优链（适合高噪音或蛋白比对）。

### 3.4. 工作流示例 (Workflows)

#### 场景 A: 优化 `hnsm synt dna` 结果
```bash
# 1. 生成初步块（可能包含碎片）
hnsm synt dna genome1.fa genome2.fa -o raw_blocks.tsv

# 2. 缝合碎片
hnsm synt merge raw_blocks.tsv --max-gap 50000 -o final_synteny.tsv
```

#### 场景 B: 蛋白共线性 (未来扩展)
```bash
# 1. 蛋白比对
diamond blastp -q proteome1.fa -d proteome2 -f 6 ... > hits.tsv

# 2. 转换为坐标锚点 (假设有工具或脚本提取基因位置)
# ... hits.tsv + gff -> anchors.tsv ...

# 3. 构建共线性
hnsm synt merge anchors.tsv --algo dag -o gene_synteny.tsv
```

## 4. 开发计划

1.  **重构 CLI 结构** (已完成):
    *   创建 `src/cmd/synt/` 模块。
    *   将 `src/cmd/synteny.rs` 移动至 `src/cmd/synt/dna.rs`。
    *   将 `src/cmd/das.rs` 移动至 `src/cmd/synt/das.rs`。
    *   将 `src/cmd/chain.rs` 重构为 `src/cmd/synt/merge.rs`。
    *   更新 `src/hnsm.rs` 注册新命令组。

2.  **实现 `merge` 功能**:
    *   定义 `Anchor` 和 `Chain` 结构体。
    *   实现 `simple` 合并算法。
    *   集成到 `src/cmd/synt/merge.rs`。

3.  **集成测试**:
    *   使用 `hnsm synt dna` 的输出作为测试数据，验证能否正确合并人为分割的块。
