# `hnsm synt` - 共线性分析命令组设计方案

## 1. 概述 (Overview)

为了统一管理共线性分析相关的功能，我们将原有的 `synteny` 命令和新设计的 `chain` 命令整合到 `hnsm synt` 子命令组下。

新的命令结构如下：

*   `hnsm synt dna`: (原 `hnsm synteny`) 基于 Minimizer 图的 DNA 共线性分析。
*   `hnsm synt merge`: (原 `hnsm chain`) 通用的共线性块合并/缝合工具。

## 2. `hnsm synt merge` 设计详情

`hnsm synt merge` 旨在作为一个通用的、独立的后处理工具，用于将碎片化的同源锚点（Anchors）或初步共线性块（Synteny Blocks）“缝合”成更长、更完整的共线性链（Chains）。

该设计遵循 "Do one thing and do it well" 的 Unix 哲学，将**特征发现**（Feature Detection）与**共线性构建**（Synteny Chaining）解耦。

### 2.1. 设计动机 (Motivation)

#### 2.1.1. 解决边缘碎片化问题
在使用分层迭代策略（如 `hnsm synt dna` 的 rounds 机制）时，真实的共线性区域往往因为不同分辨率的窗口采样而被切割成多个碎片。
*   **现象**：一个大块的主体在 Round 1 被发现，但边缘部分因为 minimizer 竞争而在 Round 2 或 3 被独立发现。
*   **需求**：需要一个后处理步骤，基于物理距离和线性关系，将这些碎片重新合并。

#### 2.1.2. 支持多种数据源 (Universality)
共线性分析不应局限于 DNA minimizer。通过标准化输入接口，`merge` 命令可以支持：
*   **DNA Synteny**: `hnsm synt dna` 的输出。
*   **Protein Synteny**: 基于 BLASTP, Diamond, MMseqs2 等比对结果。
*   **Gene Orthology**: 基于基因位置的同源性分析（如 MCScanX 流程）。
*   **Generic Anchors**: 任何具有坐标信息的同源片段（GFF, BED）。

#### 2.1.3. 参数调优灵活性
共线性构建通常需要尝试不同的距离阈值（Gap Size）和惩罚参数。将其独立后，用户可以在不重新运行昂贵的序列比对/图构建过程的情况下，快速调整参数并查看结果。

### 2.2. 核心功能 (Core Features)

#### 2.2.1. 输入格式 (Input)
支持标准的 TSV 格式（类似 BLAST m8 或 BEDPE），至少包含以下字段：
*   Query ID (`q_name`)
*   Query Start (`q_start`)
*   Query End (`q_end`)
*   Target ID (`t_name`)
*   Target Start (`t_start`)
*   Target End (`t_end`)
*   Strand (`strand`: +/-)
*   Score (`score`: 可选，用于加权链)

#### 2.2.2. 算法逻辑 (Algorithm)
1.  **预处理**：读取输入锚点，按 Query/Target ID 分组并排序。
2.  **简单合并 (Simple Merge)**：
    *   针对 `hnsm synt dna` 的主要场景。
    *   逻辑：如果两个 Block 属于同一对序列、同一方向，且它们在 Query 和 Target 上的距离都小于 `max_gap`，则直接合并。
    *   适用：解决同一线性路径上的碎片化问题。
3.  **动态规划链 (DAG Chaining)** (进阶功能)：
    *   类似 DAGChainer 或 MCScanX。
    *   构建有向无环图，寻找得分最高的路径。
    *   适用：处理复杂的重排、过滤噪音、从散乱的蛋白比对中构建共线性。

#### 2.2.3. 输出 (Output)
*   **Chained Blocks**: 合并后的超级块，格式与输入类似，但坐标范围更大。
*   **Statistics**: 报告合并前后的块数量、覆盖度变化。

### 2.3. 接口设计 (CLI Design)

```bash
hnsm synt merge [OPTIONS] <INPUT>
```

#### 参数选项
*   `--max-gap <INT>`: 允许合并的最大间隙（bp），默认如 100000。
*   `--input-fmt <STR>`: 输入格式预设 (auto, synteny, blast, gff)。
*   `--algo <STR>`: 算法选择
    *   `simple`: 仅基于距离的贪婪合并（默认，速度快）。
    *   `dag`: 基于动态规划的全局最优链（适合高噪音或蛋白比对）。

### 2.4. 工作流示例 (Workflows)

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

## 3. 开发计划

1.  **重构 CLI 结构**:
    *   创建 `src/cmd/synt/mod.rs`。
    *   将 `src/cmd/synteny.rs` 移动/重命名为 `src/cmd/synt/dna.rs`。
    *   将 `src/cmd/chain.rs` 移动/重命名为 `src/cmd/synt/merge.rs`。
    *   更新 `src/hnsm.rs` 注册新命令。

2.  **实现 `merge` 功能**:
    *   定义 `Anchor` 和 `Chain` 结构体。
    *   实现 `simple` 合并算法。
    *   集成到 `src/cmd/synt/merge.rs`。

3.  **集成测试**:
    *   使用 `hnsm synt dna` 的输出作为测试数据，验证能否正确合并人为分割的块。
