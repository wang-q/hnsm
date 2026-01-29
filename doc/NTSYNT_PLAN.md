# hnsm 的 ntSynt 功能实现计划

## 1. 概述
本文档概述了将 **ntSynt**（使用动态 minimizer 图的多基因组宏观共线性检测）的核心功能实现到 `hnsm` 中的计划。目标是在 `hnsm` 工具包内提供高性能、纯 Rust 实现的共线性检测，以替代原始 ntSynt 的混合 Python/C++ 流程。

## 2. 功能差距分析

| 功能 | ntSynt (原版) | hnsm (当前) | 所需行动 |
| :--- | :--- | :--- | :--- |
| **差异度估算** | `Mash` (外部工具) | `hnsm distance` (Minimizer Jaccard/Mash) | **就绪**。使用 `hnsm distance` 估算差异度并设置参数。 |
| **序列过滤** | Bloom Filter (公共/重复 k-mers) | `hnsm::hash` (HashSets) | **完成**。实现了 Bloom Filter 优化内存使用，用于过滤单次出现的 k-mers。 |
| **Minimizer 生成** | `indexlr` (btllib) | `hnsm::hash` | **已完成**。基于 `minimizer_iter` 实现了 `seq_sketch`，支持位置、链信息及多种哈希算法。 |
| **图构建** | Minimizer Graph (igraph) | `hnsm::synteny::graph` | **已完成**。基于 `petgraph` 实现了 `SyntenyGraph`，支持多基因组图构建、边权重过滤及路径查找。 |
| **共线性检测** | 路径查找 / 组件分析 | `hnsm::loc` | **已完成**。基于图遍历算法，寻找共线性路径并转换为 SyntenyBlock。 |
| **细化 (Refinement)** | 迭代减小 `w` (窗口大小) | `hnsm::io`, `hnsm::nt` | **已完成**。实现了多轮迭代细化（rounds），利用 IntSpan 掩盖已覆盖区域。 |

## 3. 已实现的库功能 (`src/libs`)
这些功能已在 `hnsm` 内部库中实现，可直接用于 `ntSynt` 的开发：

*   **`hnsm::hash`**:
    *   支持 `FxHash`, `MurmurHash3`, `RapidHash` 三种快速哈希算法。
    *   实现了 `JumpingMinimizer` 算法，用于生成序列的最小哈希值。
*   **`hnsm::loc`**:
    *   高效的 FASTA 序列索引 (`.loc` 文件)。
    *   支持对大文件进行快速随机访问和序列长度提取。
*   **`hnsm::linalg` & `hnsm::mds`**:
    *   `linalg`: 提供 SIMD 优化的欧几里得距离、点积等线性代数运算。
    *   `mds`: 实现主坐标分析 (PCoA)，用于多基因组关系的降维可视化。
*   **`hnsm::io` & `hnsm::nt`**:
    *   `io`: 处理 `AsmEntry` 格式（名称+向量）及序列格式检测。
    *   `nt`: 提供核苷酸 (A/C/G/T/N) 的高效编码和字符映射。

## 4. 现有的命令行工具 (`src/cmd`)
`hnsm` 已包含多个强大的命令行工具，可辅助 `ntSynt` 流程或作为其组件：

*   **序列处理**:
    *   `count`: 计算碱基统计，用于基因组概览。
    *   `filter`: 序列过滤（长度、N含量、去重），对应 `ntSynt` 的预处理步骤。
    *   `order`: 提取指定序列，用于提取共线性块序列。
    *   `dedup`: 高级去重，确保输入数据的质量。
*   **高级分析**:
    *   `prefilter`: 氨基酸 minimizer 预过滤，提供了一种快速筛选同源序列的方法。
    *   `das`: 结构域架构相似性，为基于功能的共线性提供思路。
    *   `hv` & `similarity`: 向量化序列比较，提供了除 minimizer 图之外的另一种全局相似性视角。

## 5. 实施路线图

### 第一阶段：增强的 Minimizer 索引
*   **状态**：已实现 (`src/libs/hash.rs`)
*   **目标**：生成包含完整位置信息的 minimizer。
*   **已完成任务**：
    1.  将功能集成到 `hnsm::hash` 模块。
    2.  定义了 `MinimizerInfo` 结构体，包含 `hash`, `seq_id`, `pos`, `strand`。
    3.  实现了 `seq_sketch` 函数，支持规范化 k-mer 哈希（Canonical Hashing）和滑动窗口最小值（利用 `minimizer_iter` crate）。
    4.  添加了 k-mer 过滤回调接口，支持后续的公共/重复 k-mer 过滤。

### 第二阶段：动态 Minimizer 图
*   **状态**：已实现 (`src/libs/synteny/graph.rs`)
*   **目标**：构建用于共线性检测的图结构。
*   **已完成任务**：
    1.  创建了 `hnsm::synteny::graph` 模块。
    2.  定义了 `Node` (Minimizer), `Edge` (邻接关系), `Occurrence` (基因组位置) 结构。
    3.  实现了 `SyntenyGraph` 结构及多基因组图构建逻辑（`add_minimizers`）。
    4.  实现了边权重过滤 (`prune_low_weight_edges`) 和线性路径查找 (`get_linear_paths`)。
    5.  实现了基于频率的 k-mer 过滤（在 `hnsm synteny` 命令中通过两遍扫描实现）。
    6.  针对超大基因组优化过滤内存（引入 Bloom Filter）。

### 第三阶段：共线性块识别
*   **状态**：已完成 (`src/libs/synteny/algo.rs`)
*   **目标**：从图中提取共线性块。
*   **已完成任务**：
    1.  整合流程：在 `SyntenyFinder::run` 中实现了控制逻辑。
    2.  块处理：实现了 `SyntenyBlock::from_path`。
    3.  输出：CLI 输出了 TSV 格式。

### 第四阶段：迭代细化
*   **状态**：已完成 (`src/libs/synteny/algo.rs`)
*   **目标**：细化块边界并填充空隙。
*   **已完成任务**：
    1.  实现了基于 `rounds` 参数的迭代循环。
    2.  利用 `intspan::IntSpan` 记录和更新 `coverage`。
    3.  在后续轮次中仅处理未覆盖区域。
    4.  集成 Bloom Filter 优化内存。

### 第五阶段：CLI 集成
*   **状态**：已完成 (`src/cmd/synteny.rs`)
*   **目标**：用户友好的命令 `hnsm synteny`。
*   **已完成任务**：
    1.  在 `src/cmd/mod.rs` 中添加 `synteny` 子命令。
    2.  参数：
        *   `infiles`：FASTA 文件列表。
        *   `-k`：k-mer 大小。
        *   `--rounds`：窗口大小列表（如 "100,10"）。
        *   `--min-weight`：最小边权重。
        *   `--max-freq`：最大 k-mer 频率。
    3.  输出：TSV 格式输出。

## 5. 验证与测试
*   **单元测试**:
    *   `src/libs/synteny/tests.rs`: 涵盖图构建、路径查找、循环处理和 Bloom Filter 集成。
    *   `src/libs/bloom.rs`: Bloom Filter 基本功能测试。
*   **功能测试**:
    *   使用 `tests/fasta` 中的数据进行端到端测试。
    *   验证了 `rounds` 迭代和内存优化效果。

## 6. 参考资料
*   **ntSynt Paper**: Coombe et al., 2025.
*   **Codebase**: `c:\Users\wangq\Scripts\hnsm\ntSynt-main`
*   **Existing hnsm modules**: `src/cmd/distance.rs`.
