> **Language:** [English (canonical)](QUICKSTART.md) · 中文（本文档）

# PhyloDater 快速入门指南

> **最后更新**：2026-07-12

本指南将帮助你在 5 分钟内运行第一次 PhyloDater 分析。

## 1. 安装（1 分钟）

### 使用 pip

```bash
pip install phylodater
```

### 安装外部软件（使用 conda）

```bash
conda create -n phylodater python=3.12
conda activate phylodater
conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install phylodater
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git     # r8s 后端（仅 NPRS）
pip install git+https://github.com/uym2/MD-Cat.git          # MD-Cat
# treePL 与 PATHd8 必须源码构建：bioconda/conda-forge/PyPI 已不再发布。
# 见 docs/EXTERNAL_TOOLS_CN.md §3-§4。
```

> `conda install -c bioconda r8s treepl pathd8 pyr8s` 在今天**不能工作**：这四个包已
> 不在该渠道上。请不要把这一行拷进流水线或 Dockerfile。

### 验证安装

```bash
phylodater check
```

## 2. 准备输入文件（2 分钟）

### 示例文件

你可以使用 `examples/` 目录中的示例数据：

```
examples/
├── tree.nwk           # 系统发育树
├── alignment.fasta    # 序列比对
└── calibrations.yaml  # 校准点配置
```

> **输入树必须为有根树**：校验阶段会直接拒绝无根树。运行前请先完成定根（例如在合适的外类群上定根）。

### 最小配置文件

**calibrations.yaml** 示例：

```yaml
calibrations:
  - name: "Primates"
    mrca_pair: ["human", "cat"]
    constraint:
      type: "uniform"
      min: 70.0
      max: 90.0

  - name: "Rodents"
    mrca_pair: ["mouse", "rat"]
    constraint:
      type: "soft_lower"
      min: 60.0

  - name: "Root"
    is_root: true
    constraint:
      type: "maximum"
      max: 150.0
```

> **注意**：使用 `mrca_pair` 通过两个叶节点名称指定 MRCA，或使用 `is_root: true` 指定根节点校准。

## 3. 运行分析（1 分钟）

### 基本命令（使用子命令架构）

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method pathd8
```

### 使用 MCMCTree 并检验收敛性

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method mcmctree \
    --method-args "num_runs=4"
```

### 定根选项

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method pathd8 \
    --reroot midpoint
```

### 使用 PAML 4.8

```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o my_results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.8,paml_version=4.8"
```

### 指定软件路径

PAML 版本和二进制文件路径通过 `--method-args` 传递，而非作为顶层 CLI 选项。

由于 PAML 4.8 与 4.10.8 的二进制默认同名（`mcmctree`、`codeml`、`baseml`），推荐按版本分目录安装，并指定 `paml_path`：

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.10.8,paml_version=4.10.8"
```

PhyloDater 会在运行前解析 `mcmctree` 的 banner 版本，并与 `paml_version` 比对**发布系列**：探测到 4.10.10 可满足默认值 `4.10.8`（同一 4.10 系列内控制文件语义一致），而 4.8 二进制对 4.10 期望则会被明确拒绝——两者的控制文件确实不同。指定 `paml_path` 后，PhyloDater 还会把 `<paml_path>/bin` 临时前置到 `PATH`，确保 `mcmctree` 内部调用的 `baseml` 与 `codeml` 和 `mcmctree` 属于同一版本。若 `PATH` 中已有同名命令，可通过 `mcmctree_bin=/opt/paml-4.10.8/bin/mcmctree` 强制指定完整路径。

## 4. 查看结果（1 分钟）

运行完成后，检查输出目录：

```bash
ls my_results/
```

主要输出文件：

- `summary.txt` - 执行摘要
- `comparison_table.tsv` - 节点年龄比较表
- `comparison_table.txt` - 文本格式比较表
- `pathd8/timetree.nex` - 时间树（PATHd8 NEXUS 格式）
- `mcmctree/run1/FigTree.tre` - 时间树（MCMCTree）
- `runtime_metadata.json` - 运行时元数据

## 常用命令示例

### 运行多种方法比较

```bash
# 分别运行不同方法
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method lsd2
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method treepl
phylodater dating -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ --method pathd8
```

> Each method writes to its own subdirectory under `results/`.

### 指定替换模型（LSD2）

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method lsd2 \
    --method-args "model=WAG+G"
```

### 指定初始平滑（treePL）

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method treepl \
    --method-args "initial_smooth=100"
```

### 使用配置文件

创建 `config.yaml`：

```yaml
software_paths:
  paml_path: /path/to/paml
  iqtree_bin: /path/to/iqtree2

software:
  mcmctree:
    clock: 2
    num_runs: 2
    burnin: 10000
    nsample: 20000

  lsd2:
    model: GTR+G   # 不设则由比对类型决定（核酸 GTR+G，蛋白 LG+G）
```

运行：

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --config config.yaml \
    --method mcmctree
```

## 校准点约束类型

PhyloDater 支持多种校准类型：

| 类型 | 说明 |
|------|------|
| `fixed` | 固定年龄 |
| `uniform` | 均匀分布 `[min, max]` |
| `soft_lower` | 软下界 |
| `maximum` | 上界约束 |
| `soft_bounds` | 软边界 |
| `gamma` | Gamma 先验 |
| `skew_normal` | 偏正态分布 |
| `skew_t` | 偏 t 分布 |

> **注意**：MCMCTree 的 `RootAge` 只能设置为上界，不能真正固定根节点。若根节点使用 `fixed`，PhyloDater 会在执行前报错并提示改为 `uniform` 或 `maximum`，或改用其他方法。不同方法可使用不同的校准文件。

## 序列格式

PhyloDater 支持 FASTA 和 FASTQ 格式用于分析。PHYLIP、Stockholm、Clustal 格式仅可检测，不支持分析，请先转换为 FASTA。

## 下一步

- 阅读 [README_CN.md](README_CN.md) 了解完整功能
- 查看 `examples/` 目录中的示例数据
- 参阅 `docs/api/` 了解 API 使用

## 常见问题

### 为什么提示找不到外部软件？

确保外部软件已安装在 PATH 中，或使用 `--method-args`（如 `paml_path=/path/to/paml`、`iqtree_bin=/path/to/iqtree2`）等参数指定路径。

### 校准点文件是什么格式？

PhyloDater 支持 8 种校准类型：`fixed`、`uniform`、`soft_lower`、`maximum`、`soft_bounds`、`gamma`、`skew_normal` 和 `skew_t`。使用 `uniform` 约束时，请同时提供 `min` 与 `max`。

### 支持哪些序列格式？

PhyloDater 支持 FASTA 和 FASTQ 格式用于分析。PHYLIP、Stockholm、Clustal 格式仅可检测，建议转换为 FASTA。

### 如何检验收敛性？

使用 `--method-args "num_runs=4"` 参数运行多次 MCMCTree，软件会自动计算 PSRF 和 ESS，用于收敛性诊断。

### 如何定根？

使用 `--reroot` 参数指定定根策略之一：`outgroup`、`midpoint`、`mad`、`none`。

### 如何指定 PAML 版本？

使用 `--method-args "paml_version=4.8"` 指定 PAML 版本。PhyloDater 据此选择控制文件格式和速率估算流程。
