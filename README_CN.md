> **Language:** [English (canonical)](README.md) · 中文（本文档）

# PhyloDater

**多软件并行系统发育定年平台 / A Multi-Software Parallel Platform for Phylogenetic Molecular Dating**

[![Python Version](https://img.shields.io/badge/python-3.9%20--%203.14-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **版本**：0.1.0 | **许可证**：MIT License

[English](README.md)

## 概述

PhyloDater 是一个统一的系统发育定年平台，整合多种主流定年软件，提供一致的接口和标准化的结果输出。

> **注意**：本软件处理系统发育树（非基因组坐标），不涉及 0-based/1-based 坐标体系问题。

### 遵循的标准与格式

- **树文件格式**：[Newick](https://en.wikipedia.org/wiki/Newick_format)（`.nwk`）、[Nexus](https://en.wikipedia.org/wiki/Nexus_file)（`.nex`）
- **序列格式**：[FASTA](https://en.wikipedia.org/wiki/FASTA_format)、[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)
- **配置格式**：YAML 1.2

### 主要特性

- **多方法支持**：MCMCTree、LSD2、r8s、treePL、PATHd8、wLogDate、MD-Cat
- **子命令架构**：`dating`（定年分析）、`check`（系统自检）、`fold`（树折叠）
- **分类学解析**：支持 GTDB 嵌入式格式和表格分号式格式
- **深度输入验证**：括号平衡、负分支长度、重复节点名、恶意字符检测
- **实时日志**：ISO 8601 时间戳，流式输出，支持文件持久化
- **单系群判定**：自动识别 LUCA、LACA、LBCA 等特殊标识符
- **自定义软件路径**：通过配置文件 `software_paths` 指定各外部软件二进制路径

### 支持的定年方法

| 方法 | 描述 | `--method-args` 参数 |
|------|------|---------------------|
| **MCMCTree** | PAML 包中的贝叶斯定年工具 | `clock`、`num_runs`、`burnin`、`nsample`、`sampfreq`、`skip_rate_estimation`、`paml_version`、`paml_path`、`mcmctree_bin` |
| **LSD2** | IQ-TREE2 中的快速定年方法 | `model`、`date_ci`、`clock_sd`、`date_options`、`date_outlier` |
| **r8s** | 经典分子钟定年工具 | `method`、`smoothing`、`algorithm` |
| **treePL** | 基于惩罚似然的定年方法 | `initial_smooth`、`cvstart`、`cvstop`、`cviter`、`cvmultstep` |
| **PATHd8** | 分子钟校准方法 | （无额外参数） |
| **wLogDate** | 日志日期定年 | `backward_time`、`num_replicates`、`root_time`、`leaf_time`、`sequence_length`、`max_iter`、`pseudocount`、`zero_branch_len` |
| **MD-Cat** | 矩阵校准日期 | `ncat`、`nrep`、`max_iter`、`backward_time`、`as_date`、`use_direct_import`、`ci_nboots`、`ci_plower`、`ci_pupper`、`annotate_level` |

### 分析流程

```
输入文件                验证                   分析                   输出
──────────────────────────────────────────────────────────────────────────────
树文件 (.nwk)  ─┐
                ├─→ QC Processing(可选) ─→ 树定根 ─→ 输入验证 ─┐
序列文件 (.fa)  ─┤                            │           │         │
                │                            ↓           ↓         ↓
校准点 (.yaml) ─┘                       深度验证       MRCA 定位    定年方法执行
                                        括号平衡       单系群判定    MCMCTree
                                        负分支长度     特殊标识符    treePL
                                        重复节点名     分类学解析    r8s/PATHd8/...
                                        交叉验证                     │
                                                                     ↓
                                                              结果比较与报告
```

#### 方法专有参数默认值

| 方法 | 参数 | 默认值 | 说明 |
|------|------|--------|------|
| MCMCTree | `clock` | `2` | 1=严格钟，2=独立松弛钟，3=自相关松弛钟 |
| | `num_runs` | `2` | MCMC 并行运行次数（用于收敛性检验） |
| | `burnin` | `20000` | MCMC burn-in 代数 |
| | `nsample` | `50000` | MCMC 采样数 |
| | `sampfreq` | `5` | MCMC 采样频率 |
| | `skip_rate_estimation` | `false` | 是否跳过阶段 0 速率估算（`true`=跳过并使用默认 rgene_gamma 先验）；版本行为差异见“PAML 4.8 支持”一节 |
| | `paml_version` | `4.10.8` | 期望的 PAML 发布**系列**；探测到 4.10.x 可通过，4.8 不行 |
| | `paml_path` | - | PAML 安装目录（推荐），会自动使用 `<paml_path>/bin/mcmctree`，并隔离同版本的 `baseml` 与 `codeml` |
| | `mcmctree_bin` | `mcmctree` | MCMCTree 可执行文件完整路径或命令名 |
| r8s | `method` | `NPRS` | 算法：`NPRS`（默认）/ `PL`（仅原生 r8s 二进制） |
| | `smoothing` | `None` | 平滑参数（None = 不写入 NEXUS） |
| | `algorithm` | `TN` | 树搜索算法 |
| | `r8s_bin` | `r8s` | r8s 或 pyr8s 可执行文件的名称或路径 |
| treePL | `initial_smooth` | `100.0` | 初始平滑参数 |
| | `cvstart` | `1000.0` | 交叉验证起始值 |
| | `cvstop` | `0.1` | 交叉验证终止值 |
| | `cviter` | `10` | CV 迭代次数 |
| | `cvmultstep` | `0.1` | CV 乘子步长 |
| LSD2 | `model` | _自动_ | 替换模型；未指定时按比对类型选择（核酸 `GTR+G`，蛋白 `LG+G`）；可显式覆盖，如 `HKY+G` |
| | `date_ci` | `100` | 置信区间自助重复次数 |
| | `clock_sd` | `0.2` | 分子钟标准差 |
| MD-Cat | `ncat` | `50` | 速率类别数 |
| | `nrep` | `100` | 随机初始化次数 |
| | `max_iter` | `100` | 最大 EM 迭代次数 |
| | `ci_nboots` | `0` | 置信区间 bootstrap 次数；`0` = 关闭（与上游一致）。单次开销很大，请按需显式开启 |
| wLogDate | `backward_time` | `True` | 使用反向时间方向 |
| | `num_replicates` | `1` | 重复次数 |
| | `max_iter` | `50000` | 最大迭代次数 |
| | `pseudocount` | `0.01` | 分支长度伪计数 |

## 安装

### 通过 pip 安装

```bash
pip install phylodater
```

### 从源码安装

```bash
git clone https://github.com/ZengZichao/phylodater.git
cd phylodater
pip install -e .
```

### 前提条件

每个定年方法都依赖一个外部引擎，需自行安装（不在 pip 依赖中）。下表于 2026-09-29
根据在线 conda / PyPI 渠道元数据核实；没有现成包的方法会给出源码获取途径，而不是一
条注定失败的命令。

| 软件 | 方法 | 来源 | 安装方式 |
|------|------|------|----------|
| PAML（`mcmctree`、`baseml`） | MCMCTree | https://abacus.gene.ucl.ac.uk/software/paml.html | bioconda `paml` |
| IQ-TREE 2（内置 LSD2） | LSD2 | https://github.com/iqtree/iqtree2 | bioconda `iqtree` |
| wLogDate | wLogDate | https://pypi.org/project/wlogdate/ | PyPI `wlogdate` 或 bioconda `wlogdate` |
| pyr8s | r8s | https://github.com/iTaxoTools/pyr8s | 源码（`pip install git+…`）——不在 PyPI/conda |
| r8s（原版） | r8s | https://sourceforge.net/projects/r8s/ | 源码构建——不在 PyPI/conda |
| treePL | treePL | https://github.com/blackrim/treePL | **源码构建**——bioconda/conda-forge/PyPI 已不再发布 |
| PATHd8 | PATHd8 | https://www2.math.su.se/PATHd8/ | **源码构建**（`cc PATHd8.c -lm -O3`）——conda/PyPI 已不再发布 |
| MD-Cat | MD-Cat | https://github.com/uym2/MD-Cat | 源码（`pip install git+…`）——不在 PyPI/conda |

各引擎在发布测试中实测到的版本、校验和与具体构建命令：
[docs/EXTERNAL_TOOLS_CN.md](docs/EXTERNAL_TOOLS_CN.md)。

> **r8s 与 pyr8s 说明**：若 `pyr8s` 可导入，PhyloDater 会优先通过其原生 API 调用——
> 它**只实现 NPRS**，因此配置里的 `method: PL` 会被忽略并告警；且其 `CONSTRAIN` 区间
> 只是参考值而不会被强制执行（实测，见
> [docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md) §2.2）。若没有 `pyr8s`，则回退到
> 由 `r8s_bin` 指定的命令行 `r8s` 二进制。

### 使用 Conda

```bash
conda create -n phylodater-env python=3.12
conda activate phylodater-env
conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install phylodater
# 其余引擎（已核实的路径，见 docs/EXTERNAL_TOOLS_CN.md）：
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git
pip install git+https://github.com/uym2/MD-Cat.git
# treePL 与 PATHd8：源码构建（docs/EXTERNAL_TOOLS_CN.md §3-§4）
```

### 使用 Docker

```bash
# 构建镜像
docker build -t phylodater .

# 运行容器
docker run -v $(pwd):/data phylodater dating \
    -t /data/tree.nwk \
    -s /data/alignment.fasta \
    -c /data/calibrations.yaml \
    -o /data/results/ \
    --method mcmctree
```

> **镜像到底含有什么。** 镜像只构建了 **7 个引擎中的 3 个** —— MCMCTree（PAML）、
> LSD2（IQ-TREE 2）与 PATHd8 —— 因此不能把它说成"支持全部 7 种方法"的平台。
> 逐引擎的覆盖表见 [docs/EXTERNAL_TOOLS_CN.md](docs/EXTERNAL_TOOLS_CN.md) §7。
> 请求一个引擎缺失的方法**不会静默跳过**，也**不会编造结果**：运行会依次打出
> `<engine> not found`（含安装建议）、`... environment not available, skipping`、
> `N methods failed: <method> (Environment not available)` 与 `No results to report`。
> （2026-09-30 实测：在干净 sdist 安装、`PATH` 上没有任何引擎的情况下跑 CLI；
> 与 `phylodater check` 可选依赖行所报告的行为一致。）本轮验证**没有执行**
> `docker build` —— 验证机上没有 Docker/Podman 运行时 —— 所以镜像应按"此处未验证"
> 对待，而不是"已测试"。

## 快速开始

### 1. 自检（验证安装）

```bash
phylodater check
```

输出示例（为发布环境上 2026-09-30 的原样实测输出）：
```
检查项                                      状态       信息
----------------------------------------------------------------------
依赖: biopython                            [PASS]   v1.88
可选依赖: logdate (wLogDate)                 [PASS]   vunknown
可选依赖: emd (MD-Cat)                       [PASS]   vunknown
可选依赖: pyr8s                              [PASS]   vunknown
可选依赖: ete3                               [PASS]   v3.1.3
...
总计: 26 项检查, 26 通过, 0 失败

✓ 所有检查通过！
```

### 2. 准备输入文件

需要准备三种输入文件：

#### （1）树文件

支持格式：`.nwk`、`.newick`、`.tree`、`.tre`、`.nex`、`.nexus`、`.nhx`（带注释的 Newick）

> **注意**：输入树必须为**有根树（rooted）**。校验阶段（`phylodater check` 或定年运行前）会直接拒绝无根树。若输入为无根树，请先完成定根（例如在合适的外类群上定根），再运行 PhyloDater。

#### （2）序列比对文件

主要支持格式：FASTA（`.fasta`、`.fa`、`.fna`、`.faa`）、FASTQ（`.fastq`、`.fq`）

可检测格式（仅验证，不用于分析）：PHYLIP（`.phy`）、Clustal（`.aln`）、Stockholm（`.sto`）

> **注意**：PhyloDater 可检测 PHYLIP、Stockholm 等对齐格式，但不支持将这些格式用于分析。如需使用，请先转换为 FASTA。

#### （3）校准点配置文件（YAML 格式）

校准点配置示例（calibrations.yaml）

```yaml
calibrations:
  - name: "Primates"
    node:
      type: "mrca_pair"
      mrca_pair: ["human", "cat"]
    constraint:
      type: "uniform"
      min: 70.0
      max: 90.0

  - name: "Rodents"
    node:
      type: "mrca_pair"
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

> **节点定义类型**：`auto`（根据名称或分类学推断）、`mrca`（多个分类单元列表）、`mrca_pair`（两个分类单元定义节点）。

#### 支持的约束类型

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

#### 方法专属校准注意事项

- **MCMCTree 根节点**：MCMCTree 的 `RootAge` 只能设置为上界，不能真正固定根节点年龄。若根节点使用 `fixed` 约束，PhyloDater 会在执行前报错并提示修改。建议将 MCMCTree 的根节点改为 `uniform` 或 `maximum`，或改用 treePL、PATHd8、LSD2 等方法固定根节点。
- **不同方法可使用不同校准文件**：各定年工具对约束类型的支持并不相同，用户可为不同方法准备不同的 YAML 校准文件，通过 `-c/--calibrations` 分别传入。

### 3. 运行定年分析

```bash
phylodater dating \
    -t tree.nwk \
    -s alignment.fasta \
    -c calibrations.yaml \
    -o results/ \
    --method mcmctree
```

### 4. 使用方法专有参数

```bash
# MCMCTree: 指定 PAML 4.8 版本（建议同时指定 paml_path，避免 mcmctree 及其调用的 baseml/codeml 版本混合）
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method mcmctree \
    --method-args "paml_path=/opt/paml-4.8,paml_version=4.8,num_runs=4"

# treePL: 指定平滑参数
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method treepl \
    --method-args "initial_smooth=100"
```

### 5. 定根选项

```bash
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method pathd8 \
    --reroot midpoint
```

可用策略：`outgroup`、`midpoint`、`mad`、`none`

### 6. 自动添加校准点

```bash
phylodater dating \
    -t tree.nwk -s alignment.fasta -c cal.yaml -o results/ \
    --method mcmctree \
    --auto-calibrate "Cyanobacteriota:2500" "Archaea:3800"
```

## 命令行架构

```
phylodater [通用参数] <子命令> [子命令参数]

子命令:
  dating            执行定年分析
  check (self-test) 系统自检与依赖验证
  fold              按分类群折叠树
```

### 通用参数（所有子命令共享）

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--log-level` | 日志级别（DEBUG/INFO/WARNING/ERROR/CRITICAL） | `INFO` |
| `--log-file` | 日志文件路径（UTF-8，实时 flush） | 无 |
| `--taxonomy-file` | 分类学映射文件 | 无 |
| `--taxonomy-levels` | 自定义分类级别前缀（如 `k:kingdom ss:subspecies`） | 无 |
| `--taxonomy-delimiter-mode` | 格式 A 解析策略（`reverse`/`greedy`/`segment`） | `reverse` |
| `--taxonomy-source-priority` | 分类学来源优先级（`embedded`/`table`） | `table` |
| `--table-sep` | 强制指定分类学表格分隔符 | 自动检测 |
| `--ignore-malformed` | 跳过畸形输入（如不支持的格式，发出警告） | 否 |
| `--version` | 显示版本信息（含 Git 哈希、依赖版本） | - |

### dating 子命令参数

| 参数组 | 参数 | 说明 | 默认值 |
|--------|------|------|--------|
| 输入文件 | `-t/--tree` | 树文件路径 | 必填 |
| | `-s/--sequence` | 序列文件路径 | 必填 |
| | `-c/--calibrations` | 校准点配置文件 | 必填 |
| | `-o/--output` | 输出目录 | 必填 |
| | `--config` | 配置文件路径（YAML） | 无 |
| | `--generate-config` | 生成默认配置文件到指定路径并退出 | 无 |
| 定年方法 | `--method` | 定年方法名称 | 必填 |
| | `--method-args` | 方法专有参数（键值对） | 无 |
| 树处理 | `--reroot` | 定根策略（`outgroup`/`midpoint`/`mad`/`none`） | `none` |
| | `--outgroup` | 外群名称 | 无 |
| | `--multi-tree-mode` | 多棵树模式（`split`/`first`/`last`/`random`/`ask`/`error`） | `error` |
| 校准 | `--auto-calibrate` | 自动添加校准点 | 无 |
| 执行控制 | `-j/--threads` | 并行线程数 | `4` |
| | `--seed` | 全局随机种子 | 自动生成 |
| | `--dry-run` | 干运行模式（仅验证） | 否 |
| | `--mol-type` | 序列分子类型（`DNA`/`RNA`/`protein`） | 自动检测 |
| | `--no-cross-check` | 关闭树与序列的交叉验证 | 否 |
| | `--skip-length-check` | 跳过负分支长度检查 | 否 |
| | `--low-memory` | 低内存模式（单线程，减少 MCMC 参数） | 否 |
| 检查点 | `--resume` | 从检查点恢复（跳过已完成步骤） | 否 |
| | `--clear-checkpoint` | 清除检查点文件后运行 | 否 |
| 输出控制 | `--force` | 覆盖已存在的输出文件 | 否 |
| | `--no-clobber` | 跳过已存在的输出文件 | 否 |

#### fold 子命令参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-t/--tree` | 输入树文件 | 必填 |
| `-o/--output` | 输出目录 | 必填 |
| `--fold-rank` | 折叠到的分类级别（`domain`/`phylum`/`class`/`order`/`family`/`genus`） | `family` |

> **注意**：`fold` 子命令目前为**占位实现（stub）**。它可接收参数并记录目标级别，但不会真正执行树折叠，仅输出“功能尚未实现”警告。请勿用于实际分析。

> **自定义软件路径**：可通过 YAML 配置文件的 `software_paths` 节点为各工具指定二进制路径（如 `mcmctree_bin`、`iqtree_bin`、`r8s_bin`、`treepl_bin`、`pathd8_bin`），或依赖环境 `PATH`。
>
> **多版本 PAML 共存**：PAML 4.8 与 4.10.8 的二进制默认同名（`mcmctree`、`codeml`、`baseml`）。为避免版本混淆，建议：
> 1. 将不同版本安装到独立目录（如 `/opt/paml-4.8`、`/opt/paml-4.10.8`），并在 `software_paths` 中指定 `paml_path` 指向目标版本的根目录；适配器会自动使用 `<paml_path>/bin/mcmctree`。
> 2. 或直接通过 `mcmctree_bin` 给出完整路径，例如 `/opt/paml-4.10.8/bin/mcmctree`。
> 3. 运行前 PhyloDater 会解析可执行文件 banner 中的实际 PAML 版本，并与 `paml_version` 比对**发布系列**：任意 4.10.x 都能满足默认值 `4.10.8`，而跨系列不一致（例如 4.8 二进制）会被明确报错跳过，从而避免控制文件格式错版。

## 分类学解析

PhyloDater 支持两种分类学名称格式：

### 格式 A（嵌入式）

```
GB_GCA_000252485.1_d_Bacteria_p_Cyanobacteriota_c_Cyanobacteriia_o_Cyanobacteriales_f_Prochloraceae_g_Prochloron
```

分隔符：`_d_`（域）、`_p_`（门）、`_c_`（纲）、`_o_`（目）、`_f_`（科）、`_g_`（属）

解析策略（通过 `--taxonomy-delimiter-mode` 控制）：
- `reverse`（默认）：从右向左逆序匹配，避免分类名本身含下划线导致的歧义
- `greedy`：贪婪匹配
- `segment`：从第一个 `_d_` 开始分段截取

### 格式 B（表格分号式）

```
d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales;f__Synechococcaceae;g__Synechococcus;s__
```

**容错规则**：
- 分类值不得包含 `;` 或 `__`，否则报 ERROR 并跳过该行
- `s__` 后无内容解析为 `None`，记录 DEBUG 日志

### 外部分类学表格

当末端分支缺乏完整分类名称时，可通过表格文件补充：

```
# 两列格式（制表符或逗号分隔）
GB_GCA_000252485.1	d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;...
```

使用 `--taxonomy-file` 参数指定表格文件。

**标签一致性检查**：表格中出现的末端若在树中不存在，发出 WARNING 并列出多余项；树中末端在表格中没有记录且 ID 中无法解析时，该末端分类信息置为“未知”。

**来源优先级**：当同一节点同时具有嵌入式和表格两种来源信息时，通过 `--taxonomy-source-priority` 控制优先级（默认 `table` 表格优先）。如有冲突，发出 WARNING 并记录冲突详情。

## 单系群判定

PhyloDater 支持通过指定类群名称自动判定单系群：

### 一般类群判定

给定一个分类标签（如 "Cyanobacteriota"）：
1. 收集所有属于该类群的末端节点
2. 计算这些末端的最近共同祖先（MRCA）
3. 获取 MRCA 的所有后代末端
4. 若后代集合等于步骤 1 的集合 → 单系群；否则 → 非单系群

非单系群时，软件会发出 ERROR 并拒绝操作，同时列出不属于该类群的额外末端。

### 特殊标识符

- **LUCA**：所有 Bacteria 和 Archaea 的 MRCA
- **LBCA**：所有 Bacteria 的 MRCA
- **LACA**：所有 Archaea 的 MRCA

特殊标识符大小写敏感，仅接受大写。计算前首先确认类群存在性，不存在则 ERROR。

## PAML 4.8 支持

PhyloDater 对 PAML 4.8 与 4.10.8 采用统一的三阶段分析流程：

| 阶段 | 功能 | 说明 |
|------|------|------|
| 阶段 0 | 速率估算（baseml 或 codeml，clock=1） | 默认执行。以全局钟拟合数据，从输出中提取 "Substitution rate is per time unit" 的总体替换速率，并以该速率作为 `rgene_gamma` 先验均值。结果保存在输出目录的 `mcmctree_rate_estimation.txt`。 |
| 阶段 1 | Hessian 矩阵计算（usedata=3） | 生成分支长度、梯度、Hessian |
| 阶段 2 | MCMC 采样（usedata=2） | 估算分歧时间 |

### 跳过速率估算（可选）

如需沿用旧的“默认弱信息先验”行为（`rgene_gamma` 均值 = rate_alpha/20，即默认 G(2, 20)），
可设置 `skip_rate_estimation: true`：

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --method-args "skip_rate_estimation=true"
```

也可在 YAML 配置文件的 `software.mcmctree` 块中设置。默认 `false`（执行速率估算）。

### 各 PAML 版本在速率估算上的行为差异

| 版本 | "Substitution rate is per time unit" 输出 | 说明 |
|------|-------------------------------------------|------|
| ≤ 4.9h（含 4.8） | 输出 | 正常 |
| 4.9i | **缺失** | 官方确认的 bug（见 PAML changelog v4.9j 条目）。PhyloDater 不支持 4.9i，此版本下速率估算会回退为默认先验。 |
| 4.9j / 4.10.x（含 4.10.8） | 输出 | 4.9j 已修复该 bug，行为与 4.8 一致。仅格式上有细微差别：4.10.8 在速率标题与数值之间多打印一个空行，PhyloDater 的提取逻辑已兼容。 |

实测验证：同一数据下 PAML 4.8 与 4.10.8 的 baseml（clock=1）速率估算结果一致
（差异仅出现在小数第 5 位，例如 4.695519 与 4.695467）。

**注意**：该速率行只有在树中存在绝对年龄标注（`@`，即化石或根年龄约束）时才会输出。
PhyloDater 会自动把根节点校准的代表性年龄（单位 Ma，内部转换为 Ga，与树内校准一致）
以 `@age` 形式写入速率估算树；若既无根校准也未配置 `root_age`，速率估算无法进行，
将回退为默认先验并记录 WARNING。

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --method-args "paml_version=4.8"
```

## MCMC 收敛性检验

对于 MCMCTree，当 `num_runs >= 2` 时自动执行收敛性诊断：

- **PSRF**（Potential Scale Reduction Factor）：Gelman-Rubin 统计量，理想值 < 1.2
- **ESS**（Effective Sample Size）：有效样本量，建议 ≥ 100

收敛性结果会记录在日志中。若 PSRF ≥ 2.0，程序报错；若 1.2 ≤ PSRF < 2.0，发出警告。

## 输入验证

### 树文件验证

- 括号平衡检查
- 分支长度非负（负值 → CRITICAL 终止）
- 重复节点名检测（末端重复 → ERROR，内部重复 → WARNING）
- 恶意字符注入检测（控制字符、Unicode 双向文本）

### 序列文件验证

- 格式自动识别（FASTA、FASTQ）
- 序列字母表验证（DNA、RNA、protein）
- 序列 ID 唯一性检查
- 序列长度一致性检查

### 交叉验证

当同时提供树文件和序列文件时，默认执行末端标签双向匹配，确保两个集合完全一致。可通过 `--no-cross-check` 关闭。

### 对抗性输入防护

- **恶意符号注入**：检测控制字符和 Unicode 双向文本覆盖字符
- **分类学循环依赖**：检测表格中的循环引用（如 `d__A;p__B` 和 `d__B;p__A`）
- **空文件**：树文件或序列文件大小为 0 字节时终止

## 多棵树处理

默认认为输入的树文件仅含一棵树。若检测到多棵树：

| 模式 | 行为 |
|------|------|
| `split` | 分别处理每棵树 |
| `first` | 仅处理第一棵树 |
| `last` | 仅处理最后一棵树 |
| `random` | 随机选择一棵树处理 |
| `ask` | 交互选择提示，然后以退出码 2 退出 |
| `error`（默认） | 报错并终止（非交互式默认行为） |

## 日志系统

### 日志格式

```
2025-03-21T10:15:30.123 | INFO     | module/func 开始解析树文件
2025-03-21T10:15:30.456 | WARNING  | module/func 检测到大型树 (15000 个末端节点)
```

### 日志文件

使用 `--log-file` 参数将日志同时写入文件：

```bash
phylodater dating -t tree.nwk -s aln.fasta -c cal.yaml -o out/ \
    --method mcmctree \
    --log-file phylodater.log
```

## 输出结果

运行完成后，在输出目录中会生成：

```
results/
├── summary.txt             # 执行摘要
├── comparison_table.txt    # 文本格式比较表
├── comparison_table.tsv    # TSV 格式比较表
├── comparison_plot.png     # 比较图（单方法运行也会生成；多方法比较时信息量更大）
├── phylodater.log          # 运行日志
├── runtime_metadata.json   # 运行时元数据快照
├── runtime_config.yaml     # 运行时配置快照
├── mcmctree/               # MCMCTree 结果
│   ├── run1/
│   │   ├── FigTree.tre     # 定年树（含校准标注）
│   │   └── mcmc.txt        # MCMC 样本
│   └── run2/
├── lsd2/                   # LSD2 结果
├── treepl/                 # treePL 结果
└── pathd8/                 # PATHd8 结果
```

#### 校准标注格式

定年树文件中，校准节点以 `[&height=AGE,height_95%_HPD={LOWER,UPPER}]` 格式标注。

## 配置文件

使用 `--config` 参数指定 YAML 配置文件，可替代命令行参数：

```yaml
software_paths:
  # 多版本 PAML 共存时，建议只指定 paml_path，适配器会自动使用
  # <paml_path>/bin/mcmctree、<paml_path>/bin/codeml、<paml_path>/bin/baseml
  paml_path: /opt/paml-4.10.8
  # 如需精确指定某个版本的二进制（例如 PATH 中有同名命令），可写完整路径：
  # mcmctree_bin: /opt/paml-4.10.8/bin/mcmctree
  iqtree_bin: /path/to/iqtree2
  r8s_bin: /path/to/r8s
  treepl_bin: /path/to/treepl
  pathd8_bin: /path/to/PATHd8

software:
  common:
    nthreads: 8
    seed: 42

  mcmctree:
    clock: 2
    num_runs: 4
    burnin: 20000
    nsample: 50000

  r8s:
    method: pl
    smoothing: 1000

  treepl:
    initial_smooth: 100

  lsd2:
    model: LG+G
```

配置优先级：默认值 < YAML 文件 < 命令行参数。PhyloDater 会按各方法的模式检查并校验所有值。

## 错误码

| 退出码 | 含义 | 说明 |
|--------|------|------|
| `0` | 成功 | 分析正常完成 |
| `1` | 运行错误 | 依赖缺失、自检失败等 |
| `2` | 参数错误 | 命令行参数无效 |
| `3` | 数据错误 | 输入文件格式或内容错误 |
| `130` | 用户中断 | Ctrl+C 终止 |

## 大规模数据处理

对于超过 10,000 个末端的大型树：

- 软件会在 INFO 日志中提示可能需要的计算资源
- 使用 `--low-memory` 减少资源占用（单线程，减少 MCMC 参数）
- 使用 `--threads` 增加并行线程数
- 建议至少 16 GB 内存

## 用户中断与资源清理

按 `Ctrl+C` 可优雅终止程序：

- 输出当前进度（已处理的树数量等）
- 关闭所有打开的文件句柄
- 删除未完成的临时输出
- 以退出码 `130` 退出

## 跨平台与编码

- **编码**：所有文件读写使用 UTF-8，遇到非法字节时尝试 utf-8-sig
- **路径处理**：使用 `pathlib.Path`，兼容 Windows、macOS 和 Linux
- **换行符**：读取时自动归一化，写入统一使用 `\n`（Unix 风格）
- **临时文件**：使用 `tempfile` 模块，确保异常退出时自动清理

## 版本管理

- 版本号遵循语义化版本 `MAJOR.MINOR.PATCH`；发布版本静态锁定为 `0.1.0`（写于 `pyproject.toml` 与 `phylodater/_version.py`）——包内不再烘焙由 Git 派生的 `.dev`/哈希等版本迭代串
- 在 Git 工作区内，`--version` 会额外报告当前短提交哈希，仅用于支持/排障
- CI 自动运行 `phylodater check` 验证所有测试用例

`--version` 输出示例：
```
PhyloDater 0.1.0
License: MIT License
Updated: 2026-07-12
git: abc1234
deps: BioPython=1.87, ETE3=3.1.3, DendroPy=5.0.8, NumPy=2.4.6
```

#### 测试通过环境

发布验证在 Python 3.12.14 的 conda 环境中完成，该环境装有全部七个引擎；操作系统、
库版本、引擎版本及其来源见 [docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md) §1。
套件具备可移植性：在未装齐引擎的机器上运行 `python -m pytest tests` 仍然通过，因为
依赖外部引擎的测试会带明确原因跳过而非失败。

> **覆盖率说明**：在验证环境上实测 `pytest tests --cov=phylodater` 为 **67.87%**
> （加 `-m "not slow"` 为 67.89%，可见该数字不依赖引擎是否装齐）。`pyproject.toml`
> 里的 CI 门槛因此从占位值 40 上调到 `fail_under = 60`——相对实测值留出约 8 个百分点
> 的余量以吸收 CI 平台差异；它是下限，不是目标值。

## 示例数据

`examples/` 目录包含示例文件供测试使用：

```
examples/
├── tree.nwk            # 示例树文件（Newick 格式）
├── alignment.fasta     # 示例序列比对文件（FASTA 格式）
└── calibrations.yaml   # 示例校准点配置文件
```

> **数据来源**：示例数据为模拟数据，仅供功能测试使用，不涉及任何真实生物样本或隐私数据。

运行示例：
```bash
phylodater dating \
    -t examples/tree.nwk \
    -s examples/alignment.fasta \
    -c examples/calibrations.yaml \
    -o results/ \
    --method mcmctree
```

## 项目结构

```
phylodater/
├── phylodater/
│   ├── adapters/            # 定年方法适配器（每个引擎一个，DatingMethod[ConfigT]）
│   ├── cli/                 # 命令行接口（子命令架构）
│   ├── core/                # 核心接口（method_interface、注册表）、异常、流水线、
│   │                        #   批处理、结果比较报告
│   ├── infrastructure/      # 日志、配置、进程执行、检查点、信号处理、
│   │                        #   插件、自检
│   ├── models/              # 数据模型（树、约束、校准点、结果）
│   ├── services/            # 分类学解析、校准加载/解析、验证、定根、名称映射
│   ├── viz/                 # 可视化（定年图、树/地理图、后端）
│   └── viz/data/            # 随包发布的地质年代表（CSV package data）
├── examples/                # mini 示例数据集（含 truth.yaml 与校验和）
├── test-data/
│   └── benchmark/           # 带真值的参考基准数据集 + 生成脚本
│       ├── tree.nwk / tree_ml.nwk / alignment.fasta
│       ├── truth.yaml / calibrations.yaml / checksums.sha256
│       └── results/         # 归档的 scorecard.tsv、node_ages.tsv、
│                            #   benchmark_results.json（含 ml_tree/ 变体）、
│                            #   PYTHON_SUPPORT_MATRIX.txt
├── tests/
│   ├── unit/                # 快速、全 mock
│   ├── functional/          # CLI / 流水线行为
│   ├── integration/         # 端到端流程场景
│   └── validation/          # 真实引擎端到端验证（标记：slow、
│                            #   requires_mcmctree、requires_lsd2 …）
├── scripts/
│   ├── run_reference_benchmark.py   # 唯一官方基准：按真值打分
│   ├── check_docs.py                # 双语政策与链接门禁
│   ├── generate_checksums.py        # 数据集校验和清单
│   ├── probe_pyr8s_constraints.py   # pyr8s 能力探测（结论进 docs）
│   ├── tier1_manual_comparison.py   # 手动 vs 包装层 对照实验
│   └── degradation_experiment.py    # 约束降级实验
├── tutorials/               # Jupyter 快速上手（英文 + _CN 两份 notebook）
└── docs/                    # 文档（英文权威 + _CN 中文）
    ├── ARCHITECTURE.md / ARCHITECTURE_DECISIONS.md / ERROR_CODES.md
    ├── EXTERNAL_TOOLS.md / TESTING.md / TEST_RESULTS.md / CONDA_PACKAGING.md
```

## 测试与验证

| 内容 | 位置 |
|------|------|
| 如何测试以及为何这样测（方法、门槛、复现命令） | [docs/TESTING_CN.md](docs/TESTING_CN.md) |
| 实测结果：七个引擎、对仿真真值的精度、确定性、发现并修复的缺陷 | [docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md) |
| 带真值的数据集 + 生成脚本 + 归档结果 | [`test-data/benchmark/`](test-data/benchmark/) |
| 复现整套基准 | `python scripts/run_reference_benchmark.py` |
| 每个外部引擎的获取方式（版本、校验和、源码构建） | [docs/EXTERNAL_TOOLS_CN.md](docs/EXTERNAL_TOOLS_CN.md) |
| 发布验证套件（真实引擎、端到端） | [`tests/validation/`](tests/validation/) |

```bash
phylodater check                     # 26 项自检，无需外部引擎
python -m pytest tests -q            # 全部自动化测试（1489 项）
python scripts/run_reference_benchmark.py   # 七个引擎对真值评分
```

## 依赖库

| 库 | 版本 | 许可证 |
|----|------|--------|
| BioPython | >=1.79,<2.0 | BSD-like |
| ETE3 | >=3.1.2,<4.0（可选） | GPLv3 |
| DendroPy | >=4.5.0 | BSD 3-Clause |
| NumPy | >=1.21.0 | BSD 3-Clause |
| Pandas | >=1.3.0 | BSD 3-Clause |
| PyYAML | >=6.0,<7.0 | MIT |
| SciPy | >=1.9.0,<2.0 | BSD 3-Clause |
| Matplotlib | >=3.5.0 | PSF-based |
| psutil | >=5.9.0,<6.0 | BSD 3-Clause |

> **注意**：ETE3 为可选依赖（声明于 `extras_require["ete3"]`），默认通过 pip 安装。其余库均为核心依赖。

## 文档

- [快速入门指南](QUICKSTART_CN.md)
- [文档目录](docs/README_CN.md)
- [测试指南](docs/TESTING_CN.md)与[测试结果](docs/TEST_RESULTS_CN.md)
- [外部定年引擎](docs/EXTERNAL_TOOLS_CN.md)
- [API 参考](docs/api/)
- [架构](docs/ARCHITECTURE_CN.md) · [设计决策](docs/ARCHITECTURE_DECISIONS_CN.md) · [错误码](docs/ERROR_CODES_CN.md)
- [贡献指南](CONTRIBUTING_CN.md)

所有文档均有中英两个版本，其中英文为权威版本（并列时排在前面）。

## 贡献

欢迎贡献代码，请参阅 [CONTRIBUTING_CN.md](CONTRIBUTING_CN.md) 了解如何参与开发。

- **问题报告**：[GitHub Issues](https://github.com/ZengZichao/phylodater/issues)
- **功能建议**：[GitHub Discussions](https://github.com/ZengZichao/phylodater/discussions)
- **联系邮箱**：曾子超（Zichao Zeng）— zengzichao@sjtu.edu.cn · ORCID: [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

## 许可证

本项目采用 MIT 许可证，详见 [LICENSE](LICENSE) 文件。

## 引用

如果你在研究中使用了 PhyloDater，请引用：

> 曾子超（Zeng, Zichao）. (2026). PhyloDater: A multi-software parallel platform for phylogenetic molecular dating. GitHub repository. https://github.com/ZengZichao/phylodater （ORCID: 0000-0001-6553-970X）· DOI: [10.5281/zenodo.23067037](https://doi.org/10.5281/zenodo.23067037)

BibTeX 格式：
```bibtex
@software{zeng2026phylodater,
  author = {Zeng, Zichao},
  title = {PhyloDater: A multi-software parallel platform for phylogenetic molecular dating},
  year = {2026},
  url = {https://github.com/ZengZichao/phylodater},
  doi = {10.5281/zenodo.23067037},
  version = {0.1.0},
  note = {ORCID: 0000-0001-6553-970X}
}
```

完整的引用指引（机器可读 CFF、引擎引用、本软件的 Zenodo DOI
10.5281/zenodo.23067037，以及"目前没有论文 DOI"这一事实）
见 [CITATION_CN.md](CITATION_CN.md)。
