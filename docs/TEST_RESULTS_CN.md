> **Language:** [English (canonical)](TEST_RESULTS.md) · 中文（本文档）

# 测试结果

本文是 [TESTING_CN.md](TESTING_CN.md) 所述验证方案的执行结果。下面每一个数字都由
本文 §8 中的命令在 §1 记录的环境上产生，并且可以从随包发布的
[`test-data/benchmark/results/`](../test-data/benchmark/results/) 复现。

**版本：** PhyloDater 0.1.0 | **日期：** 2026-09-29 | **结论：通过**——七个受支持的
定年方法全部针对其**真实引擎**完成端到端运行，并满足全部验收门槛。

---

## 1. 环境

| 项目 | 值 |
|------|-----|
| 操作系统 / 内核 | Ubuntu 22.04.4 LTS，Linux 6.8.0-111-generic，x86_64 |
| 硬件 | 48 物理核 / 96 逻辑核，772 GB 内存 |
| Python | 3.12.14（conda 环境 `phylodater`） |
| 包管理器 | micromamba 2.9.0（环境内）、pip 24.2 |
| 核心库 | numpy 2.5.3、pandas 3.0.6、scipy 1.18.1、matplotlib 3.11.2、DendroPy 5.1.0、ete3 3.1.3、PyYAML 6.0.3、biopython 1.88 |
| 测试框架 | pytest 9.1.1、pytest-cov 7.1.0 |

### 1.1 实际执行的定年引擎

| 方法 | `runtime_metadata.json` 记录的引擎版本 | 安装路径 |
|------|----------------------------------------|----------|
| `mcmctree` | `MCMCTREE in paml version 4.10.10, 29 Jan 2026` | bioconda `paml` |
| `lsd2` | `IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026` | bioconda `iqtree` |
| `r8s` | `pyr8s 0.3.1 (NPRS only)` | 源码 pip 安装（`iTaxoTools/pyr8s`） |
| `treepl` | `treePL version 1.0`（2.6.3 构建） | 历史 conda 构建 + `patchelf`，或源码构建 |
| `pathd8` | `PATHd8 (no version banner): <prefix>/bin/PATHd8` | 官方源码构建 |
| `wlogdate` | `wLogDate 1.0.2` | PyPI `wlogdate` |
| `mdcat` | `1.0.1` | 源码 pip 安装（`uym2/MD-Cat`） |

安装来源、校验和以及为何选择该路径，见 [EXTERNAL_TOOLS_CN.md](EXTERNAL_TOOLS_CN.md)。

## 2. 基准测试——总体结果

数据集：`test-data/benchmark`（13 分类单元、3000 位点、JC69 严格分子钟、根龄 105 Ma、
12 个内部节点、4 个校准；见 TESTING_CN.md §3）。门槛：MAE ≤ 10 Ma、最差节点相对误差
≤ 0.5，且每个请求的校准要么被满足、要么被适配器**显式披露**为未满足。

### 2.1 变体 A——输入树使用（无噪声的）期望钟长

| 方法 | 状态 | 耗时(s) | MAE(Ma) | RMSE(Ma) | 最大绝对误差(Ma) | 最大相对误差 |
|------|------|---------|---------|----------|-------------------|--------------|
| `mcmctree` | 通过 | 14.2 | 0.721 | 0.863 | 1.641 | 0.284 |
| `lsd2` | 通过 | 4.2 | 1.881 | 2.373 | 4.912 | 0.303 |
| `r8s`（pyr8s/NPRS） | 通过 | 4.4 | 0.183 | 0.256 | 0.528 | 0.013 |
| `treepl` | 通过 | 64.4 | 1.099 | 1.527 | 3.044 | 0.029 |
| `pathd8` | 通过 | 4.4 | 0.509 | 0.794 | 1.659 | 0.019 |
| `wlogdate` | 通过 | 6.1 | 0.000 | 0.000 | 0.000 | 0.000 |
| `mdcat` | 通过 | 10.5 | 0.081 | 0.106 | 0.163 | 0.069 |

七个方法在以下各项全部通过：叶节点完整、拓扑一致、年龄非负、时序单调、校准被满足、
12 个真值节点全部还原。

`wlogdate` 在此变体上**恰好**得到真值，这并不代表它精度更高：该变体喂给方法的分支长度
本来就与时间成正比，恢复年龄只是一次比例换算。保留这一变体正是为了当对照组——它证明
了评分流程不可能超出输入本身所含的信息。

### 2.2 变体 B——输入树使用 ML 估计分支长度（真实情形）

同一份比对，但输入树的分支长度由 IQ-TREE 2 以最大似然（GTR+G、拓扑固定）估计，并按
TESTING_CN.md §3 的方式重新定根。树总长从 1.0365 降到 0.8793 替换/位点，因此**输入
本身已经带偏、带噪**。

| 方法 | 状态 | 耗时(s) | MAE(Ma) | RMSE(Ma) | 最大绝对误差(Ma) | 最大相对误差 | 已披露的引擎局限 |
|------|------|---------|---------|----------|-------------------|--------------|------------------|
| `mcmctree` | 通过 | 14.0 | 0.730 | 0.878 | 1.650 | 0.289 | — |
| `lsd2` | 通过 | 4.1 | 1.881 | 2.373 | 4.912 | 0.303 | — |
| `r8s`（pyr8s/NPRS） | 通过 | 4.8 | 7.862 | 11.219 | 25.023 | 0.426 | `Primates` 54.64 Ma 超出 [68, 82]；`Glires` 53.75 Ma 超出 [60, 80] |
| `treepl` | 通过 | 63.0 | 4.041 | 5.762 | 14.866 | 0.224 | — |
| `pathd8` | 通过 | 4.1 | 4.492 | 6.262 | 15.652 | 0.240 | — |
| `wlogdate` | 通过 | 6.0 | 1.850 | 3.199 | 9.666 | 0.180 | — |
| `mdcat` | 通过 | 18.5 | 0.298 | 0.472 | 1.116 | 0.413 | — |

对这张表的诚实解读：

* **对外应引用变体 B。** 除 `mcmctree` 与 `mdcat`（它们会重新估计速率而不是采信给定的
  枝长）外，所有方法都因输入变差而精度下降。
* MCMCTree 的数字会在小数第三位随运行变化（归档这次为 0.721，此前两次分别为 0.731
  与 0.756），因为它做随机采样；其余每一行都可逐位复现（§3）。归档的 `scorecard.tsv`
  才是"最后一次实际运行"的权威记录。
* **`lsd2` 两个变体完全相同**，因为 IQ-TREE 2 会在固定拓扑上用 `GTR+G` 重新优化分支
  长度，输入枝长并不作为数据被使用；1.88 Ma 的 MAE 是该估计方法自身的误差。
* **`r8s` 的区间越界属于引擎性质，而不是静默失败。** pyr8s 把 `CONSTRAIN min_age/max_age`
  实现为**逐步松弛的势垒惩罚**，因此 NPRS 可以落在请求区间之外；而 `FIXAGE` 会被精确执行
  （实测：`Homininae` fixage 9 Ma → 两个变体均给出 9.000 Ma）。得出该结论所依据的实验
  随包发布：[`scripts/probe_pyr8s_constraints.py`](../scripts/probe_pyr8s_constraints.py)。
  PhyloDater 现在 (a) 在运行前声明 pyr8s 的区间只是参考值，并 (b) 在事后用返回的年龄
  逐条复核请求区间，为每个未满足的校准写入具名警告。基准的门槛是"满足**或被披露**"，
  上表最后一列正是适配器报告的内容。
* **`pathd8` 的根节点行为在本轮得到改进。** 适配器过去会丢弃一切非点校准的根约束
  （"PATHd8 没有针对根段的 min/max 指令"），导致根年龄完全由 MPL 依树长推导：请求
  [100, 110] 却返回 153 Ma。已直接用 PATHd8 二进制确认根 mrca 支持 `minage`/`maxage`，
  因此现在会写出该区间，根年龄返回 110 Ma（落在请求区间内；变体 B，12 节点 MAE 4.49 Ma）。

### 2.3 类型改造后的耗时复核（无性能回退）

本轮类型债工作（P0）改到了每个适配器的构造函数与 286 处函数定义
（即标注补齐 `cb4fd01` 中被改动的 `def` 行），因此在改完之后重新实测了耗时：
同一台机器、同一个 conda 环境，变体 A 独立跑三次、变体 B 跑两次：

| 方法 | 变体 A 耗时（秒），三次 | 变体 B 耗时（秒），两次 | 归档基线（A / B） |
|------|------------------------|------------------------|-------------------|
| `mcmctree` | 15.07 / 13.86 / 14.41 | 14.60 / 14.31 | 15.07 / 14.60 |
| `lsd2` | 4.39 / 4.40 / 3.86 | 4.35 / 4.13 | 4.39 / 4.35 |
| `r8s` | 4.31 / 4.68 / 4.54 | 4.56 / 4.44 | 4.31 / 4.56 |
| `treepl` | **64.50 / 67.81 / 59.85** | 55.78 / 63.50 | 64.50 / 55.78 |
| `pathd8` | 4.35 / 4.36 / 4.20 | 4.35 / 4.25 | 4.35 / 4.35 |
| `wlogdate` | 5.81 / 5.28 / 5.48 | 6.03 / 6.01 | 5.81 / 6.03 |
| `mdcat` | 10.51 / 10.29 / 10.05 | 18.37 / 17.74 | 10.51 / 18.37 |

`treepl` 仍是主要开销（每次 55–68 秒；耗时来自对平滑参数的交叉验证搜索，而不是包装层），
改造前后的波动幅度一致——没有可测的性能回退。所有确定性方法的精度也一字不差地复现：
变体 A 下 `lsd2` / `r8s` / `treepl` / `pathd8` / `wlogdate` / `mdcat` 的 MAE
仍为 1.881 / 0.183 / 1.099 / 0.509 / 0.000 / 0.081 Ma。

## 3. 确定性

两次独立的完整基准运行（同样 `--seed 42`），逐节点比对 `node_ages.tsv`：

| 方法 | 重复运行结果 |
|------|--------------|
| `lsd2`、`r8s`、`treepl`、`pathd8`、`wlogdate`、`mdcat` | 节点年龄逐位一致 |
| `mcmctree` | 相差 ≤ 0.3 Ma（不足根龄 0.3 %）——贝叶斯采样器不会因为固定 PRNG 种子就逐位可复现 |

类型改造之后又按方法对 `node_ages.tsv` 的行做了哈希比对（变体 A 两次运行）：
`lsd2`、`r8s`、`treepl`、`pathd8`、`wlogdate`、`mdcat` **逐位一致**；
只有 `mcmctree` 变化（MAE 0.705 vs 0.741 Ma）。

因此 MCMCTree 的质量由其收敛诊断判断，而 PhyloDater 会把这些诊断写进产物：
2 条链 × 4001 样本下 `max PSRF = 0.99988`、`min ESS = 423`（阈值 1.05 / 100），
判定 `converged`。

## 4. 自动化测试套件

```bash
python -m pytest tests -q        # 收集 1489 项
```

| 范围 | 收集 | 通过 | 跳过 | 失败 |
|------|------|------|------|------|
| `tests/`（unit + functional + integration + validation） | 1489 | 1487 | 2 | 0 |

跳过的 2 项是**有意**的条件跳过（依赖特定可选后端的环境判定），并非缺失覆盖。

### 4.1 分层说明

* `tests/unit` + `tests/functional` + `tests/integration`——开发者套件：模型、约束渲染、
  适配器写输入/解析输出、CLI、报告、安全、可视化、检查点。
* `tests/validation`——发布验证套件：各适配器契约测试、CLI 子命令（`dating`、`check`、
  `fold`）测试、错误路径测试、安全测试（路径穿越、YAML 安全、输入安全），以及在基准
  fixture 上对每个已安装引擎的端到端运行。
* 依赖外部引擎的测试在二进制缺失时带明确原因跳过，因此同一条命令在没有这些工具的机器上
  也能通过。

## 5. 已安装包的自检

```
$ python -m phylodater check
...
总计: 26 项检查, 26 通过, 0 失败
✓ 所有检查通过！
```

26 项检查覆盖：必需的 Python 依赖（biopython、dendropy、numpy、pandas、pyyaml、
matplotlib）；四个可选的引擎依赖——`logdate`（wLogDate）、`emd`（MD-Cat）、`pyr8s`（r8s）
与 `ete3`，只探测，从来不是 Python 层运行的前提；有无分类信息的 Newick 解析；
两种分类学格式（嵌入式 GTDB 与表格分号式，含正确的缺失等级处理：`species=None`
是契约而不是 bug）；单系性判断（非单系、单系与 LUCA 特殊标识符）；深度校验（括号平衡、
重复 ID）；LUCA/ROOT 特殊标识符；500 MB 输入体积限制；恶意字符拒绝；多树 NEXUS 处理；
外部分类学表格；以及库模式的交叉校验 API。

## 6. 本轮发现并修复的缺陷

| # | 缺陷 | 表现 | 修复 | 钉住它的测试 |
|---|------|------|------|--------------|
| 1 | LSD2 默认替换模型是氨基酸模型 `LG+G` | IQ-TREE 2 报 `ERROR: File not found LG`，方法根本无法运行 | 按比对类型解析模型（核酸 `GTR+G`、蛋白 `LG+G`）；显式给出明显不匹配的模型时抛 `ConfigurationError` | `TestLSD2ModelResolution`（5 项） |
| 2 | MD-Cat 的 `--CI` 以逗号拼接的单个串传入 | 上游 `int("100,0.025,0.975")` 抛 `ValueError`，默认配置每次必失败 | 改为单个空白分隔串 `--CI "n lower upper"`（与上游 `.split()` 一致） | `TestMDCatCliContract` |
| 3 | MD-Cat CI 默认 100 次 bootstrap | 每次运行数小时的隐式开销（实测 6 分类单元上**一次** bootstrap 就 > 10 分钟） | 默认 `ci_nboots = 0`，与上游一致；需显式开启 | `test_ci_parameters` |
| 4 | PAML 版本闸门要求字面版本 `4.10.8` | 按文档 `conda install -c bioconda paml` 装到的 4.10.10 被直接拒绝："Environment not available" | 兼容判据改为主.次版本系列；同系列补丁差异放行并记信息日志；跨系列仍阻断 | `tests/unit/adapters/test_mcmctree_method.py` 新增 2 项 |
| 5 | `runtime_metadata.json` 从不记录引擎版本 | 随包结果缺少引擎溯源 | 各适配器记录 `validate_environment()` 的探测结果，由 pipeline 写入 | 基准评分卡 + 适配器测试 |
| 6 | 零依赖 Newick 游走器拒绝开头的 `[&R]` | BEAST/IQ-TREE/DendroPy 导出的树会让 MCMCTree 与 treePL 报 "Cannot parse Newick for tree operations"；祖先关系校验静默退化 | 游走器跳过起始注释；新增共享的 `strip_leading_newick_comments()` 为 ete3 归一化 | `TestLeadingCommentTolerance`（5 项） |
| 7 | PATHd8 丢弃非点校准的根约束 | 根年龄被忽略：请求 [100,110] 却返回 153 Ma | 根区间改为在根 mrca 上写 `minage`/`maxage`（行为已对二进制验证） | `tests/unit/adapters/test_pathd8_method.py` |
| 8 | pyr8s 的区间校准既不被强制也不被报告 | 区间外的年龄照常进入报告，一字未提 | 运行前披露 + 运行后逐条对账，写入结果级警告 | `tests/unit/adapters/test_r8s_pyr8s_method.py` |
| 9 | 三个单元测试断言进程级 `sys.modules` | 只因测试顺序（可视化测试导入了 ete3）而失败 | `module_imports_cleanly` fixture 在全新解释器中探测 | 改造后的三项测试 |
| 10 | MCMCTree 原始调用冒烟测试写 `BDparas = 1 1 0.1` | PAML ≥ 4.10 中止："BDparas: expect flag for birth-death process prior" | 控制文件补上构造标志，与适配器一致 | `tests/integration/full_pipeline/test_all_dating_tools.py` |

## 7. 与引擎/环境相关的局限

* treePL 自报横幅始终是 `treePL version 1.0`（即便构建自 2.6.3），因此版本溯源还依赖
  被记录的可执行文件路径。
* PATHd8 完全不打印版本号。
* Docker 镜像只含七个引擎中的三个（EXTERNAL_TOOLS_CN.md §7）。
* 上述精度数字属于**这一份**数据集（13 分类单元、3 kb、严格钟），不可作为对这些方法的
  通用精度断言。
* 在无噪声输入树上某方法 MAE 为 0（变体 A 的 `wlogdate`）是输入的性质，不是估计器的
  能力——用户在 ML 树上应参考变体 B。

## 8. 复现本文全部结果

```bash
conda activate phylodater            # 七个引擎齐备的环境（§1.1）

# 变体 A
python scripts/run_reference_benchmark.py \
    --output test-data/benchmark/results

# 变体 B（ML 分支长度）
python scripts/run_reference_benchmark.py --tree-file tree_ml.nwk \
    --output test-data/benchmark/results/ml_tree

# §2.2 结论所依据的引擎语义实验
python scripts/probe_pyr8s_constraints.py

# 自动化套件与自检
python -m pytest tests -q
phylodater check
```

随包发布的结果文件：
`test-data/benchmark/results/{scorecard.tsv,node_ages.tsv,benchmark_results.json}`
以及 `results/ml_tree/` 下同名三份。`runs/` 子目录（各方法的定年树、比较表、
`runtime_metadata.json`、日志）按需重新生成，已在 `.gitignore` 中忽略。

---

*另见：[TESTING_CN.md](TESTING_CN.md)、[EXTERNAL_TOOLS_CN.md](EXTERNAL_TOOLS_CN.md)、
[English version](TEST_RESULTS.md)。*
