> **Language:** [English (canonical)](CHANGELOG.md) · 中文（本文档）

# 变更日志

**PhyloDater** 的所有值得关注的变更都记录在本文件中。

本文件采用 [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) 格式，
并且项目遵循[语义化版本](https://semver.org/spec/v2.0.0.html)——这一点在
[`README.md`](README_CN.md)（"版本管理"）与 [`pyproject.toml`](pyproject.toml)
（`Development Status :: 4 - Beta`）中均有声明。

> **本文件是如何写成的。** 它于 2026-10-08 依据仓库自身已经记载的内容重建——
> [`docs/TEST_RESULTS.md`](docs/TEST_RESULTS_CN.md)、
> [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE_CN.md)、
> [`docs/EXTERNAL_TOOLS.md`](docs/EXTERNAL_TOOLS_CN.md)、
> [`docs/TESTING.md`](docs/TESTING_CN.md) 以及
> [`README.md`](README_CN.md) 中的发布说明。本仓库公开的历史只有一个被压缩过的
> 提交，因此**这里不声称任何中间发布日期或版本递增**：下述内容只描述 0.1.0
> 这一版，以及它自己的验证活动记录在案的缺陷修复。

## [0.1.0] - 2026-09-30

首个公开发布版本。日期取自 `CITATION.cff`（`date-released: 2026-09-30`、
`version: 0.1.0`）；下述数字所依据的发布验证流程执行于 2026-09-29
（`docs/TEST_RESULTS.md:10`，结论 **PASS**），其归档命令输出捕获于 2026-09-30
（`test-data/benchmark/results/validation_log.txt:3`）。

### 新增

- **一个接口之下的七个引擎适配器**（`mcmctree`、`lsd2`、`treepl`、`pathd8`、
  `r8s`、`wlogdate`、`mdcat`），都通过 `DatingMethodRegistry` 注册，并共享同一
  抽象契约：`validate_environment` / `prepare_inputs` / `execute` /
  `parse_results`（`docs/ARCHITECTURE.md` 的 "Adapter Layer"）。`r8s` 与其 Python
  移植版 `pyr8s` 由注册名为 `r8s` 的同一个适配器服务。
- **八类校准约束**加其抽象基类（`models/constraints.py`）：固定年龄、均匀分布、
  软下界、最大年龄、软双界、gamma 先验、偏正态与偏 t——并具备逐引擎渲染
  （`to_software_format`），以及当目标引擎无法精确表达某个先验时的显式降级报告
  （`docs/ARCHITECTURE.md` 的 "Constraint Type Hierarchy"；`docs/TESTING.md` §8）。
- **校准分类群的自动节点解析**：特殊标识符（LUCA/LBCA/LACA/ROOT）之后是五级递进
  匹配，从尖端名精确匹配到基于分类学、再到多样性最大化的候选选择
  （`services/calibration_resolver.py`）。
- **并行多方法流水线**：同一输入在线程池中并行执行多个方法
  （`core/pipeline.py`），上限 8 个工作线程，可用 `--threads` /
  `PHYLODATER_THREADS` 配置；各方法种子由单个 `--seed` 确定性派生，使随机引擎
  可复现（`docs/ARCHITECTURE.md` 的 "Data Flow"；`docs/TESTING.md` §7）。
- **检查点与续跑**，覆盖流水线与 MCMCTree 两阶段工作流（`--resume`、
  `--clear-checkpoint`），并有优雅的 `SIGINT` 处理与资源清理（退出码 130）。
- **三个 CLI 子命令**——`dating`、`check`（26 项自检）与 `fold`——配套成文的
  退出码 0/1/2/3/130（`README.md` 的 "CLI Architecture"、"Exit Codes"）。
- **标准化的报告产物**：`comparison_table.txt`、`comparison_table.tsv`、
  `summary.txt`、`comparison_plot.png`、各方法的定年树
  （`{method}_dated_tree.nexus` / `.nhx`）、`phylodater.log`，以及
  `runtime_config.yaml` / `runtime_metadata.json`。运行时元数据现在会记录每个
  适配器探测到的引擎版本，于是定年树携带产出它的二进制文件的溯源信息
  （`docs/EXTERNAL_TOOLS.md` §6）。
- **可视化层**（`phylodater.viz`）：`DatingFigure` 绘制定年树并按地质年代做
  颜色/大小映射、支持支系高亮与坍缩，后端注册表可插拔——默认 `matplotlib`，
  另有交互式 `plotly` 后端（通过 `interactive` extra 选用，见 `pyproject.toml`）——
  并随包附带 ICS 地质年代表数据。
- **输入校验与安全加固**：括号平衡与重名检测、序列字母表与长度一致性检查、编码
  披露与 BOM/GBK 修复提示、零宽字符与 Unicode 双向文本拒绝、500 MB 输入体积守卫、
  输出目录的路径穿越约束，以及分类学表循环依赖检测——全部给出可执行的提示
  （`README.md` 的 "Input Validation"、[`SECURITY.md`](SECURITY_CN.md)、
  `test-data/README.md` 的 "boundary" 输入、`docs/TEST_RESULTS.md` §4.1 与 §5）。
- **PAML 版本感知**：4.8 与 4.10.x 两套控制文件方言，兼容性按主.次版本系列判定，
  跨系列使用会被拒绝而不是静默降级（`README.md` 的 "PAML 4.8 Support"）。
- **MCMC 收敛诊断**：多次运行的 MCMCTree 自动计算 Gelman–Rubin PSRF 与 ESS，并
  对照 1.05 / 100 阈值（`README.md` 的 "MCMC Convergence Diagnostics"）。
- **分发与治理**：MIT 许可证、中英双语的 README、QUICKSTART、CONTRIBUTING、
  CODE_OF_CONDUCT、SECURITY 与 THIRD_PARTY_LICENSES、`docs/api/` 下的 Sphinx
  文档源、教程、Conda 用的 `environment.yml`、四阶段
  `Dockerfile`（`builder`/`runtime`/`complete`/`tutorial`）与
  `docker-compose.yml`，以及 `CITATION.cff` / `CITATION.md`。
- **测试与基准设施**：收集 1,489 个测试（1,487 通过、2 跳过），覆盖 unit、
  functional、integration 与 release-validation 四层；26 项 `phylodater check`
  自检；覆盖率实测 67.87%，门槛 `fail_under = 60`；13 个分类单元 / 3 kb 参考
  基准及其归档记分卡（`test-data/benchmark/results/`）；以及运行 lint、文档策略、mypy
  与测试矩阵的 CI 工作流（`docs/TEST_RESULTS.md` §4-§5、`.github/workflows/ci.yml`）。

### 修复

0.1.0 的验证活动自身发现并修复了十个缺陷；每一个都在
[`docs/TEST_RESULTS.md`](docs/TEST_RESULTS_CN.md) §6 中记录了症状、修法与如今
钉住它的回归测试：

1. LSD2 默认使用氨基酸模型 `LG+G`，导致 IQ-TREE 2 在 DNA 上根本跑不起来；现在
   模型按比对类型解析。
2. MD-Cat 的 `--CI` 参数被当成一个逗号拼接的整串传入，上游无法解析；现在改为
   单个空白分隔的参数。
3. MD-Cat 的 CI 默认 100 次 bootstrap 重复（数小时的未告知工作）；现在默认
   `ci_nboots = 0`，与上游一致。
4. PAML 闸门死盯字面版本 `4.10.8`，拒绝了 bioconda 的 `paml` 4.10.10；现在按
   主.次版本系列判定兼容性。
5. `runtime_metadata.json` 不记录引擎版本；现在每个适配器都保存
   `validate_environment()` 探测到的结果。
6. 零依赖的 Newick 游走器拒绝开头的 `[&R]`，使 BEAST/IQ-TREE/DendroPy 导出的树
   无法处理；现在跳过并规范化前导注释。
7. PATHd8 丢弃非点式根校准（请求 [100, 110] Ma 却返回 153 Ma）；现在根的区间会
   以 `minage`/`maxage` 写出。
8. pyr8s 的区间校准既不执行也不报告；现在适配器在运行前披露 pyr8s 区间仅为建议，
   并在运行后逐条校验请求的区间。
9. 三个单元测试断言在进程级 `sys.modules` 上，纯粹因测试执行顺序而失败；现在改为
   探测一个全新的解释器。
10. 原始 MCMCTree 冒烟测试写出的 `BDparas` 行会被 PAML ≥ 4.10 拒绝；控制文件现在
    包含构造标志位。

### 已知限制（已记录，不隐藏）

- Docker 的 `complete` 阶段只内置 **7 个引擎中的 3 个**（经 PAML 4.10.8 的
  MCMCTree、经 IQ-TREE 2.3.4 的 LSD2、PATHd8）；treePL、r8s/pyr8s、wLogDate 与
  MD-Cat 需自行加装（`Dockerfile`、`docs/EXTERNAL_TOOLS.md` §7）。
- `treePL` 即便在 2.6.3 构建里也打印 `treePL version 1.0`，而 PATHd8 完全不打印
  版本，因此引擎溯源还依赖被记录下来的可执行文件路径
  （`docs/TEST_RESULTS.md` §7）。
- 未做跨平台的引擎性能测量；记录的墙钟时间与精度属于 `docs/TEST_RESULTS.md` §1
  的那台 Linux 主机和那份参考数据集，不可当作普适结论引用。
- `pyr8s` 只实现 NPRS，且把区间校准当作建议值处理
  （`docs/EXTERNAL_TOOLS.md` §5）。
