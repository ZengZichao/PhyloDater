> **Language:** [English (canonical)](TESTING.md) · 中文（本文档）

# 测试指南

本文说明 **PhyloDater 是如何被测试的**、每一项测试为何这样设计，以及如何复现
[TEST_RESULTS_CN.md](TEST_RESULTS_CN.md) 中引用的全部结果。本文属于发布内容的一部分：
测试数据、测试工装（harness）、预期判据与实测结果一并随包发布，评审者无需自行重建
即可重跑整套验证。

范围：对全部七个受支持定年方法的功能验证，以及保护这些方法的单元 / 集成 / 安全测试层。

---

## 1. "受支持"意味着什么，以及必须测什么

PhyloDater 声称支持七个 `dating` 方法。每个方法都是不同引擎之上的一层标准化外壳，
因此测试矩阵必须覆盖三个互相独立的层面：

1. **外壳本身** —— PhyloDater 是否能构造出合法的引擎输入、成功调用引擎，并把输出
   解析进统一的结果模型？
2. **语义** —— 校准约束在被翻译成引擎自己的语法之后是否仍然保持原意（以及每一次语义
   损失是否都被报告）？
3. **估算值** —— 得到的时间树在**数值上**是否接近真值？

第 1、2 点可以用 mock 检查；第 3 点不行，它需要一份已知分歧时间的数据集，这正是本包
附带一套"仿真真值"基准（§3）的原因。七个方法中的每一个都针对**真实引擎**执行，
没有任何一个方法只通过桩（stub）来测试。

## 2. 测试分层

```
tests/
├── unit/          29 个文件 - 纯 Python：模型、约束、配置、服务、适配器逻辑
├── functional/    12 个文件 - CLI 行为、输入校验、报告、安全、可视化
├── integration/    4 个文件 - 多步流水线、绘图、引擎原始调用冒烟测试
└── validation/    29 个文件 - 发布验证套件（真实引擎、端到端）
    ├── unit/ adapters/ functional/ integration/ security/
    ├── fixtures/        静态输入（树、比对、校准、配置）
    ├── helpers.py       引擎可用性探测 + CLI 运行器
    └── conftest.py      marker：slow、requires_mcmctree、requires_treepl 等
```

| 层 | 是否调用外部软件 | 默认 `pytest tests` 运行时 |
|----|------------------|-----------------------------|
| `unit`、`functional` | 从不调用 | 始终执行 |
| `integration` | 仅"引擎原始调用"测试需要 | 二进制存在则执行，否则跳过 |
| `validation` | 真实端到端运行 | 引擎已安装则执行，否则跳过并给出明确原因 |

需要引擎的测试必须显式声明，于是一台没有 PAML 的机器上报出的是 *skipped*，
而绝不是 *passed*：

```python
@pytest.mark.slow
@pytest.mark.requires_mcmctree
class TestMCMCTreeBenchmark: ...
```

`helpers.external_available()` 把方法名映射到真正需要存在的东西：
`lsd2 → iqtree2/iqtree`、`r8s → pyr8s 包或 r8s 二进制`、
`wlogdate → logdate 包或 launch_wLogDate.py`、`mdcat → emd 包或 md_cat.py`、
`pathd8 → PATHd8`、`treepl → treePL`、`mcmctree → mcmctree + baseml`。

## 3. 基准数据集

`test-data/benchmark/` 是一份**仿真**的、已知答案的数据集：

| 文件 | 内容 |
|------|------|
| `tree.nwk` | 有根二叉树，13 个分类单元，分支长度 = 每位点替换数 |
| `alignment.fasta` | 13 × 3000 nt 多序列比对 |
| `calibrations.yaml` | 位于四个互不相同节点上的四个校准 |
| `truth.yaml` | 生成参数 + 全部 12 个内部节点的真实年龄 |
| `generate_benchmark_dataset.py` | 生成脚本（确定性、仅依赖 numpy） |
| `results/` | `scorecard.tsv`、`node_ages.tsv`、`benchmark_results.json` |

设计决策，以及为什么它们才是科学上站得住脚的选择：

* **严格分子钟 + JC69 + 碱基频率相等。** 这样每个方法都在"自己的假设成立"的条件下被
  评估，偏差可以归因于估计方法本身，而不是模型违背。速率自相关、谱系间速率变异或
  不可表达的替换模型都会让"离真值多远"变得无法解释。
* **不能用带化石校准的真实数据当作真值。** 那样等于用方法被告知的校准去给方法打分，
  是循环论证。仿真数据集打破了这一循环；真实数据示例保留在 `examples/`。
* **拓扑与年龄是合成的，但合理**（有胎盘哺乳动物框架：13 个分类单元，根龄 105 Ma）。
  它们是数值参照，不构成任何演化史结论。
* **速率 `1.5e-3` 替换/位点/Myr，比对长度 3000**，根到叶约 470 个期望差异位点：既有
  足够信号支撑 ML 分支长度估计，最深节点又远未饱和。
* **随机种子 `20260929`，仅用标准库 + numpy。** 重跑生成脚本会得到逐字节一致的文件，
  因此数据集的哈希始终有效：

  ```bash
  python test-data/benchmark/generate_benchmark_dataset.py --force
  sha256sum test-data/benchmark/{tree.nwk,alignment.fasta,truth.yaml,calibrations.yaml}
  ```

* **只使用点校准与双边区间校准。** `Homininae` 固定 9 Ma、`Primates` uniform 68–82、
  `Glires` uniform 60–80、`Root` uniform 100–110（均以真值为中心）。单边约束
  （`maximum`、`soft_lower`）**不在此数据集中使用**：只接受点校准的引擎（wLogDate、
  MD-Cat）会把它降级为隐式区间 `[0, max]` 的中点，这可能与其余校准构成时间上矛盾的
  约束组，并直接在引擎内部中止。这类约束类型改由
  `examples/benchmark/degradation/` 下的"约束覆盖"数据集专门验证——在那里，
  **降级警告本身**才是被测试的对象。
* **任意两个校准不会落到同一个节点**，也没有内部校准落在根节点上：MCMCTree 把根龄写进
  控制文件（`RootAge`）而非写入树中；当它发现某个请求的校准在 `mcmctree.tree` 中缺失时，
  会直接拒绝运行而不是带着一套被悄悄削弱过的约束继续。

## 4. 评分方法

`scripts/run_reference_benchmark.py` 为每个方法执行 CLI，然后按下列顺序检查：

| # | 检查项 | 判定方式 |
|---|--------|----------|
| 1 | 完成 | CLI 退出码为 0 **且** `runtime_metadata.json` 中 `status: completed` |
| 2 | 叶节点完整 | 输出树的分类单元集合与输入完全一致，无重命名、无丢失 |
| 3 | 拓扑一致 | 输出的内部演化支集合（以叶节点集合表示）= 输入；采用演化支匹配，因此节点*标签*不参与判定 |
| 4 | 年龄非负 | 每个内部节点年龄 ≥ 0 |
| 5 | 时序单调 | 不存在子节点比其父节点更老 |
| 6 | 校准被满足 | 点校准：\|估计值 − 年龄\| ≤ 1 Ma；区间校准：落在 `[min, max]` 内外扩区间宽度的 5 % |
| 7 | 精度 | 对全部 12 个内部节点，与 `truth.yaml` 比较的 MAE、RMSE、最大绝对误差、最大相对误差 |

随包基准使用的门槛是：**MAE ≤ 10 Ma**，且**最差不超过真值 0.5 倍相对误差**。
门槛设置得足以吸收 3 kb 比对的统计噪声，又足够严格，以至于会"静默错标数量级"的方法
（例如下方第 8 节提到的单位问题）必然被拦下。

如何从各引擎输出中还原节点年龄——这一点至关重要，读错了整套基准就作废：

* **时间树解码。** 定年后的树把时间写成分支**差值**，因此
  `age(node) = 根到叶最大距离 − 该节点到根的距离`。若直接使用"到根的距离"，得到的是
  "自根以来流逝的时间"——根为 0、最年轻的节点最大，恰好是把年龄完全颠倒。
* **单位推断。** MCMCTree 输出 Ga，其余引擎输出 Ma。换算系数由根校准推断
  （在 1、10⁻³、10³、10⁶ 中取最接近者），而不是按方法硬编码；因此将来若某个引擎用了
  别的单位，会**显眼地失败**，而不是静默出错。
* **演化支匹配。** 节点靠其后代叶集合识别，而不是靠节点名：多个引擎会剥离或改写内部
  标签（`IntNode7`、`I9`，或者根本没有标签）。
* **注解剥离。** 重新解析前会移除 `[&date=…]`、`[t=…, mu=…]`、`[&95%HPD={…}]` 以及内部
  标签，并同时兼容 DendroPy 4.x 与 5.x 的标签访问接口。

## 5. 复现测试

```bash
# 1. 环境（七个引擎齐备）——见 EXTERNAL_TOOLS_CN.md
conda activate phylodater

# 2. Python 层自检（26 项，无需外部引擎）
phylodater check

# 3. 全部自动化测试（共收集 1 489 项：1 487 passed、2 skipped；
#    引擎缺失时依赖它的测试会跳过）
python -m pytest tests -q

# 4. 仅发布验证套件（真实引擎、端到端）
python -m pytest tests/validation -q

# 5. 带真值评分的基准测试（重新生成 results/）
python scripts/run_reference_benchmark.py \
    --dataset test-data/benchmark --output test-data/benchmark/results
```

常用开关：`--methods` 选择子集；`--method-args METHOD=key=value,...` 覆盖某个方法的
基准参数（下文的 MD-Cat 敏感性检查就是这么做的）；`--timeout` 单方法超时；
`--seed`（默认 42）传给 PhyloDater；`-m "not slow"` 排除依赖引擎的测试。

参考机器（48 物理核）上的耗时：七个方法的完整基准约 2 分钟，其中 treePL 约 1 分钟
（交叉验证搜索），MCMCTree 在使用下文缩减版 MCMC 参数时约 15 秒。

发布门禁不只有测试套件：`black --check`、`isort --check`、`ruff check`、
`python scripts/check_docs.py` 与 `mypy phylodater`（开启 `disallow_untyped_defs` 与
`warn_return_any` 前提下 0 error）都必须通过，两个基准变体（参考树与 ML 枝长树）
也都必须 7/7 PASS。某次全部通过的完整原样输出已归档到
[`test-data/benchmark/results/validation_log.txt`](../test-data/benchmark/results/validation_log.txt)。

## 6. 基准使用的引擎参数

为了让矩阵在几分钟内跑完，部分方法的参数低于工具默认值。所有缩减都在此处与
评分卡中记录，绝不隐藏：

| 方法 | 实际使用参数 | 默认参数 | 原因 |
|------|--------------|----------|------|
| `mcmctree` | `clock=2,num_runs=2,burnin=2000,nsample=4000,sampfreq=10` | `clock=2,num_runs=2,burnin=20000,nsample=50000,sampfreq=5` | 13 分类单元 / 3 kb 不需要 5 万次采样；收敛是**被测量的**，不是被假设的（PSRF ≤ 1.05，ESS ≥ 100） |
| `mdcat` | `ncat=10,nrep=20,max_iter=50` | `ncat=50,nrep=100,max_iter=100` | EM 随机重启主导运行时间。**实测**：在本数据集上缩减参数与默认参数给出的 12 个节点估计**完全一致**（MAE 均为 0.081 Ma，10 秒 vs 46 秒），即此处缩减不损失精度 |
| `lsd2`、`r8s`、`treepl`、`pathd8`、`wlogdate` | PhyloDater 默认值 | — | 无需缩减 |

`mdcat.ci_nboots` 保持默认 `0`（不计算置信区间），与上游一致：在本数据集上**单次**
bootstrap 就要数分钟，默认 100 次会让该方法不可用。

PATHd8 交叉核验：按上游源码构建的二进制（Dockerfile 的配方，
`gcc -O3 PATHd8.c -lm`）与环境中已有参考二进制产出的定年树**完全一致**——这样被验证的
就不只是 PhyloDater 对引擎的调用，还包括引擎本身。

## 7. 确定性

对同一 `--seed` 的两次独立基准运行逐节点比较：

* `lsd2`、`r8s`、`treepl`、`pathd8`、`wlogdate`、`mdcat`：两次运行的节点估计**逐位一致**。
* `mcmctree`：两次运行相差 ≤ 0.3 Ma（不足根龄的 0.3 %）。这是贝叶斯采样器的正常表现，
  固定 PRNG 种子并不能让 PAML 的 MCMC 达到逐位可复现。因此评判 MCMCTree 质量的正确
  闸门是 PhyloDater 记录的收敛诊断量（`max PSRF`、`min ESS`），而不是运行间完全相等。

## 8. 语义降级的报告

约束翻译是这类定年封装最容易"静默造成伤害"的地方，因此每一次降级都必须
(a) 在日志中可见，(b) 被写入结果对象，(c) 有单元测试覆盖。验证套件检查的例子：

| 输入约束 | 引擎 | 实际发生 | 警告 |
|----------|------|----------|------|
| `soft_lower 60` | lsd2 | 硬单边界 `NA:-60` | `SemanticDegradationWarning` |
| `maximum 130` | wlogdate / mdcat | 点估计 65 Ma（`[0,130]` 中点） | `SemanticDegradationWarning`，且明确"上界本身不会被强制" |
| `uniform 8–10` | mdcat | 中点 9 Ma | `SemanticDegradationWarning` |
| `gamma` / `skew_normal` / `skew_t` | r8s、pathd8、treepl、lsd2 | 正态近似区间 | `SemanticDegradationWarning`，并写明近似方式 |
| 任意软边界 | mcmctree | 精确的 `B()/RN()/ST()` 先验 | 无（忠实表达） |

若校准是被**丢失**而非降级，方法还会 fail-closed：MCMCTree 适配器会把请求的约束与
`mcmctree.tree` 中实际存在的约束逐条比对，缺少任何一条即中止运行，绝不在约束被悄悄
削弱的情况下继续计算。

## 9. 与已修复缺陷绑定的回归测试

本轮测试工作发现的每个缺陷都由对应测试"钉住"，避免问题复发：

| 缺陷 | 修复前的表现 | 钉住它的测试 |
|------|--------------|--------------|
| LSD2 默认替换模型是氨基酸模型（`LG+G`），但流程针对核酸比对 | IQ-TREE 2 报 `ERROR: File not found LG`，方法完全无法运行 | `TestLSD2ModelResolution`（5 个用例） |
| MD-Cat 的 `--CI` 以逗号拼接的单个串传入 | 上游 `int("100,0.025,0.975")` 抛 `ValueError`；默认开启的 CI 让每次运行都失败 | `TestMDCatCliContract` |
| MD-Cat CI 默认 100 次 bootstrap | 数小时的隐式运行开销 | `TestMDCatConfig::test_ci_parameters` |
| PAML 版本闸门要求逐字相等（`4.10.8`） | 按文档 `conda install -c bioconda paml` 安装的 4.10.10 被直接拒跑 | `test_validate_environment_same_series_other_patch_is_accepted`、`test_paml_versions_compatible_series_semantics` |
| `runtime_metadata.json` 从不记录引擎版本 | 随包发布的结果缺少引擎溯源信息 | 适配器版本记录 + 基准评分卡 |
| 三个测试断言的是**进程级** `sys.modules` | 它们只因测试顺序（可视化测试导入了 ete3）而失败，与产品行为无关 | `module_imports_cleanly` fixture（子进程探测） |
| MCMCTree 原始调用冒烟测试写出 `BDparas = 1 1 0.1` | PAML ≥ 4.10 中止："BDparas: expect flag for birth-death process prior" | `tests/integration/full_pipeline/test_all_dating_tools.py` 中已修正的控制文件串 |

## 10. 已知局限（明示，不隐藏）

* **Docker 镜像只覆盖七个引擎中的三个**（MCMCTree、LSD2、PATHd8），详见
  [EXTERNAL_TOOLS_CN.md](EXTERNAL_TOOLS_CN.md) §7。
* **`r8s` 经由 `pyr8s` 运行，仅支持 NPRS**；`method=PL` 会被忽略并发出
  `EnvironmentWarning`。
* **treePL 目前没有 conda/PyPI 包**，其安装与环境强相关；二进制缺失时测试会干净地跳过。
* **置信区间由各引擎自行定义。** `comparison_plot.png` / `comparison_table.tsv`
  对没有区间的单元格做显式标记，而不是画成长度为 0 的误差条。
* **精度数字只描述这份数据集。** [TEST_RESULTS_CN.md](TEST_RESULTS_CN.md) 中的 MAE
  对应这份 13 分类单元、3 kb、严格钟仿真数据，**不是**对这些方法通用精度水平的断言。

---

*另见：[TEST_RESULTS_CN.md](TEST_RESULTS_CN.md)、
[EXTERNAL_TOOLS_CN.md](EXTERNAL_TOOLS_CN.md)、[English version](TESTING.md)。*
