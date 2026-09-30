> **Language:** [English (canonical)](README.md) · 中文（本文档）

# test-data

随包发布的测试数据，目的是让 [docs/TEST_RESULTS_CN.md](../docs/TEST_RESULTS_CN.md)
记录的验证工作可以被重新运行（或重新核查），无需重建任何东西。

本目录包含两类数据，用途不同：

| 路径 | 角色 | 被谁使用 |
|------|------|----------|
| `benchmark/` | **定量**：已知分歧时间的仿真数据集，用于给每个方法的估计值与真值的接近程度打分 | `scripts/run_reference_benchmark.py`、`tests/validation/integration/` |
| `input/` | **定性**：合法 / 非法 / 边界场景输入，用于手工端到端核查工作流及其错误处理 | QA 走查、`tests/functional/` 的 fixture |

本目录内不含自动化测试；测试位于 [`tests/`](../tests/)。

## `benchmark/` —— 带真值的精度数据集

| 文件 | 内容 |
|------|------|
| `tree.nwk` | 有根二叉树，13 个分类单元，分支长度 = 严格钟下的每位点期望替换数 |
| `tree_ml.nwk` | 同一拓扑，但分支长度由比对经**最大似然**估计——即真实用户会提供的输入 |
| `alignment.fasta` | 在上述树上按 JC69 仿真的 13 × 3000 nt 比对 |
| `calibrations.yaml` | 位于四个互不相同节点上的四个校准（一个点校准、三个均匀区间） |
| `truth.yaml` | 生成参数（种子、速率、模型、长度）+ 全部 12 个内部节点的真实年龄，以 MRCA 叶对为键 |
| `generate_benchmark_dataset.py` | 确定性生成脚本；`--include-ml-tree` 会用 IQ-TREE 2 重新生成 `tree_ml.nwk` |
| `results/` | 发布运行的归档结果：`scorecard.tsv`、`node_ages.tsv`、`benchmark_results.json`，以及 ML 变体的 `ml_tree/` |
| `results/PYTHON_SUPPORT_MATRIX.txt` | 已声明 Python 支持矩阵背后的证据：逐版本的依赖解析结果，以及在 3.8–3.14 上实跑的 `tests/unit` 结论 |
| `results/validation_log.txt` | 发布门禁命令集的原样终端输出（pytest、`phylodater check`、格式化器、`check_docs.py`、mypy、两个基准变体） |
| `checksums.sha256` | 上述所有文件的 SHA-256（`results/runs/` 下的单次运行产物除外） |

为什么用仿真数据：真实化石校准不能充当真值，因为那样等于用方法被告知的数字去给方法
打分。严格钟 + JC69 仿真为每个估计器提供了明确的参照。拓扑与年龄是合成的但合理，它们
是数值参照，不构成任何演化史结论。完整设计理由与验收门槛见
[docs/TESTING_CN.md](../docs/TESTING_CN.md) §3–§4。

重新生成并重新评分：

```bash
python test-data/benchmark/generate_benchmark_dataset.py --force --include-ml-tree
sha256sum -c test-data/benchmark/checksums.sha256
python scripts/run_reference_benchmark.py --output test-data/benchmark/results
python scripts/run_reference_benchmark.py --tree-file tree_ml.nwk \
    --output test-data/benchmark/results/ml_tree
```

`results/runs/`（各方法的定年树、比较表、日志）按需重新生成并被 git 忽略；随包发布的是
汇总文件与上述两份证据文件。

`results/validation_log.txt` 不是基准工装产出的，而是 [docs/TESTING_CN.md](../docs/TESTING_CN.md)
§5 列出的发布门禁命令在验证机上一次性捕获的终端输出，原文不改动地提交（文件头记录了
它所验证的 `git HEAD`；日志文件本身及其校验和条目必然由紧随其后的提交加入）。每次发布
前都按同样方式重新捕获。

## `input/` —— 场景输入

| 路径 | 用途 |
|------|------|
| `input/valid/` | 标准工作流输入（树、比对、校准配置） |
| `input/invalid/` | 必须触发错误处理的输入（非法字符、畸形格式、重复标签等） |
| `input/boundary/` | 边界条件（空树、多棵树、BOM、GBK 编码等） |

使用方式：从对应子目录选取文件，运行 `phylodater dating`，再把观察到的行为（退出码、
日志、披露的告警）与预期比对。

## 与自动化测试的关系

- `tests/` —— pytest 驱动的自动化测试（unit / functional / integration / validation），在 CI 中运行
- `test-data/`（本目录）—— 上述测试与基准工装所读取的数据

两者互补：数据负责科学内容是否站得住，测试负责行为是否正确。
