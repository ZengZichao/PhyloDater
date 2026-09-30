> **Language:** [English (canonical)](README.md) · 中文（本文档）

# examples

供快速入门（`QUICKSTART_CN.md`、`tutorials/01_quickstart_CN.ipynb`）使用的小型、自洽
输入：10 个分类单元、800 bp 比对、3 个校准——足够小，任何引擎的首次运行都能在几秒内完成。

| 文件 | 内容 |
|------|------|
| `tree.nwk` | 有根二叉 Newick，分支长度 = 严格钟下每位点期望替换数 |
| `alignment.fasta` | 在该树上按 JC69 仿真的 10 × 800 nt 比对 |
| `calibrations.yaml` | 固定年龄（`Homininae` 9 Ma）+ 均匀区间（`Primates` 70–80）+ 根均匀区间（100–110） |
| `taxonomy.tsv` | 独立示例：表格式分类学输入（`--taxonomy-file`），GTDB 风格字符串 |
| `benchmark/` | `scripts/tier1_manual_comparison.py` 与约束降级实验使用的额外输入 |
| `checksums.sha256` | 本目录各文件的 SHA-256 |

## 来源

树、比对与校准都是**生成的**，不是手写的：

```bash
python test-data/benchmark/generate_benchmark_dataset.py \
    --profile mini --out examples --force
python scripts/generate_checksums.py --recursive examples
```

它们与 [`../test-data/benchmark/`](../test-data/benchmark/) 的真值精度数据集出自同一个
仿真器，因此分支长度隐含的年龄是已知的（根 105 Ma、`Homininae` 9 Ma、`Primates` 75 Ma）。
`taxonomy.tsv` 与上面的树无关，它只演示 `--taxonomy-file` 的输入格式，其中的分类名故意采用
accession 风格。

有两个性质对首次运行很关键，而生成器保证了它们：

* **树是二叉的**——多叉树会让若干引擎拒绝运行或发生退化；
* **至少有一个校准是固定年龄**——没有任何 `fixage` 点时 PATHd8 会拒绝运行，而
  MD-Cat/wLogDate 本来就只能表达点校准。

## 已验证

在参考环境中，七个方法用这些输入都能成功完成（引擎及其来源见
[docs/TEST_RESULTS_CN.md](../docs/TEST_RESULTS_CN.md) §1.1）：

```bash
for m in mcmctree lsd2 r8s treepl pathd8 wlogdate mdcat; do
  phylodater dating -t examples/tree.nwk -s examples/alignment.fasta \
      -c examples/calibrations.yaml -o out_$m --method "$m" --seed 42
done
```

这些文件是示例，**不是**精度参照：定量评分使用
[`../test-data/benchmark/`](../test-data/benchmark/)，那里附带真实的节点年龄。
