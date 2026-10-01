> **Language:** [English (canonical)](THIRD_PARTY_LICENSES.md) · 中文（本文档）

# 第三方依赖许可

PhyloDater 依赖下列第三方软件，各自按其自己的开源许可分发。

## 核心依赖（必需）

| 库 | 许可证（SPDX） | 网址 |
|----|----------------|------|
| BioPython | BSD-3-Clause（Biopython License） | https://biopython.org |
| DendroPy | BSD-3-Clause | https://github.com/jeetsukumaran/DendroPy |
| PyYAML | MIT | https://pyyaml.org |
| SciPy | BSD-3-Clause | https://scipy.org |
| NumPy | BSD-3-Clause | https://numpy.org |
| Pandas | BSD-3-Clause | https://pandas.pydata.org |
| Matplotlib | PSF 系（Matplotlib License） | https://matplotlib.org |
| psutil | BSD-3-Clause | https://github.com/giampaolo/psutil |

## 可选依赖

| 库 | 许可证（SPDX） | 网址 | 用途 |
|----|----------------|------|------|
| **ETE3** | **GPL-3.0-or-later** | https://etetoolkit.org | 树可视化（绘图后端之一） |
| ETE4 | GPL-3.0-or-later | https://github.com/etetoolkit/ete | ETE 的替代版本 |

## 许可兼容性

**重要：** ETE3 采用 **GPLv3**（copyleft）许可。当 PhyloDater 的可视化模块
（`phylodater.viz`）导入 ETE3 时，组合作品将落入 GPLv3 条款之下。要始终保持 MIT
许可，请按下述方式理解：

- **未安装 ETE3**：PhyloDater 核心仍为 MIT 许可，全部定年分析功能完整可用
  （可视化功能不可用）。
- **已安装 ETE3**：使用 ETE3 的 `phylodater.viz` 模块必须按 GPLv3 许可对待。

所有核心依赖（BioPython、DendroPy、NumPy、SciPy、pandas、Matplotlib、psutil、
PyYAML）均为 BSD/MIT/PSF 许可，与 PhyloDater 核心的 MIT 许可完全兼容。

树解析的**默认后端是零依赖的内置 Newick 游走器**，因此 ETE3 只影响可视化，不影响
核心正确性；相关约束由单元测试固定（见 `tests/unit/models/test_tree.py`、
`tests/unit/services/test_tree_rooting.py`）。

## 外部工具适配器（外部软件，不是 Python 依赖）

PhyloDater 集成下列外部命令行工具/引擎。它们**不**随 PhyloDater 分发，必须由用户
自行安装：

| 引擎 | 许可 | 上游 |
|------|------|------|
| PAML（MCMCTree、baseml、codeml） | 免费用于非商业研究与教学，商用需授权 | https://abacus.gene.ucl.ac.uk/software/paml.html |
| IQ-TREE 2（LSD2 定年） | GPL-2.0-or-later | https://github.com/iqtree/iqtree2 |
| r8s | 免费用于学术 | https://sourceforge.net/projects/r8s/ |
| pyr8s（r8s 的 Python 实现） | BSD-3-Clause（上游声明为准） | https://github.com/iTaxoTools/pyr8s |
| treePL | GPL-3.0-or-later（上游仓库声明） | https://github.com/blackrim/treePL |
| PATHd8 | 免费用于学术 | https://www2.math.su.se/PATHd8/ |
| wLogDate | GPL-3.0-or-later | https://pypi.org/project/wlogdate/ |
| MD-Cat | GPL-3.0 | https://github.com/uym2/MD-Cat |

这些外部工具的许可只约束其自身的安装与使用，不影响 PhyloDater 源代码的许可。
各工具的**获取途径与实测版本**见 [docs/EXTERNAL_TOOLS_CN.md](docs/EXTERNAL_TOOLS_CN.md)。

## 完整许可文本

各第三方依赖的完整许可文本请查看其发行包自带的 `LICENSE` 文件，以及本项目的
`LICENSE`。

如果您发现任何许可方面的问题或冲突，请在
https://github.com/ZengZichao/phylodater/issues 提交 issue。

## PhyloDater 核心许可说明

以下说明原先附于 `LICENSE` 文件末尾，为保证 `LICENSE` 与标准 MIT 文本完全一致
（也便于 GitHub 的许可证识别）而迁移至此。说明内容本身未做任何修改。

PhyloDater 依赖的第三方库各自按其开源许可分发。关键的许可兼容性说明：

1. ETE3（树可视化的可选依赖）采用 GNU GPLv3 许可。当运行时导入 ETE3 时，
   直接调用 ETE3 API 的可视化模块随之适用 GPLv3 条款。不含 ETE3 的
   PhyloDater 核心包仍为 MIT 许可，可按 MIT 条款使用、修改与再分发。

2. 所有核心、非可选依赖（BioPython、DendroPy、NumPy、SciPy、Pandas、
   Matplotlib、PyYAML、psutil）均以宽松的 BSD/MIT/PSF 兼容许可分发，
   与 PhyloDater 核心的 MIT 许可完全兼容。
