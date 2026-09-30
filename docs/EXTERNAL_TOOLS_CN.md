> **Language:** [English (canonical)](EXTERNAL_TOOLS.md) · 中文（本文档）

# 外部定年引擎

PhyloDater 是一个编排层：每一个 `dating` 方法都把真正的估算交给一个外部引擎。
本文记录了 [TEST_RESULTS.md](TEST_RESULTS.md) 中随包发布的测试结果所对应的
**每一个引擎的确切获取方式与验证过程**——定年结果的可复现性，取决于产出它的
那个二进制文件。

下列命令均在 conda 环境 `phylodater` 内、使用该环境自带的 `micromamba`
（`micromamba 2.9.0`）执行。各软件包的渠道可用性于 2026-09-29 针对在线渠道元数据
做过核对；对于**已经不再**通过 conda 或 PyPI 分发的引擎，本文直接写明这一点，
而不是保留一条注定失败的命令。

## 1. 引擎对照表

| `--method` | 实际执行的引擎 | 分发状态（2026-09-29） | 已验证的安装方式 |
|-----------|----------------|------------------------|------------------|
| `mcmctree` | PAML `mcmctree`（速率估算阶段还需 `baseml`） | bioconda 有 | `micromamba install -c bioconda paml`（4.10.10） |
| `lsd2` | IQ-TREE 2（`--date`，最小二乘定年） | bioconda 有 | `micromamba install -c bioconda iqtree`（3.1.3） |
| `wlogdate` | `logdate` 包提供的 `launch_wLogDate.py` | PyPI **与** bioconda 均有 | `pip install wlogdate`（1.0.2） |
| `r8s` | `pyr8s` Python API（仅 NPRS）或 `r8s` 二进制 | `pyr8s`：仅 GitHub（不在 PyPI、不在 bioconda）；`r8s`：仅 SourceForge/GitHub | `pip install git+https://github.com/iTaxoTools/pyr8s.git` |
| `mdcat` | MD-Cat 的 `md_cat.py` | 仅 GitHub（不在 PyPI、不在 bioconda） | `pip install git+https://github.com/uym2/MD-Cat.git` |
| `treepl` | `treePL` 二进制 | bioconda/conda-forge/PyPI **均已不再发布** | 源码构建，或复用历史构建产物 + `patchelf`（见 §3） |
| `pathd8` | `PATHd8` 二进制 | bioconda/conda-forge/PyPI **均已不再发布** | 从官方网站源码构建（见 §4） |

> **关于渠道。** `treepl`、`pathd8`、`r8s` 如今在 `bioconda`/`conda-forge` 上都无法
> 解析（`micromamba search` 返回 *no entries*，`https://api.anaconda.org/search?name=<pkg>`
> 也没有 bioconda 命中）。因此旧文档以及大量流水线中常见的
> `conda install -c bioconda treepl pathd8 r8s` 是**不能工作**的，不应把它当作受支持的
> 安装路径。

## 2. Python 层依赖

```bash
micromamba install -n phylodater -c conda-forge -c bioconda \
    python=3.12 paml iqtree nlopt=2.7.1 patchelf
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git
pip install git+https://github.com/uym2/MD-Cat.git
```

其中 `nlopt=2.7.1` 只有 `treePL` 二进制需要（它链接的是 `libnlopt.so.0`）；
NLopt ≥ 2.10 提供的是 `libnlopt.so.1`，无法满足该二进制。

## 3. treePL

`treePL` 是唯一一个如今无法用一条包管理命令完成安装的引擎。支持两种方式：

**方式 A —— 从上游源码构建（推荐）：**

```bash
git clone https://github.com/blackrim/treePL.git
cd treePL/deps && tar xzf nlopt-*.tar.gz   # 也可直接使用 §2 中 conda 安装的 nlopt
cd ../src && ./configure && make && make install
```

treePL 的优化过程依赖 NLopt；ADOL-C 是可选的（仅用于交叉验证并行化）。

**方式 B —— 复用最后一个 conda 构建产物并修复其库搜索路径：**

`treepl 2.6.3 (h3fd9d12_0)` 包把可执行文件放在 **prefix 根目录**，其
`RPATH = $ORIGIN/lib64:$ORIGIN/lib`；如果只是把它拷贝进 `<prefix>/bin`，库查找就会失败：

```
error while loading shared libraries: libnlopt.so.0
```

修复方法是把 `RPATH` 指向环境的 `lib` 目录：

```bash
patchelf --set-rpath '$ORIGIN/../lib' "$CONDA_PREFIX/bin/treepl"
```

在本测试环境中的核验结果：

| 项目 | 值 |
|------|-----|
| 可执行文件 | `$CONDA_PREFIX/bin/treepl` |
| SHA-256（patchelf 之前） | `50d2894a61c3f2be78c7f5bbfa9986f506496d854754d1b14986eafd6eaf5962` |
| 来源核验 | 与 conda 包 `treepl-2.6.3-h3fd9d12_0` 的 `info/paths.json` 中记录的 `treePL` 文件逐字节一致 |
| 自报版本号 | `treePL version 1.0`（上游即便在 2.6.3 构建中仍打印 1.0） |
| 依赖的共享库 | `libnlopt.so.0`（nlopt 2.7.1） |

PhyloDater 以 `treepl` 之名调用该二进制（`software.treepl.treepl_bin`，可用
`--method-args treepl_bin=/abs/path/treepl` 覆盖）。

## 4. PATHd8

PATHd8 以可移植 C 源码形式发布于官方网站 <https://www2.math.su.se/PATHd8/>；
整个程序是一个翻译单元，通过 `headers/` 里的文件被包含进来：

```bash
wget https://www2.math.su.se/PATHd8/PATHd8.zip && unzip PATHd8.zip
cc PATHd8.c -lm -O3 -o PATHd8
install -m 0755 PATHd8 "$CONDA_PREFIX/bin/"
```

该构建已在 gcc 11（Ubuntu 22.04）上验证：编译无警告，且构建出的二进制在基准数据集上
复现了环境中参考二进制的全部估算值（定年树完全一致，见
[TESTING_CN.md](TESTING_CN.md) §6）。

PATHd8 不打印版本横幅、也不接受 `--version`，因此 PhyloDater 改为记录实际解析到的
可执行文件路径；如需严格溯源，请自行固定版本。

## 5. r8s 与 pyr8s 的区别

当 `pyr8s` 可导入时 PhyloDater 优先使用它，否则回退到 `r8s` 二进制。引用结果前必须
理解其后果：

* `pyr8s` **只实现 NPRS**。它作为后端时，配置里的 `method`（`PL`、`LNL` 等）会被忽略，
  并发出一条 `EnvironmentWarning`。
* 原版 `r8s`（Sanderson 2003）发布于 <https://sourceforge.net/projects/r8s/>，并有镜像
  <https://github.com/R8S/r8s>。本仓库此前给出的地址
  （`phymaster.massey.ac.nz`）已无法解析。
* `pyr8s` 未发布于 PyPI 或任何 conda 渠道，只能从源码安装。

## 6. PhyloDater 自动记录的内容

现在每次运行的 `runtime_metadata.json` 都会写入 `validate_environment()` 探测到的
引擎标识：

```json
"software_versions": {"mcmctree": "MCMCTREE in paml version 4.10.10, 29 Jan 2026"}
```

随包发布的基准测试得到的实际记录：

| 方法 | 记录的版本字符串 |
|------|------------------|
| `mcmctree` | `MCMCTREE in paml version 4.10.10, 29 Jan 2026` |
| `lsd2` | `IQ-TREE version 3.1.3 for Linux x86 64-bit built Jul 26 2026` |
| `r8s` | `pyr8s 0.3.1 (NPRS only)` |
| `treepl` | `treePL version 1.0` |
| `wlogdate` | `wLogDate 1.0.2` |
| `mdcat` | `1.0.1` |
| `pathd8` | `PATHd8 (no version banner): /…/bin/PATHd8` |

## 7. Docker 镜像覆盖范围

[`Dockerfile`](../Dockerfile) **并未**提供全部七个引擎，它从上游压缩包构建：

| 引擎 | 是否在镜像内 | Dockerfile 使用的来源 |
|------|--------------|------------------------|
| MCMCTree | 是 | GitHub（`abacus-gene/paml`）的 PAML 4.10.8 发布压缩包 |
| LSD2 | 是 | IQ-TREE 2.3.4 Linux 压缩包 |
| PATHd8 | 是 | `www2.math.su.se` 的 `PATHd8.zip`，用 `gcc -O3` 编译 |
| treePL | 否 | — |
| r8s / pyr8s | 否 | — |
| wLogDate | 否 | — |
| MD-Cat | 否 | — |

因此该镜像只覆盖七个方法中的三个。若要在容器里使用其余引擎，请按 §2 与 §5 补齐
（镜像中已安装 `libnlopt0` 与 `libnlopt-dev`，可为后续添加 treePL 构建复用）。
每次修改 Dockerfile 都必须同步更新本表。

---

*另见：[TESTING_CN.md](TESTING_CN.md)（引擎如何被验证）、
[TEST_RESULTS_CN.md](TEST_RESULTS_CN.md)（实测结果）、
[English version](EXTERNAL_TOOLS.md)。*
