> **Language:** [English (canonical)](CONDA_PACKAGING.md) · 中文（本文档）

# PhyloDater Conda 打包指南

本指南说明如何将 PhyloDater 发布到 conda（bioconda）。

## 方案一：通过 Bioconda 投稿

### 1. 准备 recipe

创建 `recipes/phylodater/meta.yaml`：

```yaml
package:
  name: phylodater
  version: "0.1.0"

source:
  path: ..

build:
  number: 0
  script: python -m pip install . -vv

requirements:
  host:
    - python >=3.9
    - pip
  run:
    - python >=3.9
    - biopython >=1.79
    - ete3 >=3.1.2
    - ete4 >=3.0.0
    - dendropy >=4.5.0
    - pyyaml >=6.0
    - scipy >=1.9.0
    - numpy >=1.21.0
    - psutil >=5.9.0
    - matplotlib >=3.5.0
    - pandas >=1.3.0

test:
  imports:
    - phylodater
  requires:
    - pytest

about:
  home: https://github.com/yourusername/phylodater
  license: MIT
  license_family: MIT
  summary: Multi-software parallel platform for phylogenetic molecular dating
  description: |
    PhyloDater is a unified phylogenetic molecular dating platform
    that integrates multiple dating software tools.
```

### 2. 投稿步骤

1. Fork https://github.com/bioconda/bioconda-recipes
2. 添加你的 recipe 到 `recipes/phylodater/`
3. 提交 Pull Request
4. 等待 review 和合并

### 3. 注意事项

- Bioconda 只接受开源（OSI approved）许可证
- 依赖必须在 conda-forge 或 bioconda 中可用
- 外部软件（PAML、IQ-TREE2 等）需要单独安装

## 方案二：通过 conda-forge

### 1. 创建 feedstock

1. 访问 https://github.com/conda-forge/phylodater-feedstock
2. 点击 "Fork"
3. 修改 `recipe/meta.yaml`
4. 提交 PR

### 2. 自动构建

合并后，conda-forge 会自动：
- 为多个平台构建包
- 上传到 conda-forge channel
- 更新构建状态

## 方案三：个人 conda channel

### 1. 使用 conda build

```bash
# 安装 conda build
conda install conda-build

# 构建包
conda build recipes/phylodater -c conda-forge

# 上传到 channel
anaconda upload /path/to/build/package.tar.bz2
```

### 2. 分享你的 channel

```bash
# 添加 channel
conda config --add channels your_channel_name

# 安装
conda install phylodater
```

## 推荐方案

**强烈推荐使用 Bioconda**，原因：

1. 生物学社区广泛使用
2. 自动构建和测试
3. 依赖解析完善
4. 文档齐全：https://bioconda.github.io/contributor/index.html

## 外部软件安装建议

conda 不能直接包含商业软件，也不能包含结构复杂的软件，建议：

### 方案 A：通过 conda 安装依赖软件

```bash
conda install -c bioconda paml iqtree r8s treepl pathd8
```

### 方案 B：通过环境变量配置路径

用户安装后在配置文件中指定：

```yaml
software_paths:
  # 若 conda 包仅包含一个 PAML 版本，指定 paml_path 即可；
  # 多版本共存时建议分目录安装并分别指定 paml_path
  paml_path: /path/to/paml
  iqtree_bin: /path/to/iqtree2
  # 当 PATH 中存在同名 mcmctree 时，使用完整路径强制选择版本
  # mcmctree_bin: /path/to/paml/bin/mcmctree
```

## 测试 conda 构建

```bash
# 本地测试
conda install conda-verify
conda build recipes/phylodater --test

# 完整构建流程
conda build recipes/phylodater --output --python 3.10
```

## 支持的定年方法

PhyloDater 支持以下定年软件，每种软件都需要单独安装：

| 方法 | 依赖包 | 说明 |
|------|--------|------|
| MCMCTree | paml | PAML 包中的贝叶斯定年工具 |
| LSD2 | iqtree | IQ-TREE2 中的快速定年方法 |
| r8s | r8s 或 pyr8s | 经典分子钟定年工具 |
| treePL | treepl | 基于惩罚似然的定年方法 |
| PATHd8 | pathd8 | 分子钟校准方法 |
| wLogDate | wlogdate | 日志日期定年 |
| MD-Cat | mdcat | 矩阵校准日期 |
