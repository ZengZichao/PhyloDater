> **Language:** English (canonical) · [中文](CONDA_PACKAGING_CN.md)

# Conda Packaging Guide for PhyloDater

This guide explains how to publish PhyloDater to conda (bioconda).

## option 1: Contributing via Bioconda

### 1. Prepare recipe

Create `recipes/phylodater/meta.yaml`:

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

### 2. Contribution Steps

1. Fork https://github.com/bioconda/bioconda-recipes
2. Add your recipe to `recipes/phylodater/`
3. Submit a Pull Request
4. Wait for review and merge

### 3. Notes

- Bioconda only accepts open source (OSI-approved) licenses
- Dependencies must be available in conda-forge or bioconda
- External software (PAML, IQ-TREE2, etc.) requires separate installation

## Option 2: conda-forge

### 1. Create a Feedstock

1. Visit https://github.com/conda-forge/phylodater-feedstock
2. Click "Fork"
3. Edit `recipe/meta.yaml`
4. Submit a PR

### 2. Automatic build

After merge, conda-forge will automatically:
- Build packages for multiple platforms
- Upload to conda-forge channel
- Update the build status

## Option 3: Personal conda channel

### 1. Using conda build

```bash
# Install conda build
conda install conda-build

# Build package
conda build recipes/phylodater -c conda-forge

# upload to channel
anaconda upload /path/to/build/package.tar.bz2
```

### 2. Share your channel

```bash
# Add channel
conda config --add channels your_channel_name

# Install
conda install phylodater
```

## Recommended Option

**Bioconda is strongly recommended**, for the following reasons:

1. Widely used in the biology community
2. Automatic build and testing
3. Robust dependency resolution
4. Comprehensive documentation: https://bioconda.github.io/contributor/index.html

## External Software Installation Suggestions

Since conda cannot directly bundle commercial/complex software, the following approaches are recommended:

### Option A: Install dependency software via conda

```bash
conda install -c bioconda paml iqtree r8s treepl pathd8
```

### Option B: Configure paths via environment variables

Users specify paths in the configuration file after installation:

```yaml
software_paths:
  # If the conda package ships a single PAML version, set paml_path only.
  # For coexisting versions, install them in separate directories and set paml_path accordingly.
  paml_path: /path/to/paml
  iqtree_bin: /path/to/iqtree2
  # When a same-named mcmctree exists in PATH, force a specific version with the full path:
  # mcmctree_bin: /path/to/paml/bin/mcmctree
```

## Testing conda builds

```bash
# Local test
conda install conda-verify
conda build recipes/phylodater --test

# Complete build workflow
conda build recipes/phylodater --output --python 3.10
```

## Supported Dating Methods

PhyloDater supports the following dating software; each one requires separate installation:

| Method | Dependency | Note |
|--------|-----------|------|
| MCMCTree | paml | Bayesian dating tool in the PAML package |
| LSD2 | iqtree | Fast dating method in IQ-TREE2 |
| r8s | r8s or pyr8s | Classical molecular clock dating tool |
| treePL | treepl | Penalized likelihood-based dating method |
| PATHd8 | pathd8 | Molecular clock calibration method |
| wLogDate | wlogdate | Log-date dating |
| MD-Cat | mdcat | Matrix-calibrated date |
