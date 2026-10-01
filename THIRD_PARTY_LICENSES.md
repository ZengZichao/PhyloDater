> **Language:** English (canonical) · [中文](THIRD_PARTY_LICENSES_CN.md)

# Third-Party Dependency Licenses

PhyloDater depends on the following third-party software. Each is distributed
under its own open-source license.

## Core dependencies (required)

| Library      | License (SPDX)                        | URL                                          |
|--------------|---------------------------------------|----------------------------------------------|
| BioPython    | BSD-3-Clause (Biopython License)      | https://biopython.org                        |
| DendroPy     | BSD 3-Clause                          | https://github.com/jeetsukumaran/DendroPy    |
| PyYAML       | MIT                                   | https://pyyaml.org                           |
| SciPy        | BSD 3-Clause                          | https://scipy.org                            |
| NumPy        | BSD 3-Clause                          | https://numpy.org                            |
| Pandas       | BSD 3-Clause                          | https://pandas.pydata.org                    |
| Matplotlib   | PSF-based (Matplotlib License)        | https://matplotlib.org                       |
| psutil       | BSD 3-Clause                          | https://github.com/giampaolo/psutil          |

## Optional dependencies

| Library      | License (SPDX)                        | URL                                          | Purpose                        |
|--------------|---------------------------------------|--------------|--------------------------------|
| **ETE3**     | **GPL-3.0-or-later**                  | https://etetoolkit.org                       | Tree visualisation (one plotting backend) |
| ETE4         | GPL-3.0-or-later                      | https://github.com/etetoolkit/ete            | Alternative ETE version        |

## License compatibility

**Important:** ETE3 is licensed under **GPLv3** (copyleft). When PhyloDater's
visualization module (`phylodater.viz`) imports ETE3, the combined work falls
under GPLv3 terms. To remain fully MIT-licensed:

- **Without ETE3 installed**: PhyloDater core remains MIT-licensed and fully
  functional for all dating analyses (visualization features disabled).
- **With ETE3 installed**: The `phylodater.viz` module that uses ETE3 must be
  treated as GPLv3-licensed.

All core dependencies (BioPython, DendroPy, NumPy, SciPy, Pandas, Matplotlib,
psutil, PyYAML) are BSD/MIT/PSF-licensed and are fully compatible with the MIT
license of the PhyloDater core.

Tree parsing uses a **zero-dependency built-in Newick walker by default**, so
ETE3 affects visualisation only, never core correctness; that invariant is
pinned by unit tests (`tests/unit/models/test_tree.py`,
`tests/unit/services/test_tree_rooting.py`).

## Third-party tool adapters (external software, not Python dependencies)

PhyloDater integrates the following external engines. They are NOT distributed
with PhyloDater and must be installed separately by the user:

| engine | licence | upstream |
|--------|---------|----------|
| PAML (MCMCTree, baseml, codeml) | free for non-commercial research and teaching; commercial use requires a licence | https://abacus.gene.ucl.ac.uk/software/paml.html |
| IQ-TREE 2 (LSD2 dating) | GPL-2.0-or-later | https://github.com/iqtree/iqtree2 |
| r8s | free for academic use | https://sourceforge.net/projects/r8s/ |
| pyr8s (Python implementation of r8s) | BSD-3-Clause (per upstream) | https://github.com/iTaxoTools/pyr8s |
| treePL | GPL-3.0-or-later (per upstream repository) | https://github.com/blackrim/treePL |
| PATHd8 | free for academic use | https://www2.math.su.se/PATHd8/ |
| wLogDate | GPL-3.0-or-later | https://pypi.org/project/wlogdate/ |
| MD-Cat | GPL-3.0 | https://github.com/uym2/MD-Cat |

The licenses of these external tools govern their own use and distribution. They
do not affect the license of PhyloDater source code. Where each engine was
obtained from, and the version actually measured in the release tests, is
documented in [docs/EXTERNAL_TOOLS.md](docs/EXTERNAL_TOOLS.md).

## Full license texts

For the complete license text of each third-party dependency, see the project's
`LICENSE` file and the deprecation information in `pyproject.toml`.

If you identify a licensing concern or incompatibility, please open an issue at
https://github.com/ZengZichao/phylodater/issues

## PhyloDater core license note

These notes previously lived at the bottom of the `LICENSE` file. They were
moved here so that `LICENSE` stays byte-identical to the standard MIT text
(which also lets GitHub's license detection recognize the repository as MIT).
The wording of the notes themselves is unchanged.

PhyloDater depends on third-party libraries, each distributed under their own
open-source licenses. Key license compatibility notes:

1. ETE3 (optional dependency for tree visualization) is licensed under GNU
   GPLv3. When ETE3 is imported at runtime, the visualization modules that
   directly invoke ETE3 APIs fall under GPLv3 terms. The core PhyloDater
   package (without ETE3) remains under MIT license and can be used,
   modified, and redistributed under MIT terms.

2. All core, non-optional dependencies (BioPython, DendroPy, NumPy, SciPy,
   Pandas, Matplotlib, PyYAML, psutil) are distributed under permissive
   BSD/MIT/PSF-compatible licenses and are fully compatible with the MIT
   license of the PhyloDater core.
