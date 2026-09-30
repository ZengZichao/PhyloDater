> **Language:** English (canonical) · [中文](CONTRIBUTING_CN.md)

# Contributing Guide

> 📚 This document is also available in [中文 (Chinese)](CONTRIBUTING.md).
> 本文档的权威中文版见 [CONTRIBUTING.md](CONTRIBUTING.md)。
> **Last Updated**: 2026-04-09

Thank you for considering contributing code to PhyloDater!

## How to Contribute

### Reporting Issues

If you find a bug or have a new feature request, please:

1. Search existing issues to confirm the problem has not already been reported
2. Create a new issue containing:
   - A clear title and description
   - Steps to reproduce
   - Expected behavior vs actual behavior
   - Environment information (Python version, operating system, etc.)

### Submitting Code

1. Fork this repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Develop and commit your changes
4. Ensure tests pass: `pytest tests/`
5. Submit a Pull Request

### Code Style

- Follow PEP 8
- Use type annotations
- Add docstrings
- Write unit tests

## Development Environment Setup

### 1. Clone Repository

```bash
git clone https://github.com/yourusername/phylodater.git
cd phylodater
```

### 2. Create Development Environment

```bash
conda create -n phylodater-dev python=3.12
conda activate phylodater-dev

pip install -e ".[dev]"

conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git   # r8s backend (NPRS only)
pip install git+https://github.com/uym2/MD-Cat.git        # MD-Cat
# treePL and PATHd8: build from source (docs/EXTERNAL_TOOLS.md sections 3-4);
# neither is published on bioconda/conda-forge/PyPI any more.
```

### 3. Run Tests

```bash
python -m pytest tests -q                 # developer + release-validation suites
python -m pytest tests -m "not slow" -q   # skip tests that run real engines
python -m pytest tests/validation -q      # release validation only
phylodater check                          # 26-item environment self-test
```

Engine-dependent tests skip themselves (with a reason) when the binary is absent, so
the suite stays runnable on a machine that has none of the seven engines installed.

Adding a dating method means adding a ground-truth run for it too:
`python scripts/run_reference_benchmark.py` scores every method against the
simulated truth in `test-data/benchmark/` — see
[docs/TESTING.md](docs/TESTING.md) for the gates and
[docs/TEST_RESULTS.md](docs/TEST_RESULTS.md) for the current numbers.

### 4. Code formatting, linting and type checking

```bash
black phylodater tests scripts test-data
isort  phylodater tests scripts test-data
ruff check phylodater tests scripts test-data
mypy   phylodater            # blocking gate; must stay at 0 errors
```

**CI gate contract** (`.github/workflows/ci.yml`):

| tool | blocking? | state today |
|------|-----------|-------------|
| `black --check`, `isort --check` | yes | clean; configuration lives in `pyproject.toml` (`line-length = 88`, isort `profile = "black"`) |
| `ruff check` | yes | clean for `phylodater`, `tests`, `scripts` and `test-data` |
| `mypy phylodater` | yes | **0 errors** (53 source files) |

The mypy step is a blocking gate again. `disallow_untyped_defs` and
`warn_return_any` are switched **on** in `pyproject.toml`; the debt they used to
expose (397 errors, 759 with those two flags) has been paid off. Practical rules
for keeping it at zero:

- Annotate, don't suppress. `# type: ignore`, `-> Any` and `cast()` to make an
  error disappear are all off-limits; `cast` is acceptable only for a third-party
  library with no annotations, with a comment naming the real type's source.
- Adapter configuration is typed via the generic base —
  `class MyMethod(DatingMethod[MyConfig])` — and the constructor normalises the
  `ToolConfig | MyConfig` entry before calling `super().__init__`. Never widen
  `self.config` to `Any`.
- Shared value types live in `phylodater/core/method_interface.py`
  (`ConfigT`, `InputFileInfo`); taxonomy maps are `TaxonomyMap = Dict[str, Optional[str]]`
  in `services/taxonomy_parser.py` because absent ranks are deliberately stored as `None`.
- `[tool.mypy] python_version = 3.12` is the *checker* target (numpy's stubs need
  it); the *runtime* floor is `requires-python >= 3.9`, so new code must stay
  3.9-syntax compatible (`typing.TypeGuard` and similar 3.10+ constructs break it).

### 5. Documentation

Every user-facing document exists in English (canonical, bare file name) and
Chinese (`_CN` suffix), and its language switch line lists English first — see the
policy in [docs/README.md](docs/README.md). A change to documented behaviour is not
complete until both files are updated in the same commit.

## Project Structure

```
phylodater/
├── phylodater/              # Main Package
│   ├── adapters/            # Dating method adapters
│   │   ├── lsd2_method.py
│   │   ├── mcmctree_method.py
│   │   ├── mdcat_method.py
│   │   ├── pathd8_method.py
│   │   ├── r8s_pyr8s_method.py
│   │   ├── treepl_method.py
│   │   └── wlogdate_method.py
│   ├── cli/                 # CLI
│   │   └── main.py
│   ├── core/                # Core interfaces and exceptions
│   │   ├── comparison_reporter.py
│   │   ├── exceptions.py
│   │   ├── method_interface.py
│   │   └── pipeline.py
│   ├── infrastructure/      # Infrastructure
│   │   ├── alignment_metadata.py
│   │   ├── auto_naming.py
│   │   ├── checkpoint.py
│   │   ├── configuration.py
│   │   ├── logging.py
│   │   └── process_runner.py
│   ├── models/              # Data models
│   │   ├── calibration.py
│   │   ├── constraints.py
│   │   ├── results.py
│   │   └── tree.py
│   └── services/            # Service layer
│       ├── calibration_loader.py
│       ├── calibration_resolver.py
│       ├── name_mapping.py
│       ├── taxonomy_parser.py
│       ├── tree_rooting.py
│       └── tree_validator.py
├── docs/                   # Documentation
│   ├── api/                # API documentation
│   └── dev/                # Development documentation
├── examples/                # Example data
└── tests/                   # Test files
```

## Adding a New Dating Method

1. Create a new adapter class in `phylodater/adapters/`
2. Inherit the `DatingMethod` base class
3. Implement all abstract methods
4. Register at the end of the adapter file: `DatingMethodRegistry.register("method_name", YourClass)`
5. Add unit tests

### Adapter Template

```python
from pathlib import Path
from typing import List, Optional, Dict, Union

from ..core import DatingMethod, DatingMethodRegistry
from ..models import PhylogeneticTree, CalibrationPoint, DatingResult
from ..infrastructure import ProcessRunner, get_logger
from ..infrastructure.configuration import (
    CommonConfig,
    SoftwarePaths,
    ToolConfig,
    YourMethodConfig,          # 本适配器真正读取的那一份子配置
)


class YourMethod(DatingMethod[YourMethodConfig]):
    """基类按 ConfigT 参数化：self.config 静态视角下就是 YourMethodConfig。"""

    def __init__(
        self,
        config: Union[ToolConfig, YourMethodConfig],
        output_dir: Path,
        software_paths: Optional[SoftwarePaths] = None,
        common_config: Optional[CommonConfig] = None,
    ) -> None:
        # 入口保留双形态（传整个 ToolConfig 或只传本方法子配置都可以），但在
        # 交给基类之前先归一，所以后面 self.config.<参数> 都是类型精确的。
        common: Optional[CommonConfig] = common_config
        if isinstance(config, ToolConfig):
            method_config = config.your_method
            if common is None:
                common = config.common
        else:
            method_config = config
        super().__init__(
            method_config, output_dir, software_paths, common_config=common
        )
        self.logger = get_logger()

    @property
    def method_name(self) -> str:
        return "your_method"

    def validate_environment(self) -> bool:
        runner = ProcessRunner()
        your_bin = self.get_software_path('your_bin') or 'your_binary'
        available = runner.check_executable(your_bin)
        if available:
            version = runner.get_version(your_bin, '--version')
            self.logger.info(f"Your method detected: {version}")
        else:
            self.logger.warning(f"Your method not found: {your_bin}")
        return available

    def prepare_inputs(self,
                       tree: PhylogeneticTree,
                       calibrations: List[CalibrationPoint],
                       alignment_path: Optional[Path] = None) -> Dict:
        pass

    def execute(self) -> bool:
        pass

    def parse_results(self) -> DatingResult:
        pass


DatingMethodRegistry.register("your_method", YourMethod)
```

## Documentation

- User documentation updates [README.md](README.md) and [QUICKSTART.md](QUICKSTART.md)
- API documentation is in `docs/api/`

## Core Module Descriptions

### phylodater.core

The core layer contains the exception hierarchy and the `DatingMethod` interface:

- `method_interface.py`: defines the `DatingMethod` abstract base class and the `DatingMethodRegistry` registry
- `exceptions.py`: all custom exception classes
- `pipeline.py`: `DatingPipeline` and `ParallelPipeline` workflow engines
- `comparison_reporter.py`: result comparison report generation

### phylodater.models

The data model layer contains:

- `constraints.py`: age constraint types (FixedAgeConstraint, UniformAgeConstraint, etc.)
- `calibration.py`: calibration-point models (CalibrationPoint, FossilMetadata, etc.)
- `results.py`: dating result models (DatingResult, NodeAgeEstimate)
- `tree.py`: phylogenetic tree model (PhylogeneticTree)

### phylodater.infrastructure

The infrastructure layer provides:

- `configuration.py`: three-tier configuration management (defaults < YAML < CLI)
- `logging.py`: logging system
- `process_runner.py`: process execution and version detection
- `alignment_metadata.py`: alignment file metadata extraction
- `checkpoint.py`: checkpoint management (MCMCTree support)
- `auto_naming.py`: automatic output file naming

### phylodater.services

The service layer contains:

- `calibration_loader.py`: calibration-point configuration loading
- `calibration_resolver.py`: calibration-point resolution (MRCA localization)
- `name_mapping.py`: taxonomic name mapping management
- `taxonomy_parser.py`: taxonomic data parsing
- `tree_validator.py`: tree structure validation
- `tree_rooting.py`: tree rooting strategies

## License

By submitting code, you agree that your contributions will be licensed under the MIT License.
