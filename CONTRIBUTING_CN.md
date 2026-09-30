> **Language:** [English (canonical)](CONTRIBUTING.md) · 中文（本文档）

# 贡献指南

感谢你考虑为 PhyloDater 贡献代码。

## 如何贡献

### 报告问题

如果你发现 bug 或有新功能请求，请：

1. 搜索现有 issues，确认该问题还没有人报告
2. 创建新的 issue，包含：
   - 清晰的标题和描述
   - 复现步骤
   - 预期行为与实际行为
   - 环境信息（Python 版本、操作系统等）

### 提交代码

1. Fork 本仓库
2. 创建特性分支：`git checkout -b feature/your-feature-name`
3. 完成开发并提交更改
4. 确保测试通过：`pytest tests/`
5. 提交 Pull Request

### 代码规范

- 遵循 PEP 8
- 使用类型注解
- 添加 docstrings
- 编写单元测试

## 开发环境设置

### 1. 克隆仓库

```bash
git clone https://github.com/yourusername/phylodater.git
cd phylodater
```

### 2. 创建开发环境

```bash
conda create -n phylodater-dev python=3.12
conda activate phylodater-dev

pip install -e ".[dev]"

conda install -c bioconda -c conda-forge paml iqtree   # MCMCTree + LSD2
pip install wlogdate
pip install git+https://github.com/iTaxoTools/pyr8s.git   # r8s backend (NPRS only)
pip install git+https://github.com/uym2/MD-Cat.git        # MD-Cat
# treePL 与 PATHd8：源码构建（见 docs/EXTERNAL_TOOLS_CN.md §3-§4）；
# 两者均已不在 bioconda/conda-forge/PyPI 上发布。
```

### 3. 运行测试

```bash
python -m pytest tests -q                 # 开发者套件 + 发布验证套件
python -m pytest tests -m "not slow" -q   # 跳过调用真实引擎的测试
python -m pytest tests/validation -q      # 仅发布验证套件
phylodater check                          # 26 项环境自检
```

依赖外部引擎的测试在二进制缺失时会带原因自行跳过，因此即使机器上七个引擎一个也没装，
套件仍然可以运行。

新增定年方法时必须同时给它配上真值运行：`python scripts/run_reference_benchmark.py`
会拿 `test-data/benchmark/` 的仿真真值为每个方法打分——验收门槛见
[docs/TESTING_CN.md](docs/TESTING_CN.md)，当前数值见
[docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md)。

### 4. 代码格式化、Lint 与类型检查

```bash
black phylodater tests scripts test-data
isort  phylodater tests scripts test-data
ruff check phylodater tests scripts test-data
mypy   phylodater            # 阻断门禁，必须保持 0 error
```

**CI 关卡约定**（`.github/workflows/ci.yml`）：

| 工具 | 是否阻断 | 当前状态 |
|------|----------|----------|
| `black --check`、`isort --check` | 是 | 已全部通过；配置在 `pyproject.toml`（`line-length = 88`，isort `profile = "black"`） |
| `ruff check` | 是 | `phylodater`、`tests`、`scripts`、`test-data` 已全部通过 |
| `mypy phylodater` | 是 | **0 个错误**（53 个源文件） |

mypy 已重新成为阻断门禁。`pyproject.toml` 里 `disallow_untyped_defs` 与
`warn_return_any` 均为**开启**；它们曾经暴露的欠债（397 个错误；打开两个开关则 759 个）
已清偿。保持为 0 的实操规则：

- 标注，而不是压制。用 `# type: ignore`、把返回值写成 `Any`、或靠 `cast()` 让报错消失
  都在禁止之列；`cast` 只允许用于"第三方库无标注"处，且要用注释写明真实类型的来源。
- 适配器配置通过泛型基类定型——`class MyMethod(DatingMethod[MyConfig])`——构造函数在
  调用 `super().__init__` 之前先把 `ToolConfig | MyConfig` 的入口归一。绝不把
  `self.config` 放宽成 `Any`。
- 共享值类型放在 `phylodater/core/method_interface.py`（`ConfigT`、`InputFileInfo`）；
  分类学字典统一用 `services/taxonomy_parser.py` 里的
  `TaxonomyMap = Dict[str, Optional[str]]`，因为缺失的等级是**故意**存成 `None` 的。
- `[tool.mypy] python_version = 3.12` 是**检查器**目标（numpy 的桩需要它）；
  **运行期**下限是 `requires-python >= 3.9`，所以新代码必须保持 3.9 语法兼容
  （`typing.TypeGuard` 一类 3.10+ 的构造会破坏它）。

### 5. 文档

所有面向用户的文档都有英文（权威、裸文件名）与中文（`_CN` 后缀）两个版本，且顶部的
语言切换行将英文排在首位——完整约定见 [docs/README_CN.md](docs/README_CN.md)。对已
记录行为的修改，必须在同一次提交内同时更新两份文件，否则该修改不算完成。

## 项目结构

```
phylodater/
├── phylodater/              # 主包
│   ├── adapters/            # 定年方法适配器
│   │   ├── lsd2_method.py
│   │   ├── mcmctree_method.py
│   │   ├── mdcat_method.py
│   │   ├── pathd8_method.py
│   │   ├── r8s_pyr8s_method.py
│   │   ├── treepl_method.py
│   │   └── wlogdate_method.py
│   ├── cli/                 # 命令行接口
│   │   └── main.py
│   ├── core/                # 核心接口和异常
│   │   ├── comparison_reporter.py
│   │   ├── exceptions.py
│   │   ├── method_interface.py
│   │   └── pipeline.py
│   ├── infrastructure/      # 基础设施
│   │   ├── alignment_metadata.py
│   │   ├── auto_naming.py
│   │   ├── checkpoint.py
│   │   ├── configuration.py
│   │   ├── logging.py
│   │   └── process_runner.py
│   ├── models/              # 数据模型
│   │   ├── calibration.py
│   │   ├── constraints.py
│   │   ├── results.py
│   │   └── tree.py
│   └── services/            # 服务层
│       ├── calibration_loader.py
│       ├── calibration_resolver.py
│       ├── name_mapping.py
│       ├── taxonomy_parser.py
│       ├── tree_rooting.py
│       └── tree_validator.py
├── docs/                   # 文档
│   ├── api/                # API 文档
│   └── dev/                # 开发文档
├── examples/                # 示例数据
└── tests/                   # 测试文件
```

## 添加新的定年方法

1. 在 `phylodater/adapters/` 创建新的适配器类
2. 继承 `DatingMethod` 基类
3. 实现所有抽象方法
4. 在适配器文件末尾注册：`DatingMethodRegistry.register("method_name", YourClass)`
5. 添加单元测试

### 适配器模板

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

## 文档

- 用户文档更新 [README_CN.md](README_CN.md) 和 [QUICKSTART_CN.md](QUICKSTART_CN.md)
- API 文档在 `docs/api/`

## 核心模块说明

### phylodater.core

核心层包含异常层次结构和 `DatingMethod` 接口：

- `method_interface.py`：定义 `DatingMethod` 抽象基类和 `DatingMethodRegistry` 注册表
- `exceptions.py`：所有自定义异常类
- `pipeline.py`：`DatingPipeline` 和 `ParallelPipeline` 工作流引擎
- `comparison_reporter.py`：结果比较报告生成

### phylodater.models

数据模型层包含：

- `constraints.py`：年龄约束类型（FixedAgeConstraint、UniformAgeConstraint 等）
- `calibration.py`：校准点相关模型（CalibrationPoint、FossilMetadata 等）
- `results.py`：定年结果模型（DatingResult、NodeAgeEstimate）
- `tree.py`：系统发育树模型（PhylogeneticTree）

### phylodater.infrastructure

基础设施层提供：

- `configuration.py`：三层配置管理（默认值 < YAML < CLI）
- `logging.py`：日志系统
- `process_runner.py`：进程执行和版本检测
- `alignment_metadata.py`：比对文件元数据提取
- `checkpoint.py`：检查点管理（MCMCTree 支持）
- `auto_naming.py`：输出文件自动命名

### phylodater.services

服务层包含：

- `calibration_loader.py`：校准点配置加载
- `calibration_resolver.py`：校准点解析（MRCA 定位）
- `name_mapping.py`：分类学名称映射管理
- `taxonomy_parser.py`：分类学数据解析
- `tree_validator.py`：树结构验证
- `tree_rooting.py`：树定根策略

## 许可证

通过提交代码，你同意你的贡献将遵循 MIT 许可证。
