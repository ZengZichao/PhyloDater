> **Language:** [English (canonical)](README.md) · 中文（本文档）

# 发布验证测试套件

本套件按评审者的方式验证 PhyloDater：真实调用各定年引擎端到端运行，并检查它们返回的
数值，而不只是测试外围的 Python 管道。

所有内容都被限制在本目录内：测试不会修改 `tests/validation/` 之外的文件（临时输出写入
`tests/validation/_tmp/`，已被 git 忽略）。

## 目录结构

```
tests/validation/
├── conftest.py     共享 fixture、自定义 marker、临时目录管理
├── helpers.py      引擎可用性探测、CLI 运行器、fixture 构造器
├── unit/           校准模型/加载/解析、树模型与拓扑、检查点、配置、日志、方法注册
├── functional/     CLI 子命令（dating / check / fold）、dry-run、错误处理、自动校准
├── adapters/       适配器契约、MCMCTree 适配器、其他适配器的环境检测
├── integration/    每个已安装引擎的端到端运行、多方法调度、检查点恢复
├── security/       输出隔离、路径穿越、YAML 与输入安全
└── fixtures/       静态输入：树、比对、校准、配置、分类学、基准
```

## 与 `tests/` 其余部分的关系

| 目录 | 用途 | 是否需要外部软件 |
|------|------|------------------|
| `tests/unit`、`tests/functional`、`tests/integration` | 开发者套件：快速、大多使用 mock | 仅"引擎原始调用"冒烟测试需要 |
| `tests/validation` | 发布验证：真实引擎、端到端 | 需要；引擎缺失时带原因跳过 |
| `test-data/benchmark` + `scripts/run_reference_benchmark.py` | 对仿真真值的**精度**评分 | 需要 |

本套件的实测结论记录在 [docs/TEST_RESULTS_CN.md](../../docs/TEST_RESULTS_CN.md)——
本目录内刻意不再保留第二份容易过期的结果文件。

## 运行方式

```bash
python -m pytest tests/validation -v          # 全部（按环境可用的引擎执行）
python -m pytest tests/validation/unit -v     # 单个类别
python -m pytest tests/validation -m slow     # 只跑真实引擎的慢速测试
python -m pytest tests/validation -m "not slow"
python -m pytest tests/validation -m security
PHYLODATER_KEEP_TMP=1 python -m pytest tests/validation/integration -v   # 保留临时目录
```

`python -m pytest tests`（`pyproject.toml` 中默认的 `testpaths`）已经包含本套件。

## 自定义 Marker

| marker | 含义 |
|--------|------|
| `slow` | 会真正调用外部引擎（数秒到数分钟） |
| `security` | 路径安全、注入、输入校验相关测试 |
| `requires_mcmctree` … `requires_mdcat` | 每个定年引擎一个；引擎缺失时带明确原因跳过 |

可用性由 `helpers.external_available()` 探测，它知道某些方法的后端是 Python 包而非
可执行文件：

| 方法 | 必须存在的东西 |
|------|----------------|
| `mcmctree` | `mcmctree`（以及 `baseml`） |
| `lsd2` | `iqtree2` 或 `iqtree`（LSD2 内建于 IQ-TREE 2） |
| `r8s` | `pyr8s` 包**或** `r8s` 二进制 |
| `wlogdate` | `logdate` 包**或** `launch_wLogDate.py` |
| `mdcat` | `emd` 包**或** `md_cat.py` / `mdcat` |
| `pathd8` | `PATHd8` 二进制 |
| `treepl` | `treePL` 二进制 |

引擎安装与溯源：[docs/EXTERNAL_TOOLS_CN.md](../../docs/EXTERNAL_TOOLS_CN.md)。
测试方法与验收门槛：[docs/TESTING_CN.md](../../docs/TESTING_CN.md)。
