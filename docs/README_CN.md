> **Language:** [English (canonical)](README.md) · 中文（本文档）

# 文档

## 文档语言政策

本仓库中每份文档都有**两个**版本：

| 规则 | 含义 |
|------|------|
| 英文为权威版本 | 英文文本使用裸文件名（`README.md`、`docs/TESTING.md`），是项目视为权威的版本 |
| 中文为对照译本 | 中文文本使用 `_CN` 后缀（`README_CN.md`、`docs/TESTING_CN.md`） |
| 显示时英文在前 | 每份文档顶部的语言切换行，无论文件本身是什么语言，都把英文排在中文之前 |
| 两版同步修改 | 涉及已记录行为的改动，必须在同一次提交内同时改英文文件与其 `_CN` 对照文件 |
| 机器可读元数据 | `CITATION.cff` 受格式限制只有英文，其配套的中文说明是 `CITATION_CN.md` |
| Issue 模板 | `.github/ISSUE_TEMPLATE/*.md` 保持每个入口一个文件：英文小标题下附中文说明。GitHub 按文件列出模板，再拆一份 `_CN` 只会让选择列表重复 |

## 用户文档

| 文档 | 内容 |
|------|------|
| [../README_CN.md](../README_CN.md) | 功能总览、安装、快速示例、方法与参数参考 |
| [../QUICKSTART_CN.md](../QUICKSTART_CN.md) | 五分钟上手、常见工作流、故障排查 |
| [EXTERNAL_TOOLS_CN.md](EXTERNAL_TOOLS_CN.md) | 各定年引擎的获取方式、版本溯源、源码构建、Docker 覆盖情况 |
| [TESTING_CN.md](TESTING_CN.md) | 如何测试：分层、带真值基准的设计、评分规则、验收门槛、复现命令 |
| [TEST_RESULTS_CN.md](TEST_RESULTS_CN.md) | 七个引擎的实测结果、对仿真真值的精度、确定性、发现并修复的缺陷 |
| [ARCHITECTURE_CN.md](ARCHITECTURE_CN.md) | 包结构与模块职责 |
| [ARCHITECTURE_DECISIONS_CN.md](ARCHITECTURE_DECISIONS_CN.md) | ADR 记录：决定了什么、为什么、拒绝了什么 |
| [ERROR_CODES_CN.md](ERROR_CODES_CN.md) | 退出码与异常类含义 |
| [CONDA_PACKAGING_CN.md](CONDA_PACKAGING_CN.md) | 为 conda 打包 PhyloDater（及其引擎） |
| [../CONTRIBUTING_CN.md](../CONTRIBUTING_CN.md) | 开发环境、代码风格、测试与文档要求 |
| [../SECURITY_CN.md](../SECURITY_CN.md) | 漏洞报告与范围 |
| [../THIRD_PARTY_LICENSES_CN.md](../THIRD_PARTY_LICENSES_CN.md) | 依赖与引擎许可 |

## 随包发布的测试数据

| 路径 | 内容 |
|------|------|
| [`../test-data/benchmark/`](../test-data/benchmark/) | 带真值的精度基准（仿真的树/比对、真实节点年龄、校准、生成脚本、归档结果） |
| [`../test-data/input/`](../test-data/input/) | 合法 / 非法 / 边界场景输入，用于手工端到端核查 |
| [`../examples/`](../examples/) | 快速入门使用的小型示例输入 |
| [`../tests/`](../tests/) | 自动化测试：`unit`、`functional`、`integration`，以及 `validation` 发布验证套件 |

## API 参考（Sphinx）

API 参考由 Sphinx + numpydoc 从 docstring 生成。

```bash
pip install -e ".[dev]"
sphinx-build -b html docs/api docs/api/_build
```

来源文件：`api/index.rst`（toctree）、`api/conf.py`（配置），以及每个包模块对应的
一个 `.rst`。新增模块时需要添加其 `phylodater.<module>.rst` 并在 `api/index.rst` 中列入。

## 新增或修改文档

1. 先写英文文件。
2. 在同一次提交内，把相同的结构性改动应用到 `_CN` 文件。
3. 保留文件顶部的语言切换行（英文在前）。
4. 把新文档从上面的表格以及根 README 链接出去，保证两个入口都能找到。

上述约定由机器校验：`python scripts/check_docs.py` 会在内部链接失效、某份文档失去
对照版本、`*_EN.md` 重复文件重新出现、或语言切换行不再"英文在前"时失败。CI 已把它作为
一个 lint 步骤运行。
