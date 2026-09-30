> **Language:** [English (canonical)](ARCHITECTURE.md) · 中文（本文档）

# PhyloDater 架构文档

本文档描述 PhyloDater 的软件架构设计。

## 架构概览

PhyloDater 采用分层架构设计，遵循关注点分离原则：

```
┌─────────────────────────────────────────────────────────────┐
│                        CLI Layer                            │
│                    (命令行接口)                              │
├─────────────────────────────────────────────────────────────┤
│                      Service Layer                          │
│    (Calibration, Tree Validation, Name Mapping)             │
├─────────────────────────────────────────────────────────────┤
│                       Core Layer                            │
│    (Method Interface, Pipeline, Exceptions)                 │
├─────────────────────────────────────────────────────────────┤
│                      Adapter Layer                          │
│  (MCMCTree, LSD2, r8s, treePL, PATHd8, wLogDate, MD-Cat)    │
├─────────────────────────────────────────────────────────────┤
│                      Model Layer                            │
│     (Constraints, Calibration, Results, Tree)               │
├─────────────────────────────────────────────────────────────┤
│                   Infrastructure Layer                      │
│  (Configuration, Logging, Process Runner, Checkpoint)       │
└─────────────────────────────────────────────────────────────┘
```

## 各层详细说明

### 1. CLI Layer（命令行接口层）

**职责**：提供用户交互接口

**主要组件**：
- `cli/main.py`：命令行参数解析和主流程控制

**设计原则**：
- 单一职责：只负责用户输入和输出展示
- 依赖下层：通过 Service Layer 和 Core Layer 执行业务逻辑

### 2. Service Layer（服务层）

**职责**：封装业务逻辑，协调领域对象

**主要组件**：
- `calibration_loader.py`：校准点配置加载
- `calibration_resolver.py`：校准点解析（MRCA 定位）
- `tree_validator.py`：树结构验证
- `tree_rooting.py`：树定根策略
- `name_mapping.py`：分类学名称映射管理
- `taxonomy_parser.py`：分类学数据解析
- `input_validator.py`：输入文件格式验证
- `deep_validator.py`：深度验证（括号平衡、负分支、恶意字符等）
- `auto_calibrator.py`：自动校准点添加

**设计原则**：
- 高内聚：每个服务专注于特定业务功能
- 可测试：业务逻辑独立于外部依赖

### 3. Core Layer（核心层）

**职责**：定义核心业务逻辑和接口

**主要组件**：
- `method_interface.py`：定年方法抽象接口和注册表
- `pipeline.py`：工作流引擎（顺序和并行）
- `comparison_reporter.py`：结果比较报告生成
- `exceptions.py`：异常层次结构

**关键设计模式**：
- **策略模式**：`DatingMethod` 接口允许不同定年算法的统一调用
- **注册表模式**：`DatingMethodRegistry` 实现方法动态发现和创建
- **工厂模式**：通过注册表创建具体方法实例
- **模板方法模式**：`DatingMethod.validate_calibrations()` 等可选钩子提供默认实现，允许适配器扩展方法专属校验

### 4. Adapter Layer（适配器层）

**职责**：集成外部定年软件

**主要组件**：
- `mcmctree_method.py`: PAML MCMCTree 适配器
- `lsd2_method.py`: IQ-TREE2 LSD2 适配器
- `r8s_pyr8s_method.py`：r8s 与 pyr8s 适配器
- `treepl_method.py`: treePL 适配器
- `pathd8_method.py`: PATHd8 适配器
- `wlogdate_method.py`: wLogDate 适配器
- `mdcat_method.py`: MD-Cat 适配器

**设计原则**：
- **适配器模式**：将不同软件的接口统一为 `DatingMethod` 接口
- **单一职责**：每个适配器只负责一个外部软件
- **错误隔离**：外部软件的错误不影响其他组件

### 5. Model Layer（模型层）

**职责**：定义领域模型和数据结构

**主要组件**：
- `constraints.py`：年龄约束类型体系
- `calibration.py`：校准点相关模型
- `results.py`：定年结果模型
- `tree.py`：系统发育树模型

**设计原则**：
- **值对象**：约束类型使用 dataclass 实现不可变对象
- **富领域模型**：模型包含业务逻辑（如格式转换）

### 6. Infrastructure Layer（基础设施层）

**职责**：提供技术能力和横切关注点

**主要组件**：
- `configuration.py`：三层配置管理（默认值 < YAML < CLI）
- `logging.py`：日志系统
- `process_runner.py`：进程执行和版本检测
- `checkpoint.py`：检查点管理（MCMCTree 支持）
- `auto_naming.py`：输出文件自动命名
- `alignment_metadata.py`：比对文件元数据提取

**设计原则**：
- **依赖倒置**：上层依赖抽象接口，不依赖具体实现
- **横切关注点**：日志、配置等功能集中管理

## 数据流

```
输入文件
    │
    ▼
┌─────────────┐
│ CLI/Service │
│   Layer     │
└──────┬──────┘
       │
       ▼
┌─────────────┐     ┌─────────────┐
│   Model     │────▶│  Adapter    │
│   Layer     │     │   Layer     │
└─────────────┘     └──────┬──────┘
                           │
                           ▼
                    ┌─────────────┐
                    │  External   │
                    │  Software   │
                    └──────┬──────┘
                           │
                           ▼
                    ┌─────────────┐
                    │   Results   │
                    │   Layer     │
                    └──────┬──────┘
                           │
                           ▼
                    ┌─────────────┐
                    │  Comparison │
                    │  Reporter   │
                    └──────┬──────┘
                           │
                           ▼
                       输出文件
```

## 关键设计决策

### 1. 约束类型体系

**决策**：使用 dataclass 实现年龄约束类型层次结构

**理由**：
- 类型安全：编译时检查
- 不可变性：避免意外修改
- 自文档化：代码即文档

### 2. 适配器注册机制

**决策**：使用类级别的注册表，而非配置文件

**理由**：
- 类型安全：编译时验证适配器存在
- 自动发现：新适配器自动注册
- 零配置：无需维护配置文件

### 3. 三层配置管理

**决策**：内置默认值 < YAML 配置 < 命令行参数

**理由**：
- 灵活性：支持多种配置方式
- 可预测性：优先级明确
- 可重复性：配置文件可版本控制

### 4. 异常层次结构

**决策**：定义清晰的异常类型层次

**理由**：
- 精确错误处理：调用者可以捕获特定异常
- 错误分类：验证错误、执行错误、结果错误
- 语义降级警告：非致命问题使用警告而非异常

### 5. 适配器配置类型

**决策**：抽象基类按"它实际读取的那份配置"参数化 —— `DatingMethod[ConfigT]`，
每个适配器把该参数钉到自己的配置 dataclass 上
（`TreePLMethod(DatingMethod[TreePLConfig])`、`MCMCTreeMethod(DatingMethod[MCMCTreeConfig])` …）。

**它解决的问题**：适配器历史上同时接受**整个** `ToolConfig` **或**单个方法的子配置
（如 `TreePLConfig`），用 `hasattr(config, "treepl")` 在运行期区分两者，再把子配置
赋给 `self.config`。基类却把 `config` 声明为 `ToolConfig`，于是每一个
`self.config.treepl_bin` / `self.config.numsites` 在静态视角下都是"`ToolConfig`
没有这个属性"—— 223 条 `attr-defined`，而且无法把真正的拼写错误从这层鸭子类型里分出来。

**为什么保留双形态入口**：测试套件里两种形态都被大量使用
（`MCMCTreeMethod(MCMCTreeConfig(), …)` 与 `WLogDateMethod(ToolConfig(), …)`），
若把入口收敛成"只收 `ToolConfig`"就得重写这些调用点及其断言。所以改成在门口归一：

```python
if isinstance(config, ToolConfig):      # 与原先的 hasattr 判定等价：
    method_config = config.treepl       # ToolConfig 必有 .treepl，
                                        # TreePLConfig 必无
else:
    method_config = config
super().__init__(method_config, output_dir, software_paths, common_config=common)
```

这样对任何读者来说 `self.config` 就是精确的 `TreePLConfig`，全程没有用
`self.config: Any` 逃生。`DatingMethod.__init__` 同时接收 `common_config` 并归一成
非 `Optional` 的 `CommonConfig`，这是那十五处 `self.common_config.nthreads / .timeout
/ .seed` 读取所需要的。

**相关的共享类型**：`InputFileInfo = Union[Path, int, float]` —— `prepare_inputs()`
返回的"输入文件字典"确实同时装文件路径与每次运行的统计量（位点数、耗时、对数似然）。
`DatingMethod.cleanup()` 会跳过非路径条目；在这件事被写进类型之前，它对它们执行
`Path(value)`，而 `Path(1000)` 抛的是 `TypeError`，并不在 `cleanup()` 的 `except` 里。

**被否决的替代方案**：`self.config: Any = config`。它同样能消掉那 223 条报错，
但会把适配器配置访问上的全部检查删掉 —— 而那恰恰是"引擎参数拼错就静默失效"的现场。

## 扩展指南

### 添加新的定年方法

1. 在 `phylodater/adapters/` 创建新的适配器类
2. 继承 `DatingMethod` 抽象基类
3. 实现所有抽象方法
4. 在文件末尾注册：`DatingMethodRegistry.register("method_name", YourClass)`
5. 添加单元测试

### 添加新的约束类型

1. 在 `phylodater/models/constraints.py` 创建新的 dataclass
2. 继承 `AgeConstraint` 抽象基类
3. 实现 `to_software_format()` 和 `to_mcmctree_calib_string()`
4. 添加验证逻辑到 `__post_init__()`
5. 添加单元测试

## 性能考虑

- **并行执行**：`ParallelPipeline` 支持多方法并行运行
- **检查点**：MCMCTree 支持检查点恢复，避免重复计算
- **内存管理**：大型数据集使用流式处理
- **缓存**：配置和元数据缓存避免重复计算

## 安全考虑

- **输入验证**：所有外部输入都经过验证
- **命令注入防护**：使用列表而非字符串传递命令参数
- **资源限制**：支持超时和内存限制

## 测试策略

- **单元测试**：每个模块独立测试
- **集成测试**：测试模块间交互
- **适配器测试**：使用 mock 测试外部软件集成
- **端到端测试**：完整工作流测试
