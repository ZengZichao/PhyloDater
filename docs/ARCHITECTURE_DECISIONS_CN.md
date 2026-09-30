> **Language:** [English (canonical)](ARCHITECTURE_DECISIONS.md) · 中文（本文档）

# PhyloDater 架构决策记录（ADR）

> 本文档记录 PhyloDater 项目中的关键架构决策，包括决策背景、方案选择和后果权衡。

---

## ADR-001：采用适配器模式集成多个定年引擎

**状态**：已通过（2024-01-15）

**背景**：
用户需要使用不同定年方法（如 MCMCTree、treePL、PATHd8 等）进行比较分析，但各方法的输入格式、输出结构和参数配置差异很大。如果在主管线代码中直接处理这些差异，会导致：
- 代码复杂度随方法数量线性增长
- 新增方法需要修改核心管线
- 各方法之间的实现细节泄漏到业务逻辑层

**决策**：
引入 `DatingMethod` 抽象基类和 `DatingMethodRegistry` 注册表工厂模式。

```python
class DatingMethod(ABC):
    @abstractmethod
    def validate_environment(self) -> bool: pass

    @abstractmethod
    def prepare_inputs(self, tree, calibrations, alignment_path) -> Dict: pass

    @abstractmethod
    def execute(self) -> bool: pass

    @abstractmethod
    def parse_results(self) -> DatingResult: pass
```

**方案选择**：
- **方案 A（本文采用）**：抽象基类 + 注册表工厂
  - 优点：类型安全、编译时验证、自动发现
  - 缺点：需要继承，不够灵活

- **方案 B**：协议（Protocol）+ 依赖注入
  - 优点：无需继承、更灵活
  - 缺点：类型检查更复杂

- **方案 C**：配置驱动 + 反射
  - 优点：零代码新增方法
  - 缺点：运行时错误、性能开销

**后果**：
- **正面**：新增定年方法只需实现接口并调用 `DatingMethodRegistry.register()`，无需修改管线代码
- **负面**：适配器可能泄露实现细节到可视化层（已通过重构 `ComparisonReporter` 解决）
- **负面**：抽象基类的修改会影响所有适配器

**维护策略**：
- `DatingMethod` 接口一旦发布，只增不减
- 新增可选方法使用默认实现
- 破坏性变更通过版本升级处理

---

## ADR-002：标准化 DatingResult 数据结构

**状态**：已通过（2024-02-20）

**背景**：
不同定年引擎输出格式各异：
- MCMCTree: FigTree.tre 格式，带 HPD 区间
- treePL: Newick 格式，仅点估计
- LSD2：支持置信区间

可视化模块需要统一的方式访问结果。

**决策**：
定义 `DatingResult` 和 `NodeAgeEstimate` 标准化数据结构：

```python
@dataclass
class NodeAgeEstimate:
    mean_age: float
    median_age: Optional[float] = None
    ci_lower: Optional[float] = None
    ci_upper: Optional[float] = None
    ci_type: CIType = CIType.NONE

@dataclass
class DatingResult:
    method_name: str
    run_id: str
    dated_tree_newick: str
    node_ages: Dict[str, NodeAgeEstimate]
    is_converged: Optional[bool] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
```

**方案选择**：
- **方案 A（本文采用）**：强类型 dataclass
  - 优点：类型安全、自文档化、支持序列化
  - 缺点：扩展需要修改类定义

- **方案 B**：泛型字典 + JSON Schema
  - 优点：灵活、JSON 原生
  - 缺点：缺乏类型提示、运行时错误

**约束**：
- `metadata` 字段用于存储引擎特定信息，但可视化模块**不应**依赖其内部结构
- 各适配器通过 `ADAPTER_CAPABILITIES` 声明其支持的字段

**后果**：
- **正面**：可视化模块完全解耦于引擎实现
- **正面**：支持结果的序列化与反序列化
- **负面**：部分引擎特有信息（如 treePL 的 optimal_smooth）需要通过扩展字段或能力声明机制暴露

---

## ADR-003：分层配置管理（默认值 < YAML < CLI）

**状态**：已通过（2024-03-01）

**背景**：
用户在不同场景下需要配置参数：
- 个人研究：固定使用某套参数
- 批量分析：共享配置文件
- 特殊需求：临时覆盖特定参数

**决策**：
三层配置优先级：CLI > YAML > 默认值

```python
@dataclass
class PipelineConfig:
    methods: List[str] = field(default_factory=lambda: ['mcmctree', 'treepl', 'pathd8'])
    output_dir: Path = Path("./phylodater_output")
    threads: int = field(default_factory=lambda: get_env_threads(4))
```

环境变量用于运行时覆盖：
- `PHYLODATER_SEED`：全局随机种子
- `PHYLODATER_THREADS`：并行线程数
- `PHYLODATER_BURNIN`: MCMC burnin 代数
- `PHYLODATER_TIMEOUT`：执行超时时间

**后果**：
- **正面**：灵活的配置方式适应不同使用场景
- **正面**：配置文件可版本控制
- **正面**：环境变量支持容器化和自动化

---

## ADR-004：断点续跑机制设计

**状态**：已通过（2024-04-01）

**背景**：
MCMCTree 等贝叶斯方法运行时间可能长达数天至数周。在长时间运行过程中：
- HPC 作业可能因超时而终止
- 系统可能因资源不足崩溃
- 用户可能需要中断后修改参数重新运行

**决策**：
引入两层检查点：

1. **适配器级检查点**（`MCMCTreeCheckpointManager`）：
   - 记录 MCMCTree 内部步骤（Hessian 计算、MCMC 采样等）
   - 支持从上一个完成的步骤恢复

2. **管线级检查点**（`PipelineCheckpointManager`）：
   - 记录各方法的执行状态
   - 缓存 `DatingResult` 支持跨方法恢复
   - 结果持久化到 JSON 文件

```python
class PipelineCheckpointManager:
    def get_methods_to_run(self, requested_methods: List[str]) -> List[str]:
        # 跳过已完成的、保留失败待重试的

    def get_cached_result(self, method_name: str) -> Optional[DatingResult]:
        # 从检查点恢复已完成的结果
```

**方案选择**：
- **方案 A（本文采用）**：结果缓存 + 状态追踪
  - 优点：支持完整恢复、结果可复用
  - 缺点：需要持久化存储

- **方案 B**：仅状态追踪
  - 优点：简单
  - 缺点：恢复时需要重新计算

**后果**：
- **正面**：长时间运行任务的可靠性大幅提升
- **正面**：支持增量分析（新增方法无需重新运行已有方法）
- **负面**：检查点文件需要额外存储空间
- **负面**：结果版本管理复杂度增加

---

## ADR-005：插件系统架构

**状态**：已通过（2024-04-08）

**背景**：
用户可能需要插入自定义的质控或处理步骤，如：
- 长枝吸引序列检测
- 异常序列过滤
- 自定义校准约束

**决策**：
定义插件接口并集成到管线：

```python
class QCPlugin(Plugin):
    def process_tree(self, tree: Any, context: Dict[str, Any]) -> Any: pass

class PluginRegistry:
    _plugins: Dict[str, LoadedPlugin]

class PluginLoader:
    def discover_plugins(self, search_paths: List[Path]) -> List[Path]: pass
```

**方案选择**：
- **方案 A（本文采用）**：抽象基类 + 文件系统发现
  - 优点：类型安全、自动发现
  - 缺点：需要继承

- **方案 B**：装饰器 + 入口点
  - 优点：无需继承
  - 缺点：隐式依赖

**后果**：
- **正面**：用户可定制质控流程
- **正面**：生态扩展性
- **负面**：插件质量难以控制
- **负面**：引入安全考虑（待加载代码执行）

---

## ADR-006：错误分类与可操作性建议

**状态**：已通过（2024-04-09）

**背景**：
生物信息学管线运行中会遇到多种错误：
- 用户输入错误（文件格式、参数错误）
- 引擎运行失败（收敛失败、内存溢出）
- 系统资源不足（超时、内存不足）

用户需要可操作的错误信息来快速定位问题。

**决策**：
引入 `ErrorCategory` 枚举和结构化错误建议：

```python
class ErrorCategory(Enum):
    USER_INPUT = "user_input"
    ENGINE_FAILURE = "engine_failure"
    RESOURCE_LIMITED = "resource_limited"
    SYSTEM_ERROR = "system_error"

class ExecutionError(PhyloDaterError):
    error_category: ErrorCategory
    recoverable: bool
    suggestion: Optional[str]
```

错误分类辅助函数根据外部程序输出版本分类：

```python
def classify_execution_error(stderr: str, returncode: int) -> tuple:
    if "out of memory" in stderr_lower:
        return (ErrorCategory.RESOURCE_LIMITED, "请尝试减少数据规模...")
```

**后果**：
- **正面**：用户可快速判断错误类型和解决方案
- **正面**：支持自动化的错误处理策略
- **正面**：支持结构化日志和监控

---

## ADR-009：方法级校准兼容性预检

**状态**：已通过（2026-07-12）

**背景**：
不同定年方法对校准约束的支持存在差异。例如 MCMCTree 的 `RootAge` 只能设置为上界，无法真正固定根节点年龄；如果用户提供 `fixed` 根节点约束，MCMC 采样会收敛到低于该值的某一点，导致结果偏离预期。此前这类问题要到运行中期才会暴露，浪费计算资源且用户体验差。

**决策**：
在 `DatingMethod` 抽象基类中新增 `validate_calibrations(calibrations)` 钩子，各适配器可在 `prepare_inputs()` 之前执行方法专属的校准检查。

```python
class DatingMethod(ABC):
    def validate_calibrations(self, calibrations: List[CalibrationPoint]) -> List[str]:
        """返回错误信息列表；空列表表示通过。"""
        return []
```

管线在 `validate_environment()` 之后、`prepare_inputs()` 之前调用 `validate_calibrations()` 钩子。若返回非空，则记录 ERROR、标记该方法失败，并继续运行其他方法。

**方案选择**：
- **方案 A（本文采用）**：基类钩子 + 适配器重写
  - 优点：职责清晰，新增方法只需重写钩子即可声明自身限制
  - 缺点：需要每个适配器主动实现

- **方案 B**：全局校准规则表
  - 优点：集中管理
  - 缺点：难以表达复杂逻辑（如 MCMCTree 的 RootAge 语义），新增规则需修改核心代码

- **方案 C**：完全由外部软件自己报错
  - 优点：零额外代码
  - 缺点：错误发现较晚、提示信息不友好、浪费计算资源

**后果**：
- **正面**：PhyloDater 在执行早期发现不兼容配置，避免运行到一半才失败
- **正面**：错误信息可包含针对性的修改建议（如改为 uniform、maximum 或换用其他方法）
- **正面**：单方法失败不影响其他方法继续运行
- **负面**：需要为每种方法维护校验逻辑

**更新记录**：
- 2026-07-12：通过 ADR；已在 MCMCTree 适配器中实现 fixed root 拦截。

---

## 待决策

### ADR-007：可视化模块与引擎解耦（进行中）

**问题**：当前 `DatingFigure` 组件直接解析 FigTree.tre 格式，这违反了 ADR-002 的设计原则。

**候选方案**：
- **方案 A**：只依赖 `DatingResult`，增强 `NodeAgeEstimate` 字段
- **方案 B**：引入“可视化适配器”抽象层
- **方案 C**：扩展 `ADAPTER_CAPABILITIES` 机制

### ADR-008：批量处理架构

**问题**：当前批量处理需要外部脚本循环调用 `DatingPipeline`。

**候选方案**：
- **方案 A**：扩展 `PipelineConfig` 支持多样本配置
- **方案 B**：独立的 `BatchProcessor` 类（已实现）
- **方案 C**：工作流语言集成（Snakemake、Nextflow）

---

## 维护指南

### 如何添加新的 ADR

1. 在本文档末尾添加新章节
2. 包含：状态、背景、决策、方案选择、后果
3. 状态选项：提议中、已接受、已拒绝、已废弃

### 如何更新已有 ADR

1. 在 ADR 末尾添加“更新记录”段落
2. 说明更新原因和影响
3. 保持决策历史完整

### ADR 审查周期

- 重大版本发布前审查所有“待决策”ADR
- 新增外部依赖时评估是否需要新 ADR
- 定期清理已废弃的 ADR
