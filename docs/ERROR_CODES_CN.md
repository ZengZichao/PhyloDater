> **Language:** [English (canonical)](ERROR_CODES.md) · 中文（本文档）

# PhyloDater 错误码

| 错误码 | 名称 | 说明 |
|------|------|-------------|
| 0 | EXIT_SUCCESS | 执行成功完成 |
| 1 | EXIT_RUNTIME_ERROR | 运行时错误（缺少依赖、自检失败） |
| 2 | EXIT_USAGE_ERROR | 命令行用法错误（参数无效、缺必需选项） |
| 3 | EXIT_DATA_ERROR | 输入数据错误（格式无效、文件未找到、负分支长度、不支持 PHYLIP 格式、空校准、标签不匹配） |
| 130 | EXIT_INTERRUPTED | 用户中断（Ctrl+C） |

## 错误码 3 的常见场景

- **File not found**：`--tree`、`--sequence`、`--calibrations` 的路径不存在
- **Negative branch lengths**：树的某些分支长度为负值（可使用 `--skip-length-check` 跳过）
- **PHYLIP format**：不支持 PHYLIP 比对格式，请转为 FASTA（可使用 `--ignore-malformed` 以警告方式继续）
- **Empty calibration**：校准文件的有效约束为空
- **Label mismatch**：树与序列文件的末端标签不匹配
- **Input too large**：文件超过 500 MB 大小限制
- **Invalid tree structure**：括号不平衡、重复名称
- **Unsupported format**：PhyloDater 无法识别系统树或比对格式
- **Circular taxonomy dependencies**：分类表包含循环引用

## 方法专属校准错误

部分定年方法对校准约束有额外限制，PhyloDater 会在执行前调用各适配器的 `validate_calibrations()` 钩子检查，避免运行到一半才失败。

### MCMCTree 不支持固定根节点

**错误信息示例**：
```
方法 'mcmctree' 的校准配置不兼容：
  - MCMCTree 不支持固定根节点年龄。根节点校准 'Root' 当前为 fixed=500.0 Ma，
    但 MCMCTree 只能将其作为上界（RootAge='<500.0'）使用，导致根节点年龄无法精确固定。
  - 建议修改方案（任选其一）：
    1) 将根节点约束改为 uniform，例如 min=499.9, max=500.1；
    2) 将根节点约束改为 maximum，例如 max=500.0；
    3) 如果必须精确固定根节点，请改用 treePL/pathd8/lsd2 等方法。
```

**原因**：MCMCTree 的 `RootAge` 控制文件参数只能设置上界（如 `<500`），无法真正固定根节点年龄。若用户提供 `fixed` 约束，MCMC 采样后根节点年龄会收敛到低于该值的某一点，导致结果偏离预期。

**解决方案**：
1. 将根节点约束改为 `uniform`，例如 `min=499.9, max=500.1`；
2. 将根节点约束改为 `maximum`，例如 `max=500.0`；
3. 如需严格固定根节点，改用 treePL、PATHd8 或 LSD2 等方法。

> **提示**：不同定年方法可使用不同的校准配置文件。例如 MCMCTree 使用 `maximum` 根节点约束，而其他方法继续使用 `fixed` 根节点约束。

## 异常类

PhyloDater 使用结构化的异常类实现精确定位的错误处理：

| 异常类 | 父类 | 退出码 | 说明 |
|----------------|--------|-----------|-------------|
| `PhyloDaterError` | Exception | - | 所有 PhyloDater 错误的基类 |
| `ValidationError` | PhyloDaterError | 3 | 输入验证错误的基类 |
| `TreeValidationError` | ValidationError | 3 | 树结构校验失败 |
| `AlignmentValidationError` | ValidationError | 3 | 序列比对校验失败 |
| `ConsistencyError` | ValidationError | 3 | 树-序列标签不匹配 |
| `PhyloFormatError` | ValidationError | 3 | 不支持的文件格式或格式错误 |
| `NegativeBranchLengthError` | TreeValidationError | 3 | 树包含负分支长度 |
| `InputSizeLimitError` | ValidationError | 3 | 输入文件超出大小限制（500 MB） |
| `CalibrationError` | PhyloDaterError | 3 | 校准配置错误 |
| `CalibrationConflictError` | CalibrationError | 3 | 校准约束冲突 |
| `CalibrationResolutionError` | CalibrationError | 3 | 校准点解析失败 |
| `TaxonomyConflictError` | CalibrationError | 3 | 分类表存在循环引用 |
| `MonophylyError` | CalibrationError | 3 | 分类单元不构成单系群 |
| `ConfigurationError` | PhyloDaterError | 2 | 配置文件或参数错误 |
| `UnknownMethodError` | PhyloDaterError | 2 | 请求的定年方法未注册 |
| `ExecutionError` | PhyloDaterError | 1 | 外部软件执行失败 |
| `ExecutionTimeoutError` | ExecutionError | 1 | 外部软件超时 |
| `ResultParsingError` | PhyloDaterError | 1 | 解析外部软件输出失败 |

### 告警类

非致命问题使用告警类（仅记录，不终止运行）：

| 告警类 | 父类 | 说明 |
|---------------|--------|-------------|
| `ConvergenceWarning` | PhyloDaterWarning | MCMC 收敛诊断阈值超限 |
| `SemanticDegradationWarning` | PhyloDaterWarning | 输入语义降级（如约束类型降级） |
| `EnvironmentWarning` | PhyloDaterWarning | 环境问题（如可选依赖缺失） |
| `ResourceWarning` | PhyloDaterWarning | 资源使用告警（如检测到大型树） |
| `UnknownParameterWarning` | PhyloDaterWarning | 遇到未知的 method-args 参数 |

## 分类表格式

PhyloDater 支持两种分类学注释格式：

### 格式 A（嵌入在 tip 名称中）
```
GB_GCA_001_d_Bacteria_p_Proteobacteria_c_Gammaproteobacteria
```
前缀含义：d=域（domain），p=门（phylum），c=纲（class），o=目（order），f=科（family），g=属（genus），s=种（species）

### 格式 B（分号分隔）
```
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria
```

### 外部表（TSV 或 CSV）
```
name	taxonomy
Sample1	d__Bacteria;p__Firmicutes
Sample2	d__Archaea;p__Euryarchaeota
```

分类学相关 CLI 选项：
- `--taxonomy-delimiter-mode`：解析模式，可选 `reverse`（默认）、`greedy` 或 `segment`
- `--taxonomy-source-priority`：数据来源优先级，可选 `table` 或 `embedded`
- `--taxonomy-levels`：自定义层级前缀，例如 `d,p,c,o,f,g,s`
- `--table-sep`：表字段分隔符（默认：自动检测）
