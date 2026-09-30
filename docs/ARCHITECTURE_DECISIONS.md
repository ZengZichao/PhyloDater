> **Language:** English (canonical) · [中文](ARCHITECTURE_DECISIONS_CN.md)

# PhyloDater Architecture Decision Record (ADR)

> This document records the key architectural decisions in the PhyloDater project, including decision context, choices, and consequences trade-offs.

---

## ADR-001: Adopting Adapter Pattern to Integrate Multiple Dating Engines

**Status**: Accepted (2024-01-15)

**Context**:
Users need to use different dating methods (e.g., MCMCTree, treePL, PATHd8, etc.) for comparative analysis, but each method's input format, output structure, and parameter configuration vary greatly. Handling these differences directly in the pipeline code would lead to:
- Code complexity growing linearly with the number of methods
- Adding new methods requires modifying the core pipeline
- Implementation details between methods leak into the business logic layer

**Decision**:
Introduce `DatingMethod` abstract base class and `DatingMethodRegistry` registry factory pattern.

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

**Choices**:
- **Option A (adopted)**: Abstract base class + registry factory
  - Pros: Type safety, compile-time validation, auto-discovery
  - Cons: Requires inheritance, less flexible

- **Option B**: Protocol + dependency injection
  - Pros: No inheritance required, more flexible
  - Cons: More complex type checking

- **Option C**: Configuration-driven + reflection
  - Pros: Zero-code addition of new methods
  - Cons: Runtime errors, performance overhead

**Consequences**:
- **Positive**: Adding a new dating method only requires implementing the interface and calling `DatingMethodRegistry.register()`, no need to modify the pipeline code
- **Negative**: Adapters may leak implementation details to the visualization layer (resolved by refactoring `ComparisonReporter`)
- **Negative**: Modifications to the abstract base class will affect all adapters

**Maintenance Strategy**:
- Once released, the `DatingMethod` interface will only add methods, never remove
- New optional methods provide default implementations
- Breaking changes are handled through version upgrades

---

## ADR-002: Standardizing the DatingResult Data Structure

**Status**: Accepted (2024-02-20)

**Context**:
Different dating engines have different output formats:
- MCMCTree: FigTree.tre format with HPD intervals
- treePL: Newick format with point estimates only
- LSD2: Supports confidence intervals

The visualization module needs a uniform way to access results.

**Decision**:
Define `DatingResult` and `NodeAgeEstimate` standardized data structures:

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

**Choices**:
- **Option A (adopted)**: Strongly-typed dataclass
  - Pros: Type safety, self-documenting, serialization support
  - Cons: Extension requires class definition modification

- **Option B**: Generic dictionary + JSON Schema
  - Pros: Flexible, JSON-native
  - Cons: Lacks type hints, runtime errors

**Constraints**:
- The `metadata` field stores engine-specific information, but the visualization module **should not** depend on its internal structure
- Each adapter declares the fields it supports via `ADAPTER_CAPABILITIES`

**Consequences**:
- **Positive**: Visualization module completely decoupled from engine implementation
- **Positive**: Supports result serialization/deserialization
- **Negative**: Some engine-specific information (e.g., treePL's optimal_smooth) needs to be exposed through extension fields or capability declaration mechanisms

---

## ADR-003: Layered Configuration Management (defaults < YAML < CLI)

**Status**: Accepted (2024-03-01)

**Context**:
Users need to configure parameters in different scenarios:
- Personal research: Fixed set of parameters
- Batch analysis: Shared configuration file
- Special requirements: Temporarily override specific parameters

**Decision**:
Three-tier configuration priority: CLI > YAML > defaults

```python
@dataclass
class PipelineConfig:
    methods: List[str] = field(default_factory=lambda: ['mcmctree', 'treepl', 'pathd8'])
    output_dir: Path = Path("./phylodater_output")
    threads: int = field(default_factory=lambda: get_env_threads(4))
```

Environment variables used for runtime overrides:
- `PHYLODATER_SEED`: Global random seed
- `PHYLODATER_THREADS`: Number of parallel threads
- `PHYLODATER_BURNIN`: Number of MCMC burnin generations
- `PHYLODATER_TIMEOUT`: Execution timeout

**Consequences**:
- **Positive**: Flexible configuration adapts to different usage scenarios
- **Positive**: Configuration files can be version-controlled
- **Positive**: Environment variables support containerization and automation

---

## ADR-004: Checkpoint-Resume Mechanism Design

**Status**: Accepted (2024-04-01)

**Context**:
Bayesian methods like MCMCTree can run for days to weeks. During long-running processes:
- HPC jobs may be terminated due to timeout
- Systems may crash due to insufficient resources
- Users may need to interrupt and re-run with modified parameters

**Decision**:
Introduce two levels of checkpointing:

1. **Adapter-level checkpoint** (`MCMCTreeCheckpointManager`):
   - Records MCMCTree internal steps (Hessian computation, MCMC sampling, etc.)
   - Supports recovery from the last completed step

2. **Pipeline-level checkpoint** (`PipelineCheckpointManager`):
   - Records the execution status of each method
   - Caches `DatingResult` to support cross-method recovery
   - Results persisted to JSON files

```python
class PipelineCheckpointManager:
    def get_methods_to_run(self, requested_methods: List[str]) -> List[str]:
        # Skip completed ones, keep failed ones pending retry

    def get_cached_result(self, method_name: str) -> Optional[DatingResult]:
        # Recover completed results from the checkpoint
```

**Choices**:
- **Option A (adopted)**: Result caching + state tracking
  - Pros: Supports full recovery, results are reusable
  - Cons: Requires persistent storage

- **Option B**: State tracking only
  - Pros: Simple
  - Cons: Recovery requires recomputation

**Consequences**:
- **Positive**: Reliability of long-running tasks greatly improved
- **Positive**: Supports incremental analysis (adding new methods without re-running existing ones)
- **Negative**: Checkpoint files require additional storage space
- **Negative**: Increased complexity in result version management

---

## ADR-005: Plugin System Architecture

**Status**: Accepted (2024-04-08)

**Context**:
Users may need to insert custom quality control or processing steps, such as:
- Long branch attraction sequence detection
- Outlier sequence filtering
- Custom calibration constraints

**Decision**:
Define plugin interfaces and integrate into the pipeline:

```python
class QCPlugin(Plugin):
    def process_tree(self, tree: Any, context: Dict[str, Any]) -> Any: pass

class PluginRegistry:
    _plugins: Dict[str, LoadedPlugin]

class PluginLoader:
    def discover_plugins(self, search_paths: List[Path]) -> List[Path]: pass
```

**Choices**:
- **Option A (adopted)**: Abstract base class + file system discovery
  - Pros: Type safety, auto-discovery
  - Cons: Requires inheritance

- **Option B**: Decorator + entry points
  - Pros: No inheritance required
  - Cons: Implicit dependencies

**Consequences**:
- **Positive**: Users can customize QC workflows
- **Positive**: Ecosystem extensibility
- **Negative**: Plugin quality is hard to control
- **Negative**: Introduces security considerations (pending code execution)

---

## ADR-006: Error Classification and Actionable Suggestions

**Status**: Accepted (2024-04-09)

**Context**:
Bioinformatics pipelines encounter various errors:
- User input errors (file format, parameter errors)
- Engine run failures (convergence failure, out of memory)
- Insufficient system resources (timeout, insufficient memory)

Users need actionable error information to quickly identify problems.

**Decision**:
Introduce `ErrorCategory` enum and structured error suggestions:

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

Helper function to classify errors based on external program output:

```python
def classify_execution_error(stderr: str, returncode: int) -> tuple:
    if "out of memory" in stderr_lower:
        return (ErrorCategory.RESOURCE_LIMITED, "请尝试减少数据规模...")
```

**Consequences**:
- **Positive**: Users can quickly determine the error type and solution
- **Positive**: Supports automated error handling strategies
- **Positive**: Supports structured logging and monitoring

---

## ADR-009: Per-Method Calibration Compatibility Pre-check

**Status**: Accepted (2026-07-12)

**Context**:
Different dating methods support calibration constraints differently. For example, MCMCTree's `RootAge` can only be set as an upper bound, not a true fixed root age. If a user supplies a `fixed` root constraint, MCMC sampling will converge to a value below it, producing misleading results. Previously such problems only surfaced mid-run, wasting compute and producing poor UX.

**Decision**:
Add a `validate_calibrations(calibrations)` hook to the `DatingMethod` abstract base class. Each adapter can perform method-specific calibration checks before `prepare_inputs()` is called.

```python
class DatingMethod(ABC):
    def validate_calibrations(self, calibrations: List[CalibrationPoint]) -> List[str]:
        """Return a list of error messages; an empty list means pass."""
        return []
```

The pipeline calls this hook after `validate_environment()` and before `prepare_inputs()`. If non-empty errors are returned, the method is marked failed and the pipeline continues with the remaining methods.

**Options considered**:
- **Option A (chosen)**: Base-class hook + adapter override
  - Pros: clear responsibility; new methods declare their own constraints by overriding the hook
  - Cons: each adapter must implement it actively

- **Option B**: Global calibration rule table
  - Pros: centralized management
  - Cons: hard to express complex logic (e.g. MCMCTree RootAge semantics); new rules require core code changes

- **Option C**: Rely entirely on external software errors
  - Pros: no extra code
  - Cons: late discovery, unfriendly messages, wasted compute

**Consequences**:
- **Positive**: Incompatible configurations are caught early, avoiding mid-run failures
- **Positive**: Error messages can include actionable suggestions (e.g. switch to uniform/maximum or use another method)
- **Positive**: Failure of one method does not stop other methods
- **Negative**: Validation logic must be maintained per adapter

**Update history**:
- 2026-07-12: ADR accepted; fixed-root check implemented in MCMCTree adapter.

---

## Pending Decisions

### ADR-007: Decoupling Visualization Module from Engines (in progress)

**Issue**: The current `DatingFigure` component directly parses FigTree.tre format, violating the ADR-002 design principles.

**Candidate Options**:
- **Option A**: Rely solely on `DatingResult`, enhance `NodeAgeEstimate` fields
- **Option B**: Introduce a "visualization adapter" abstraction layer
- **Option C**: Extend the `ADAPTER_CAPABILITIES` mechanism

### ADR-008: Batch Processing Architecture

**Issue**: Current batch processing requires external scripts to loop and call `DatingPipeline`.

**Candidate Options**:
- **Option A**: Extend `PipelineConfig` to support multi-sample configuration
- **Option B**: Standalone `BatchProcessor` class (already implemented)
- **Option C**: Workflow language integration (Snakemake/Nextflow)

---

## Maintenance Guide

### How to Add a New ADR

1. Add a new section at the end of this document
2. Include: Status, Context, Decision, Choices, Consequences
3. Status options: Proposed, Accepted, Accepted, Superseded, Deprecated

### How to Update an Existing ADR

1. Add an "Update Record" paragraph at the end of the ADR
2. Explain the reasons and impacts of the update
3. Keep the decision history intact

### ADR Review Cycle

- Review all "Pending" ADRs before major version releases
- Evaluate whether a new ADR is needed when adding external dependencies
- Periodically clean up deprecated ADRs
