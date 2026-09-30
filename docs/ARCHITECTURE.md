> **Language:** English (canonical) · [中文](ARCHITECTURE_CN.md)

# PhyloDater Architecture Document

This document describes the software architecture design of PhyloDater.

## Architecture Overview

PhyloDater adopts a layered architecture design, following the separation of concerns principle:

```
┌─────────────────────────────────────────────────────────────┐
│                        CLI Layer                            │
│               (Command-Line Interface)                      │
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

## Layer Details

### 1. CLI Layer (Command-Line Interface Layer)

**Responsibility**: Provides user interaction interface

**Main Components**:
- `cli/main.py`: Command-line argument parsing and main flow control

**Design Principles**:
- Single Responsibility: Only handles user input and output display
- Depends on lower layers: Executes business logic through Service Layer and Core Layer

### 2. Service Layer

**Responsibility**: Encapsulates business logic and coordinates domain objects

**Main Components**:
- `calibration_loader.py`: Calibration point configuration loading
- `calibration_resolver.py`: Calibration point resolution (MRCA location)
- `tree_validator.py`: Tree structure validation
- `tree_rooting.py`: Tree rooting strategy
- `name_mapping.py`: Taxonomic name mapping management
- `taxonomy_parser.py`: Taxonomy data parsing
- `input_validator.py`: Input file format validation
- `deep_validator.py`: Deep validation (parenthesis balance, negative branches, malicious characters, etc.)
- `auto_calibrator.py`: Automatic calibration point addition

**Design Principles**:
- High Cohesion: Each service focuses on a specific business function
- Testability: Business logic is independent of external dependencies

### 3. Core Layer

**Responsibility**: Defines core business logic and interfaces

**Main Components**:
- `method_interface.py`: Dating method abstract interface and registry
- `pipeline.py`: Workflow engine (sequential and parallel)
- `comparison_reporter.py`: Comparison report generation
- `exceptions.py`: Exception hierarchy

**Key Design Patterns**:
- **Strategy Pattern**: `DatingMethod` interface allows unified invocation of different dating algorithms
- **Registry Pattern**: `DatingMethodRegistry` implements dynamic method discovery and creation
- **Factory Pattern**: Creates concrete method instances through the registry
- **Template Method Pattern**: Optional hooks such as `DatingMethod.validate_calibrations()` provide default implementations and allow adapters to add method-specific validation

### 4. Adapter Layer

**Responsibility**: Integrates external dating software

**Main Components**:
- `mcmctree_method.py`: PAML MCMCTree adapter
- `lsd2_method.py`: IQ-TREE2 LSD2 adapter
- `r8s_pyr8s_method.py`: r8s/pyr8s adapter
- `treepl_method.py`: treePL adapter
- `pathd8_method.py`: PATHd8 adapter
- `wlogdate_method.py`: wLogDate adapter
- `mdcat_method.py`: MD-Cat adapter

**Design Principles**:
- **Adapter Pattern**: Unifies interfaces from different software to the `DatingMethod` interface
- **Single Responsibility**: Each adapter only handles one external software
- **Error Isolation**: Errors in external software do not affect other components

### 5. Model Layer

**Responsibility**: Defines domain models and data structures

**Main Components**:
- `constraints.py`: Age constraint type hierarchy
- `calibration.py`: Calibration point related models
- `results.py`: Dating result models
- `tree.py`: Phylogenetic tree models

**Design Principles**:
- **Value Object**: Constraint types implemented as immutable objects using dataclass
- **Rich Domain Model**: Models contain business logic (e.g., format conversion)

### 6. Infrastructure Layer

**Responsibility**: Provides technical capabilities and cross-cutting concerns

**Main Components**:
- `configuration.py`: Three-tier configuration management (defaults < YAML < CLI)
- `logging.py`: Logging system
- `process_runner.py`: Process execution and version detection
- `checkpoint.py`: Checkpoint management (MCMCTree support)
- `auto_naming.py`: Output file automatic naming
- `alignment_metadata.py`: Sequence alignment metadata extraction

**Design Principles**:
- **Dependency Inversion**: Upper layers depend on abstract interfaces, not concrete implementations
- **Cross-Cutting Concerns**: Logging, configuration and other functions are centrally managed

## Data Flow

```
Input files
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
                       Output files
```

## Key Design Decisions

### 1. Constraint Type Hierarchy

**Decision**: Use dataclass to implement age constraint type hierarchy

**Rationale**:
- Type safety: Compile-time checking
- Immutability: Prevents accidental modification
- Self-documenting: Code is documentation

### 2. Adapter Registration Mechanism

**Decision**: Use class-level registry rather than configuration file

**Rationale**:
- Type safety: Compile-time validation of adapter existence
- Auto-discovery: New adapters automatically register
- Zero configuration: No need to maintain configuration files

### 3. Three-Tier Configuration Management

**Decision**: Built-in defaults < YAML configuration < Command-line arguments

**Rationale**:
- Flexibility: Supports multiple configuration methods
- Predictability: Clear priority
- Reproducibility: Configuration files can be version-controlled

### 4. Exception Hierarchy

**Decision**: Define clear exception type hierarchy

**Rationale**:
- Precise error handling: Callers can catch specific exceptions
- Error classification: Validation errors, execution errors, result errors
- Semantic degradation warnings: Non-fatal issues use warnings rather than exceptions

### 5. Adapter Configuration Typing

**Decision**: the abstract base class is generic in the configuration it reads —
`DatingMethod[ConfigT]` — and each adapter pins that parameter to its own config
dataclass (`TreePLMethod(DatingMethod[TreePLConfig])`, `MCMCTreeMethod(DatingMethod[MCMCTreeConfig])`, …).

**Problem it solves**: adapters historically accepted *either* a whole `ToolConfig`
or a single method's sub-config (e.g. `TreePLConfig`). They distinguished the two at
runtime with `hasattr(config, "treepl")`, then assigned the sub-config to
`self.config`. The base class declared `config: ToolConfig`, so every
`self.config.treepl_bin` / `self.config.numsites` access was, to a static reader,
an attribute that `ToolConfig` does not have — 223 `attr-defined` errors, and no way
to tell a real typo apart from the duck typing.

**Why the two-form entry point was kept**: both forms are used heavily by the test
suite (`MCMCTreeMethod(MCMCTreeConfig(), …)` and `WLogDateMethod(ToolConfig(), …)`),
so collapsing the entry to `ToolConfig` only would have meant rewriting those call
sites and their assertions. Instead each adapter normalises at the door:

```python
if isinstance(config, ToolConfig):      # equivalent to the old hasattr test:
    method_config = config.treepl       # ToolConfig always has .treepl,
                                        # TreePLConfig never does
else:
    method_config = config
super().__init__(method_config, output_dir, software_paths, common_config=common)
```

`self.config` is therefore exactly `TreePLConfig` for every reader, and no
`self.config: Any` escape hatch was used. `DatingMethod.__init__` also accepts
`common_config` and normalises it to a non-`Optional` `CommonConfig`, which is what
the fifteen `self.common_config.nthreads / .timeout / .seed` accesses needed.

**Related shared type**: `InputFileInfo = Union[Path, int, float]` — the
`prepare_inputs()` "input files" dictionary legitimately holds both file paths and
per-run statistics (site count, wall time, log-likelihood). `DatingMethod.cleanup()`
skips the non-path entries; before this was typed it called `Path(value)` on them,
and `Path(1000)` raises `TypeError`, which is not in `cleanup()`'s `except` clause.

**Rejected alternative**: `self.config: Any = config`. It silences the same 223
errors but deletes every check on adapter configuration access, which is exactly
the surface where a misspelled engine parameter goes silently missing.

## Extension Guide

### Adding a New Dating Method

1. Create a new adapter class in `phylodater/adapters/`
2. Inherit from `DatingMethod` abstract base class
3. Implement all abstract methods
4. Register at end of file: `DatingMethodRegistry.register("method_name", YourClass)`
5. Add unit tests

### Adding a New Constraint Type

1. Create a new dataclass in `phylodater/models/constraints.py`
2. Inherit from `AgeConstraint` abstract base class
3. Implement `to_software_format()` and `to_mcmctree_calib_string()`
4. Add validation logic in `__post_init__()`
5. Add unit tests

## Performance Considerations

- **Parallel Execution**: `ParallelPipeline` supports parallel running of multiple methods
- **Checkpoint**: MCMCTree supports checkpoint recovery, avoiding repeated computation
- **Memory Management**: Streaming processing for large datasets
- **Caching**: Configuration and metadata caching to avoid repeated computation

## Security Considerations

- **Input Validation**: All external input is validated
- **Command Injection Prevention**: Use lists instead of strings to pass command arguments
- **Resource Limits**: Support timeout and memory limits

## Testing Strategy

- **Unit Tests**: Each module tested independently
- **Integration Tests**: Test interactions between modules
- **Adapter Tests**: Use mock to test external software integration
- **End-to-End Tests**: Complete workflow testing
