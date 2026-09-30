> **Language:** English (canonical) · [中文](ERROR_CODES_CN.md)

# PhyloDater Error Codes

| Code | Name | Description |
|------|------|-------------|
| 0 | EXIT_SUCCESS | Successful completion |
| 1 | EXIT_RUNTIME_ERROR | Runtime error (missing dependency, self-test failure) |
| 2 | EXIT_USAGE_ERROR | Command-line usage error (invalid arguments, missing required options) |
| 3 | EXIT_DATA_ERROR | Input data error (invalid format, file not found, negative branch lengths, PHYLIP not supported, empty calibration, label mismatch) |
| 130 | EXIT_INTERRUPTED | User interrupt (Ctrl+C) |

## Common Error Code 3 Scenarios

- **File not found**: `--tree`, `--sequence`, `--calibrations` path does not exist
- **Negative branch lengths**: Tree contains branches with negative lengths (use `--skip-length-check` to bypass)
- **PHYLIP format**: PHYLIP alignment format is not supported; convert to FASTA (use `--ignore-malformed` to continue with warning)
- **Empty calibration**: Calibration file loaded 0 valid constraints
- **Label mismatch**: Tree and sequence file terminal labels do not match
- **Input too large**: File exceeds 500 MB limit
- **Invalid tree structure**: Unbalanced parentheses, duplicate names
- **Unsupported format**: Tree or alignment format not recognized
- **Circular taxonomy dependencies**: Taxonomy table contains circular references

## Method-Specific Calibration Errors

Some dating methods impose additional restrictions on calibration constraints. PhyloDater checks these requirements up front via each adapter's `validate_calibrations()` hook, preventing mid-run failures.

### MCMCTree Does Not Support Fixed Root Age

**Example error message**:
```
Method 'mcmctree' calibration incompatible:
  - MCMCTree does not support a fixed root age. Root calibration 'Root' is fixed=500.0 Ma,
    but MCMCTree can only use it as an upper bound (RootAge='<500.0'), so the root age cannot be exactly fixed.
  - Suggested fixes (choose one):
    1) Change the root constraint to uniform, e.g. min=499.9, max=500.1;
    2) Change the root constraint to maximum, e.g. max=500.0;
    3) If an exact fixed root is required, use treePL, PATHd8, or LSD2 instead.
```

**Cause**: MCMCTree's `RootAge` control-file parameter can only specify an upper bound (e.g. `<500`), not a true fixed age. If a `fixed` constraint is supplied, MCMC sampling will converge to a root age below that value, producing unexpected results.

**Resolution**:
1. Change the root constraint to `uniform`, e.g. `min=499.9, max=500.1`;
2. Change the root constraint to `maximum`, e.g. `max=500.0`;
3. If a strictly fixed root is required, use treePL, PATHd8, or LSD2.

> **Tip**: Different dating methods may use different calibration files. For example, MCMCTree can use a `maximum` root constraint while other methods continue to use a `fixed` root constraint.

## Exception Classes

PhyloDater uses structured exception classes for precise error handling:

| Exception Class | Parent | Exit Code | Description |
|----------------|--------|-----------|-------------|
| `PhyloDaterError` | Exception | - | Base exception for all PhyloDater errors |
| `ValidationError` | PhyloDaterError | 3 | Base class for input validation errors |
| `TreeValidationError` | ValidationError | 3 | Tree structure validation failure |
| `AlignmentValidationError` | ValidationError | 3 | Sequence alignment validation failure |
| `ConsistencyError` | ValidationError | 3 | Tree-sequence label mismatch |
| `PhyloFormatError` | ValidationError | 3 | Unsupported file format or format error |
| `NegativeBranchLengthError` | TreeValidationError | 3 | Tree contains negative branch lengths |
| `InputSizeLimitError` | ValidationError | 3 | Input file exceeds size limit (500MB) |
| `CalibrationError` | PhyloDaterError | 3 | Calibration configuration error |
| `CalibrationConflictError` | CalibrationError | 3 | Conflicting calibration constraints |
| `CalibrationResolutionError` | CalibrationError | 3 | Calibration point resolution failure |
| `TaxonomyConflictError` | CalibrationError | 3 | Circular dependencies in taxonomy table |
| `MonophylyError` | CalibrationError | 3 | Taxa do not form monophyletic clade |
| `ConfigurationError` | PhyloDaterError | 2 | Configuration file or parameter error |
| `UnknownMethodError` | PhyloDaterError | 2 | Requested dating method not registered |
| `ExecutionError` | PhyloDaterError | 1 | External software execution failure |
| `ExecutionTimeoutError` | ExecutionError | 1 | External software timed out |
| `ResultParsingError` | PhyloDaterError | 1 | Failed to parse external software output |

### Warning Classes

Non-fatal issues use warning classes (logged but do not terminate):

| Warning Class | Parent | Description |
|---------------|--------|-------------|
| `ConvergenceWarning` | PhyloDaterWarning | MCMC convergence diagnostic threshold exceeded |
| `SemanticDegradationWarning` | PhyloDaterWarning | Input semantics degraded (e.g., constraint type downgraded) |
| `EnvironmentWarning` | PhyloDaterWarning | Environment issue (e.g., optional dependency missing) |
| `ResourceWarning` | PhyloDaterWarning | Resource usage warning (e.g., large tree detected) |
| `UnknownParameterWarning` | PhyloDaterWarning | Unknown method-args parameter encountered |

## Taxonomy Table Format

PhyloDater supports two taxonomy annotation formats:

### Format A (Embedded in tip names)
```
GB_GCA_001_d_Bacteria_p_Proteobacteria_c_Gammaproteobacteria
```
Prefixes: d=domain, p=phylum, c=class, o=order, f=family, g=genus, s=species

### Format B (Semicolon-separated)
```
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria
```

### External Table (TSV/CSV)
```
name	taxonomy
Sample1	d__Bacteria;p__Firmicutes
Sample2	d__Archaea;p__Euryarchaeota
```

CLI options for taxonomy:
- `--taxonomy-delimiter-mode`: `reverse` (default), `greedy`, or `segment`
- `--taxonomy-source-priority`: `table` or `embedded`
- `--taxonomy-levels`: Custom level prefixes, e.g. `d,p,c,o,f,g,s`
- `--table-sep`: Table field separator (default: auto-detect)
