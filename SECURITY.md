> **Language:** English (canonical) · [中文](SECURITY_CN.md)

# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in PhyloDater, please report it
responsibly.

**Do NOT open a public GitHub issue for security vulnerabilities.**

Instead, please email the maintainer at **zengzichao@sjtu.edu.cn** with:

1. A description of the vulnerability
2. Steps to reproduce or a proof-of-concept
3. The potential impact
4. Any suggested fixes (optional)

You should receive an acknowledgement within 72 hours. If you do not receive
a response, please follow up.

## Scope

PhyloDater processes phylogenetic tree files and configuration YAML. The
following are **in scope**:

- Code execution through crafted input files (e.g., Newick injection)
- Path traversal via configuration parameters
- Plugin loading from untrusted paths

The following are **out of scope**:

- Vulnerabilities in third-party dependencies (report to upstream maintainers)
- Social engineering attacks
- DoS through extremely large input files (use `--max-input-size`)

## Plugin Security

PhyloDater's plugin system only loads Python files from explicitly allowed
directories (see `infrastructure/plugins.py`). Plugins from arbitrary paths
are rejected by default. Do not disable this safety check.

## Disclosure Policy

- We will acknowledge receipt of your report within 72 hours.
- We will investigate and provide an initial assessment within 7 days.
- We will coordinate with you on the disclosure timeline.
- We prefer coordinated disclosure after a fix is available.
