# CLAUDE.md - AI Assistant Guide for qPCR- Repository

## Project Overview

**Repository**: qPCR-
**Purpose**: Quantitative PCR (qPCR) data analysis project
**Status**: Initial setup - repository initialized but awaiting source code

## Repository Structure

```
qPCR-/
├── CLAUDE.md          # This file - AI assistant guidelines
└── .git/              # Git version control
```

> **Note**: This repository is currently empty. The structure above will be updated as the codebase develops.

## Project Context

qPCR (quantitative Polymerase Chain Reaction) projects typically involve:
- Processing raw fluorescence data from PCR machines
- Calculating Ct/Cq values (cycle threshold/quantification cycle)
- Performing relative or absolute quantification
- Statistical analysis of gene expression data
- Visualization of amplification curves and results

## Development Guidelines

### Getting Started

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd qPCR-
   ```

2. Set up development environment (to be defined based on chosen tech stack)

### Recommended Tech Stack (Suggestions)

For qPCR analysis projects, consider:

**Python-based**:
- `pandas` - Data manipulation
- `numpy` - Numerical computations
- `scipy` - Statistical analysis
- `matplotlib`/`seaborn` - Visualization
- `openpyxl` - Excel file handling (common qPCR export format)

**R-based**:
- `qpcR` - qPCR analysis package
- `tidyverse` - Data manipulation
- `ggplot2` - Visualization

### Code Conventions

#### File Organization
- `/src` or `/qpcr` - Main source code
- `/tests` - Unit and integration tests
- `/data` - Sample data files (add to .gitignore if sensitive)
- `/docs` - Documentation
- `/scripts` - Utility scripts
- `/notebooks` - Jupyter notebooks for analysis

#### Naming Conventions
- Use `snake_case` for Python files and functions
- Use descriptive variable names (e.g., `ct_values`, `reference_gene`, `fold_change`)
- Prefix private functions with underscore (`_internal_function`)

#### Documentation
- Include docstrings for all public functions
- Document input/output formats for data processing functions
- Maintain a changelog for version tracking

### Git Workflow

#### Branch Naming
- `main` - Production-ready code
- `develop` - Integration branch
- `feature/<description>` - New features
- `fix/<description>` - Bug fixes
- `claude/<session-id>` - AI-assisted development branches

#### Commit Messages
Follow conventional commits:
```
type(scope): brief description

- Detailed explanation if needed
- List specific changes
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

### Testing

- Write unit tests for all data processing functions
- Include test data files with known expected outputs
- Validate calculations against established qPCR analysis tools

### Data Handling

#### Input Formats (Common)
- `.csv` - Comma-separated values
- `.xlsx` - Excel spreadsheets
- `.txt` - Tab-delimited text
- Proprietary formats (Bio-Rad, Applied Biosystems, etc.)

#### Output Considerations
- Support export to common formats (CSV, Excel)
- Include metadata in output files
- Preserve traceability to raw data

## AI Assistant Instructions

### When Working on This Repository

1. **Understand the domain**: qPCR analysis requires understanding of:
   - Amplification curves and fluorescence data
   - Ct/Cq value calculation methods
   - Relative quantification (ΔΔCt method)
   - Reference gene normalization
   - Statistical significance testing

2. **Prioritize accuracy**: qPCR data is used for scientific conclusions
   - Double-check calculations
   - Include validation steps
   - Document assumptions and limitations

3. **Handle data carefully**:
   - Never commit sensitive research data
   - Use sample/mock data for testing
   - Respect data privacy requirements

4. **Code quality**:
   - Write clean, readable code
   - Include comprehensive error handling
   - Validate input data formats
   - Provide informative error messages

### Common Tasks

- **Data import**: Parse various qPCR machine output formats
- **Ct calculation**: Implement threshold-based or second derivative methods
- **Normalization**: Apply reference gene corrections
- **Statistical analysis**: Calculate fold changes, p-values, confidence intervals
- **Visualization**: Generate amplification curves, bar charts, heatmaps

### Key Formulas Reference

**Relative Quantification (ΔΔCt method)**:
```
ΔCt = Ct(target) - Ct(reference)
ΔΔCt = ΔCt(sample) - ΔCt(control)
Fold Change = 2^(-ΔΔCt)
```

**Efficiency Correction**:
```
Ratio = E(target)^ΔCt(target) / E(reference)^ΔCt(reference)
```

## Quick Reference

| Task | Command/Action |
|------|----------------|
| Run tests | `pytest tests/` (Python) or `Rscript tests/run_tests.R` (R) |
| Format code | `black .` (Python) or `styler::style_dir()` (R) |
| Lint | `flake8 src/` (Python) |
| Build docs | `sphinx-build docs/ docs/_build` |

## Resources

- [MIQE Guidelines](https://rdml.org/miqe.html) - Minimum Information for qPCR Experiments
- [RDML Format](https://rdml.org/) - Real-time PCR Data Markup Language
- qPCR analysis best practices literature

---

*This CLAUDE.md will be updated as the project develops. Last updated: 2026-01-26*
