# Transcriptomics Dashboard - Complete Setup Guide

## Phase 1 Implementation - Setup Instructions

This guide will walk you through setting up the complete Phase 1 (MVP) implementation.

## Prerequisites

Before starting, ensure you have:

- **Python 3.9+** installed
- **R 4.0+** installed
- **Git** (for cloning the repository)
- **4GB+ RAM** available
- **Command line/terminal** access

## Step-by-Step Setup

### 1. Create Project Directory Structure

```bash
mkdir transcriptomics-dashboard
cd transcriptomics-dashboard

# Create the source directory structure
mkdir -p src/transcriptomics_dashboard
mkdir -p tests
mkdir -p examples
mkdir -p docs
mkdir -p data/reference
```

### 2. Create Package Files

Create the following directory structure:

```
transcriptomics-dashboard/
├── src/
│   └── transcriptomics_dashboard/
│       ├── __init__.py
│       ├── app.py
│       ├── config.py
│       ├── validation.py
│       ├── deseq2.py
│       └── visualizations.py
├── tests/
│   ├── __init__.py
│   └── test_validation.py
├── examples/
│   └── generate_example_data.py
├── pyproject.toml
├── environment.yml
├── README.md
├── SETUP.md
└── .gitignore
```

Create `src/transcriptomics_dashboard/__init__.py`:

```python
"""Transcriptomics Dashboard - Interactive RNA-seq analysis tool."""

__version__ = "0.1.0"
__author__ = "Alec Forston"

from .config import get_config, Config
from .validation import validate_count_matrix, validate_metadata
from .deseq2 import run_deseq2

__all__ = [
    'get_config',
    'Config',
    'validate_count_matrix',
    'validate_metadata',
    'run_deseq2'
]
```

### 3. Install R and DESeq2

#### On macOS:
```bash
# Install R
brew install r

# Start R and install DESeq2
R
```

#### On Ubuntu/Debian:
```bash
# Install R
sudo apt-get update
sudo apt-get install r-base r-base-dev

# Start R
R
```

#### On Windows:
Download and install R from [r-project.org](https://www.r-project.org/)

#### Install DESeq2 (in R):
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

# Verify installation
library(DESeq2)
packageVersion("DESeq2")

# Exit R
quit()
```

### 4. Set Up Python Environment

#### Option A: Using Conda (Recommended)

```bash
# Create environment from file
conda env create -f environment.yml

# Activate environment
conda activate transcriptomics

# Install package in development mode
pip install -e .
```

#### Option B: Using venv and pip

```bash
# Create virtual environment
python -m venv venv

# Activate environment
# On macOS/Linux:
source venv/bin/activate
# On Windows:
venv\Scripts\activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

Create `requirements.txt`:
```
dash>=2.14.0
dash-bootstrap-components>=1.5.0
plotly>=5.18.0
pandas>=2.0.0
numpy>=1.24.0
scipy>=1.11.0
rpy2>=3.5.0
scikit-learn>=1.3.0
umap-learn>=0.5.0
goatools>=1.3.0
gseapy>=1.0.0
openpyxl>=3.1.0
pydantic>=2.0.0
pydantic-settings>=2.0.0
pyyaml>=6.0.0
joblib>=1.3.0
pytest>=7.4.0
pytest-cov>=4.1.0
```

### 5. Generate Example Data

```bash
# Generate example datasets
python examples/generate_example_data.py
```

This creates:
- `examples/sample_counts.csv` - Count matrix
- `examples/sample_metadata.csv` - Sample metadata
- `examples/ground_truth.csv` - True DE status

### 6. Run Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src/transcriptomics_dashboard --cov-report=html

# View coverage report
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

### 7. Launch the Dashboard

```bash
# Start the application
transcriptomics-dashboard

# Or run directly
python -m transcriptomics_dashboard.app
```

The dashboard will be available at: `http://127.0.0.1:8050`

### 8. Test the Application

1. Open your browser to `http://127.0.0.1:8050`
2. Upload `examples/sample_counts.csv` as the count matrix
3. Upload `examples/sample_metadata.csv` as the metadata
4. Select "condition" from the condition column dropdown
5. Choose comparison: "treatment" vs "control"
6. Adjust thresholds if desired (FDR: 0.05, log2FC: 1.0)
7. Click "Run Analysis"
8. Explore the visualizations in the tabs

## Troubleshooting

### rpy2 Installation Issues

**Problem**: `rpy2` fails to install

**Solution**:
```bash
# macOS
brew install r
export R_HOME=/Library/Frameworks/R.framework/Resources

# Ubuntu/Debian
sudo apt-get install r-base-dev
export R_HOME=/usr/lib/R

# Then retry
pip install rpy2
```

### DESeq2 Not Found

**Problem**: "DESeq2 not found" error

**Solution**:
```r
# In R, reinstall DESeq2
BiocManager::install("DESeq2", force = TRUE)

# Check installation path
.libPaths()
```

### Memory Issues

**Problem**: Application crashes with large datasets

**Solution**:
- Increase system memory
- Filter low-count genes before analysis
- Adjust config: `~/.transcriptomics_dashboard/config.yaml`

### Port Already in Use

**Problem**: Port 8050 is already in use

**Solution**:
```bash
# Use a different port
export TRANS_PORT=8051
transcriptomics-dashboard
```

Or edit `config.yaml`:
```yaml
port: 8051
```

## Configuration

The application creates a configuration directory at:
`~/.transcriptomics_dashboard/`

Structure:
```
~/.transcriptomics_dashboard/
├── config.yaml           # User configuration
├── data/
│   └── reference/        # GO/KEGG databases (Phase 2)
├── sessions/             # Saved sessions (Phase 2)
└── cache/                # Computation cache
```

Edit `config.yaml` to customize:
```yaml
defaults:
  fdr_threshold: 0.05
  log2fc_threshold: 1.0
  threads: -1

app_title: My Transcriptomics Dashboard
port: 8050
debug: false
```

## Development Setup

### Install Development Dependencies

```bash
# Install pre-commit hooks
pip install pre-commit
pre-commit install

# Install development packages
pip install -e ".[dev]"
```

### Code Style

The project uses:
- **Black** for code formatting
- **isort** for import sorting
- **flake8** for linting

Run formatters:
```bash
# Format code
black src/ tests/

# Sort imports
isort src/ tests/

# Check style
flake8 src/ tests/
```

### Running Tests in Development

```bash
# Run tests on file save (watch mode)
pytest-watch

# Run specific test file
pytest tests/test_validation.py -v

# Run specific test
pytest tests/test_validation.py::TestCountMatrixValidation::test_valid_count_matrix
```

## Next Steps

After successful Phase 1 setup:

1. **Validate Installation**: Run all tests and example analysis
2. **Customize**: Modify config.yaml for your needs
3. **Prepare Data**: Format your own RNA-seq data
4. **Phase 2**: Implement enrichment analysis and PCA

## Phase 1 Checklist

- [ ] R and DESeq2 installed
- [ ] Python environment set up
- [ ] All dependencies installed
- [ ] Example data generated
- [ ] Tests pass
- [ ] Dashboard launches successfully
- [ ] Example analysis runs without errors
- [ ] All visualizations display correctly

## Getting Help

### Common Resources

- **Documentation**: Check the README.md for usage examples
- **Tests**: Look at test files for API examples
- **Examples**: Review generate_example_data.py for data formats

### Debugging

Enable debug mode:
```yaml
# In config.yaml
debug: true
```

View logs:
```bash
# Run with verbose logging
python -m transcriptomics_dashboard.app --log-level DEBUG
```

## Project Statistics

**Phase 1 Completion:**
- ✓ 7 major components implemented
- ✓ ~2,500 lines of Python code
- ✓ 4 visualization types
- ✓ Full data validation pipeline
- ✓ Complete DESeq2 integration
- ✓ Interactive web interface

**Estimated Development Time:** 2-3 weeks

**Next Phase:** Core Features (enrichment analysis, dimensionality reduction, session management)