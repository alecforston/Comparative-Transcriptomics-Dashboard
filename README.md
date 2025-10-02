# Transcriptomics Dashboard

A local interactive web application for RNA-seq differential expression analysis using DESeq2.

![Version](https://img.shields.io/badge/version-0.1.0-blue)
![Python](https://img.shields.io/badge/python-3.9+-green)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

## Features

### Phase 1 (MVP) - Completed ✓

- **Easy Data Upload**: Drag-and-drop interface for count matrices and metadata files
- **Data Validation**: Comprehensive validation with helpful error messages
- **DESeq2 Integration**: Full differential expression analysis using the gold-standard DESeq2 method
- **Interactive Visualizations**:
  - Volcano plots with gene labeling
  - Expression heatmaps with hierarchical clustering
  - MA plots
  - PCA plots
- **Results Export**: Download results as CSV files

## Installation

### Prerequisites

- Python 3.9 or higher
- R (≥4.0.0) with DESeq2 installed
- At least 4GB RAM

### Step 1: Install R and DESeq2

First, install R from [r-project.org](https://www.r-project.org/)

Then install DESeq2 in R:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

### Step 2: Install Python Package

#### Using pip:

```bash
pip install transcriptomics-dashboard
```

#### From source:

```bash
git clone https://github.com/yourusername/transcriptomics-dashboard.git
cd transcriptomics-dashboard
pip install -e .
```

#### Using conda (recommended):

```bash
conda env create -f environment.yml
conda activate transcriptomics
pip install -e .
```

## Quick Start

### 1. Launch the Dashboard

```bash
transcriptomics-dashboard
```

The dashboard will open at `http://127.0.0.1:8050`

### 2. Upload Your Data

**Count Matrix Format** (CSV, TSV, or Excel):
```
           Sample1  Sample2  Sample3  Sample4
Gene1         120      150       80      110
Gene2          45       55       40       50
Gene3        1200     1100     1300     1250
...
```

**Metadata Format**:
```
           condition    batch
Sample1    control      1
Sample2    control      1
Sample3    treatment    2
Sample4    treatment    2
```

### 3. Run Analysis

1. Upload count matrix and metadata files
2. Select the condition column (e.g., "condition")
3. Choose comparison groups (e.g., treatment vs control)
4. Adjust FDR and log2FC thresholds
5. Click "Run Analysis"

### 4. Explore Results

- **Volcano Plot**: Interactive plot showing differential expression
- **Heatmap**: Top 50 differentially expressed genes
- **MA Plot**: Expression-dependent fold changes
- **PCA**: Sample clustering and quality control
- **Results Table**: Full results with statistics

## Input File Formats

### Count Matrix

- **Supported formats**: CSV, TSV, Excel (.xlsx, .xls)
- **Structure**: Genes as rows, samples as columns
- **Requirements**:
  - Non-negative integer counts
  - No missing values
  - Unique gene IDs

### Metadata

- **Supported formats**: CSV, TSV, Excel
- **Structure**: Samples as rows, conditions/variables as columns
- **Requirements**:
  - Sample names must match count matrix columns
  - At least 2 replicates per condition
  - No missing values in condition column

## Configuration

Configuration file is automatically created at: `~/.transcriptomics_dashboard/config.yaml`

```yaml
defaults:
  fdr_threshold: 0.05
  log2fc_threshold: 1.0
  min_count: 10
  threads: -1

paths:
  user_home: ~/.transcriptomics_dashboard

performance:
  max_cache_size_gb: 10.0
  parallel_backend: loky

app_title: Transcriptomics Dashboard
port: 8050
```

## Project Structure

```
transcriptomics_dashboard/
├── src/transcriptomics_dashboard/
│   ├── __init__.py
│   ├── app.py                    # Main Dash application
│   ├── config.py                 # Configuration management
│   ├── validation.py             # Data validation
│   ├── deseq2.py                 # DESeq2 wrapper
│   └── visualizations.py         # Plotting functions
├── tests/
│   ├── test_validation.py
│   ├── test_deseq2.py
│   └── test_visualizations.py
├── docs/
├── examples/
│   ├── sample_counts.csv
│   └── sample_metadata.csv
├── pyproject.toml
├── environment.yml
└── README.md
```

## Troubleshooting

### R Package Issues

If you see "DESeq2 not found":

```bash
# Check R is available
which R

# Install DESeq2 in R
R -e "BiocManager::install('DESeq2')"
```

### rpy2 Installation

If rpy2 fails to install:

```bash
# macOS
brew install r

# Ubuntu/Debian
sudo apt-get install r-base r-base-dev

# Then install rpy2
pip install rpy2
```

### Memory Issues

For large datasets (>50K genes):

1. Increase cache size in config.yaml
2. Use fewer genes (filter low counts)
3. Run on machine with more RAM

## Examples

### Example Data

The repository includes example datasets:

- `examples/sample_counts.csv`: 1000 genes × 6 samples
- `examples/sample_metadata.csv`: Sample annotations

### Basic Analysis Script

```python
from transcriptomics_dashboard.validation import read_count_matrix, read_metadata
from transcriptomics_dashboard.deseq2 import run_deseq2
from transcriptomics_dashboard.visualizations import create_volcano_plot

# Load data
counts = read_count_matrix('counts.csv')
metadata = read_metadata('metadata.csv')

# Run DESeq2
results = run_deseq2(
    counts,
    metadata,
    condition_col='condition',
    contrast=('treatment', 'control'),
    fdr_threshold=0.05,
    lfc_threshold=1.0
)

# Create volcano plot
fig = create_volcano_plot(results['results'])
fig.show()

# Save results
results['results'].to_csv('de_results.csv', index=False)
```

## Roadmap

### Phase 2: Core Features (In Progress)
- PCA and UMAP dimensionality reduction
- GO and KEGG enrichment analysis
- Session save/load functionality
- Additional diagnostic plots

### Phase 3: Polish & UX
- Batch processing CLI
- Automated HTML reports
- Performance optimizations
- Comprehensive documentation

### Phase 4: Advanced Features
- GSEA implementation
- Multi-experiment comparison
- Batch effect correction
- Gene network visualization

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Ensure all tests pass
5. Submit a pull request

## Citation

If you use this tool in your research, please also cite DESeq2:

```
Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and 
dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550.
```

## License

MIT License - see LICENSE file for details

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/transcriptomics-dashboard/issues)
- **Documentation**: [Read the Docs](https://transcriptomics-dashboard.readthedocs.io)
- **Email**: your.email@example.com

## Acknowledgments

- Built with [Dash](https://dash.plotly.com/) by Plotly
- Uses [DESeq2](https://bioconductor.org/packages/DESeq2/) for differential expression
- Inspired by tools like iDEP and DEBrowser