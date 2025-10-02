# Example Data

This directory contains synthetic RNA-seq datasets for testing the Transcriptomics Dashboard.

## Files

### Standard Dataset
- `sample_counts.csv` - Count matrix (1000 genes × 6 samples)
- `sample_metadata.csv` - Sample metadata with conditions
- `ground_truth.csv` - True differential expression status (for validation)

Available in multiple formats:
- CSV: `.csv`
- TSV: `.tsv`
- Excel: `.xlsx`

### Dataset Characteristics

**Samples:**
- 3 control samples
- 3 treatment samples

**Genes:**
- Total: 1,000 genes
- Differentially expressed: 100 genes
  - Up-regulated in treatment: 50 genes
  - Down-regulated in treatment: 50 genes
- Non-DE genes: 900 genes

**Expression Distribution:**
- Base expression follows log-normal distribution
- Fold changes range from 2x to 5x for DE genes
- Counts follow negative binomial distribution (realistic for RNA-seq)

## Usage

### In Dashboard:
1. Upload `sample_counts.csv` as count matrix
2. Upload `sample_metadata.csv` as metadata
3. Select "condition" as the condition column
4. Choose comparison: treatment vs control
5. Run analysis

### In Python:
```python
import pandas as pd

# Load data
counts = pd.read_csv('examples/sample_counts.csv', index_col=0)
metadata = pd.read_csv('examples/sample_metadata.csv', index_col=0)

# Run analysis
from transcriptomics_dashboard.deseq2 import run_deseq2

results = run_deseq2(
    counts,
    metadata,
    condition_col='condition',
    contrast=('treatment', 'control')
)
```

## Generating New Datasets

You can generate custom datasets with different parameters:

```python
from generate_example_data import generate_example_data

counts, metadata, truth = generate_example_data(
    n_genes=2000,
    n_control=4,
    n_treatment=4,
    n_de_genes=200,
    fold_change_range=(2, 6),
    output_dir='my_data'
)
```

## Validation

The `ground_truth.csv` file contains the true DE status of each gene, useful for:
- Validating the analysis pipeline
- Calculating sensitivity and specificity
- Benchmarking different methods

## Notes

- This is **synthetic data** generated for demonstration purposes
- Expression patterns are simplified compared to real biological data
- No batch effects are included by default
- All samples have similar sequencing depths
