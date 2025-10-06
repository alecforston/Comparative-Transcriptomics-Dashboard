"""Generate example datasets for testing the transcriptomics dashboard."""

import numpy as np
import pandas as pd
from pathlib import Path


def generate_example_data(
    n_genes: int = 1000,
    n_control: int = 3,
    n_treatment: int = 3,
    n_de_genes: int = 100,
    fold_change_range: tuple = (2, 5),
    output_dir: str = "examples",
    seed: int = 42
):
    """
    Generate synthetic RNA-seq count data with differential expression.
    
    Args:
        n_genes: Total number of genes
        n_control: Number of control samples
        n_treatment: Number of treatment samples
        n_de_genes: Number of differentially expressed genes
        fold_change_range: (min, max) fold change for DE genes
        output_dir: Directory to save files
        seed: Random seed for reproducibility
    """ 
    np.random.seed(seed)
    
    n_samples = n_control + n_treatment
    
    # Generate gene names
    gene_names = [f"Gene_{i:05d}" for i in range(n_genes)]
    sample_names = (
        [f"Control_{i+1}" for i in range(n_control)] +
        [f"Treatment_{i+1}" for i in range(n_treatment)]
    )
    
    # Generate base expression levels (log-normal distribution)
    base_expression = np.random.lognormal(mean=5, sigma=2, size=n_genes)
    
    # Generate count matrix
    counts = np.zeros((n_genes, n_samples))
    
    # Generate control samples
    for i in range(n_control):
        # Negative binomial with varying dispersion
        dispersion = np.random.uniform(0.05, 0.2, n_genes)
        counts[:, i] = np.random.negative_binomial(
            n=1/dispersion,
            p=1/(1 + base_expression * dispersion)
        )
    
    # Generate treatment samples with DE
    de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
    
    # Split DE genes into up and down regulated
    n_up = n_de_genes // 2
    up_indices = de_indices[:n_up]
    down_indices = de_indices[n_up:]
    
    # Generate fold changes
    fc_up = np.random.uniform(fold_change_range[0], fold_change_range[1], n_up)
    fc_down = np.random.uniform(1/fold_change_range[1], 1/fold_change_range[0], len(down_indices))
    
    treatment_expression = base_expression.copy()
    treatment_expression[up_indices] *= fc_up
    treatment_expression[down_indices] *= fc_down
    
    for i in range(n_treatment):
        dispersion = np.random.uniform(0.05, 0.2, n_genes)
        counts[:, n_control + i] = np.random.negative_binomial(
            n=1/dispersion,
            p=1/(1 + treatment_expression * dispersion)
        )
    
    # Create count matrix DataFrame
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=gene_names,
        columns=sample_names
    )
    
    # Create metadata DataFrame
    metadata_df = pd.DataFrame({
        'condition': ['control'] * n_control + ['treatment'] * n_treatment,
        'batch': [1, 2, 1, 2, 1, 2][:n_samples],
        'sequencing_depth': counts_df.sum(axis=0).values,
        'replicate': list(range(1, n_control + 1)) + list(range(1, n_treatment + 1))
    }, index=sample_names)
    
    # Create ground truth DataFrame (for validation)
    ground_truth = pd.DataFrame({
        'gene': gene_names,
        'is_de': False,
        'true_fc': 1.0,
        'direction': 'none'
    })
    
    ground_truth.loc[up_indices, 'is_de'] = True
    ground_truth.loc[up_indices, 'true_fc'] = fc_up
    ground_truth.loc[up_indices, 'direction'] = 'up'
    
    ground_truth.loc[down_indices, 'is_de'] = True
    ground_truth.loc[down_indices, 'true_fc'] = fc_down
    ground_truth.loc[down_indices, 'direction'] = 'down'
    
    # Save files
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    counts_df.to_csv(output_path / "sample_counts.csv")
    metadata_df.to_csv(output_path / "sample_metadata.csv")
    ground_truth.to_csv(output_path / "ground_truth.csv", index=False)
    
    # Also save as TSV and Excel
    counts_df.to_csv(output_path / "sample_counts.tsv", sep='\t')
    metadata_df.to_csv(output_path / "sample_metadata.tsv", sep='\t')
    
    counts_df.to_excel(output_path / "sample_counts.xlsx")
    metadata_df.to_excel(output_path / "sample_metadata.xlsx")
    
    print(f"✓ Generated example data:")
    print(f"  - Genes: {n_genes}")
    print(f"  - Samples: {n_samples} ({n_control} control, {n_treatment} treatment)")
    print(f"  - DE genes: {n_de_genes} ({n_up} up, {len(down_indices)} down)")
    print(f"  - Files saved to: {output_path.absolute()}")
    
    return counts_df, metadata_df, ground_truth


def generate_large_dataset(output_dir: str = "examples/large"):
    """Generate a larger, more realistic dataset."""
    return generate_example_data(
        n_genes=20000,
        n_control=4,
        n_treatment=4,
        n_de_genes=2000,
        fold_change_range=(1.5, 8),
        output_dir=output_dir,
        seed=42
    )


def generate_minimal_dataset(output_dir: str = "examples/minimal"):
    """Generate a minimal dataset for quick testing."""
    return generate_example_data(
        n_genes=100,
        n_control=2,
        n_treatment=2,
        n_de_genes=10,
        fold_change_range=(3, 5),
        output_dir=output_dir,
        seed=42
    )


def create_readme(output_dir: str = "examples"):
    """Create README for example data."""
    readme_content = """# Example Data

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
"""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    with open(output_path / "README.md", 'w') as f:
        f.write(readme_content)
    
    print(f"✓ Created README in {output_path}")


if __name__ == "__main__":
    # Generate standard example dataset
    print("Generating standard example dataset...")
    generate_example_data()
    
    # Generate minimal dataset for quick tests
    print("\nGenerating minimal dataset...")
    generate_minimal_dataset()
    
    # Create README
    print("\nCreating README...")
    create_readme()
    
    print("\n✓ All example data generated successfully!")