"""Data validation module for count matrices and metadata."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import numpy as np
from pydantic import BaseModel, Field, field_validator
from typing import Dict, List, Optional, Tuple, Union, Any


class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass


class ValidationWarning(BaseModel):
    """Warning message from validation."""
    message: str
    severity: str = Field(default="warning")  # warning, info


class ValidationResult(BaseModel):
    """Result of data validation."""
    valid: bool
    errors: List[str] = Field(default_factory=list)
    warnings: List[ValidationWarning] = Field(default_factory=list)
    summary: Dict[str, Any] = Field(default_factory=dict)


class CountMatrixSchema(BaseModel):
    """Schema for count matrix validation."""
    n_genes: int
    n_samples: int
    gene_ids: List[str]
    sample_ids: List[str]
    has_negative: bool
    has_non_integer: bool
    has_missing: bool
    library_sizes: Dict[str, float]
    
    
class MetadataSchema(BaseModel):
    """Schema for metadata validation."""
    n_samples: int
    sample_ids: List[str]
    columns: List[str]
    condition_column: Optional[str] = None
    n_conditions: Optional[int] = None
    replicates_per_condition: Optional[Dict[str, int]] = None


def read_count_matrix(
    filepath: Union[str, Path],
    delimiter: Optional[str] = None,
    gene_col: int = 0,
    header: int = 0
) -> pd.DataFrame:
    """
    Read count matrix from file.
    
    Args:
        filepath: Path to count matrix file
        delimiter: Column delimiter (auto-detected if None)
        gene_col: Column index for gene IDs
        header: Row index for column names
        
    Returns:
        DataFrame with genes as rows, samples as columns
    """
    filepath = Path(filepath)
    
    # Determine file type and read accordingly
    if filepath.suffix.lower() in ['.xlsx', '.xls']:
        df = pd.read_excel(filepath, index_col=gene_col, header=header)
    else:
        # Auto-detect delimiter if not provided
        if delimiter is None:
            with open(filepath, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    delimiter = '\t'
                elif ',' in first_line:
                    delimiter = ','
                else:
                    delimiter = None  # pandas will try to detect
        
        df = pd.read_csv(filepath, sep=delimiter, index_col=gene_col, header=header)
    
    # Clean up gene IDs (remove whitespace)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    
    return df


def read_metadata(
    filepath: Union[str, Path],
    delimiter: Optional[str] = None,
    sample_col: int = 0,
    header: int = 0
) -> pd.DataFrame:
    """
    Read metadata file.
    
    Args:
        filepath: Path to metadata file
        delimiter: Column delimiter (auto-detected if None)
        sample_col: Column index for sample IDs
        header: Row index for column names
        
    Returns:
        DataFrame with samples as rows, conditions/metadata as columns
    """
    filepath = Path(filepath)
    
    if filepath.suffix.lower() in ['.xlsx', '.xls']:
        df = pd.read_excel(filepath, index_col=sample_col, header=header)
    else:
        if delimiter is None:
            with open(filepath, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    delimiter = '\t'
                elif ',' in first_line:
                    delimiter = ','
                else:
                    delimiter = None
        
        df = pd.read_csv(filepath, sep=delimiter, index_col=sample_col, header=header)
    
    # Clean up sample IDs
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    
    return df


def validate_count_matrix(counts: pd.DataFrame) -> Tuple[ValidationResult, CountMatrixSchema]:
    """
    Validate count matrix.
    
    Args:
        counts: Count matrix DataFrame
        
    Returns:
        Tuple of (ValidationResult, CountMatrixSchema)
    """
    errors = []
    warnings = []
    
    # Basic structure checks
    if counts.empty:
        errors.append("Count matrix is empty")
        return ValidationResult(valid=False, errors=errors), None
    
    n_genes, n_samples = counts.shape
    
    # Check for negative values
    has_negative = (counts < 0).any().any()
    if has_negative:
        errors.append("Count matrix contains negative values")
    
    # Check for non-integer values
    has_non_integer = not np.allclose(counts.values, counts.values.astype(int), equal_nan=True)
    if has_non_integer:
        warnings.append(ValidationWarning(
            message="Count matrix contains non-integer values. They will be rounded.",
            severity="warning"
        ))
    
    # Check for missing values
    has_missing = counts.isna().any().any()
    if has_missing:
        n_missing = counts.isna().sum().sum()
        errors.append(f"Count matrix contains {n_missing} missing values")
    
    # Check gene count
    if n_genes < 5000:
        warnings.append(ValidationWarning(
            message=f"Low number of genes ({n_genes}). Typical RNA-seq has 15,000-25,000 genes.",
            severity="warning"
        ))
    
    # Calculate library sizes
    library_sizes = counts.sum(axis=0).to_dict()
    
    # Check library sizes
    for sample, size in library_sizes.items():
        if size < 1e6:
            warnings.append(ValidationWarning(
                message=f"Sample '{sample}' has low library size: {size:,.0f} reads",
                severity="warning"
            ))
        elif size > 100e6:
            warnings.append(ValidationWarning(
                message=f"Sample '{sample}' has very high library size: {size:,.0f} reads",
                severity="info"
            ))
    
    # Check for duplicate gene IDs
    if counts.index.duplicated().any():
        n_duplicates = counts.index.duplicated().sum()
        errors.append(f"Count matrix contains {n_duplicates} duplicate gene IDs")
    
    # Check for duplicate sample IDs
    if counts.columns.duplicated().any():
        n_duplicates = counts.columns.duplicated().sum()
        errors.append(f"Count matrix contains {n_duplicates} duplicate sample IDs")
    
    # Create schema
    schema = CountMatrixSchema(
        n_genes=n_genes,
        n_samples=n_samples,
        gene_ids=counts.index.tolist(),
        sample_ids=counts.columns.tolist(),
        has_negative=has_negative,
        has_non_integer=has_non_integer,
        has_missing=has_missing,
        library_sizes=library_sizes
    )
    
    # Summary
    summary = {
        "n_genes": n_genes,
        "n_samples": n_samples,
        "total_counts": int(counts.sum().sum()),
        "mean_library_size": float(np.mean(list(library_sizes.values()))),
        "median_library_size": float(np.median(list(library_sizes.values())))
    }
    
    result = ValidationResult(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        summary=summary
    )
    
    return result, schema


def validate_metadata(
    metadata: pd.DataFrame,
    count_samples: Optional[List[str]] = None,
    condition_column: Optional[str] = None
) -> Tuple[ValidationResult, MetadataSchema]:
    """
    Validate metadata.
    
    Args:
        metadata: Metadata DataFrame
        count_samples: List of sample IDs from count matrix (for matching check)
        condition_column: Name of condition column to validate
        
    Returns:
        Tuple of (ValidationResult, MetadataSchema)
    """
    errors = []
    warnings = []
    
    if metadata.empty:
        errors.append("Metadata is empty")
        return ValidationResult(valid=False, errors=errors), None
    
    n_samples = len(metadata)
    sample_ids = metadata.index.tolist()
    columns = metadata.columns.tolist()
    
    # Check if condition column exists
    replicates_per_condition = None
    n_conditions = None
    
    if condition_column:
        if condition_column not in metadata.columns:
            errors.append(f"Condition column '{condition_column}' not found in metadata")
        else:
            # Check for missing values in condition column
            if metadata[condition_column].isna().any():
                errors.append(f"Condition column '{condition_column}' contains missing values")
            
            # Count replicates per condition
            condition_counts = metadata[condition_column].value_counts()
            replicates_per_condition = condition_counts.to_dict()
            n_conditions = len(condition_counts)
            
            # Check for sufficient replicates
            for condition, count in replicates_per_condition.items():
                if count < 2:
                    errors.append(
                        f"Condition '{condition}' has only {count} replicate(s). "
                        "At least 2 replicates per condition are required."
                    )
                elif count < 3:
                    warnings.append(ValidationWarning(
                        message=f"Condition '{condition}' has only {count} replicates. "
                                "3+ replicates recommended for robust analysis.",
                        severity="warning"
                    ))
    
    # Check sample ID matching with count matrix
    if count_samples is not None:
        count_set = set(count_samples)
        meta_set = set(sample_ids)
        
        missing_in_meta = count_set - meta_set
        missing_in_counts = meta_set - count_set
        
        if missing_in_meta:
            errors.append(
                f"Samples in count matrix but not in metadata: {', '.join(sorted(missing_in_meta))}"
            )
        
        if missing_in_counts:
            warnings.append(ValidationWarning(
                message=f"Samples in metadata but not in count matrix: {', '.join(sorted(missing_in_counts))}",
                severity="info"
            ))
    
    # Check for duplicate sample IDs
    if metadata.index.duplicated().any():
        n_duplicates = metadata.index.duplicated().sum()
        errors.append(f"Metadata contains {n_duplicates} duplicate sample IDs")
    
    # Create schema
    schema = MetadataSchema(
        n_samples=n_samples,
        sample_ids=sample_ids,
        columns=columns,
        condition_column=condition_column,
        n_conditions=n_conditions,
        replicates_per_condition=replicates_per_condition
    )
    
    # Summary
    summary = {
        "n_samples": n_samples,
        "n_columns": len(columns),
        "columns": columns
    }
    
    if condition_column and replicates_per_condition:
        summary["conditions"] = replicates_per_condition
    
    result = ValidationResult(
        valid=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        summary=summary
    )
    
    return result, schema


def validate_analysis_inputs(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_column: str
) -> ValidationResult:
    """
    Validate complete analysis inputs.
    
    Args:
        counts: Count matrix
        metadata: Sample metadata
        condition_column: Column name for experimental conditions
        
    Returns:
        ValidationResult with combined validation from both inputs
    """
    all_errors = []
    all_warnings = []
    
    # Validate count matrix
    counts_result, counts_schema = validate_count_matrix(counts)
    all_errors.extend(counts_result.errors)
    all_warnings.extend(counts_result.warnings)
    
    # Validate metadata
    meta_result, meta_schema = validate_metadata(
        metadata,
        count_samples=counts.columns.tolist(),
        condition_column=condition_column
    )
    all_errors.extend(meta_result.errors)
    all_warnings.extend(meta_result.warnings)
    
    # Combined summary
    summary = {
        "counts": counts_result.summary,
        "metadata": meta_result.summary
    }
    
    return ValidationResult(
        valid=len(all_errors) == 0,
        errors=all_errors,
        warnings=all_warnings,
        summary=summary
    )


if __name__ == "__main__":
    # Example usage
    import numpy as np
    
    # Create sample data
    np.random.seed(42)
    counts = pd.DataFrame(
        np.random.poisson(100, (1000, 6)),
        index=[f"Gene_{i}" for i in range(1000)],
        columns=[f"Sample_{i}" for i in range(6)]
    )
    
    metadata = pd.DataFrame({
        'condition': ['control'] * 3 + ['treatment'] * 3,
        'batch': [1, 2, 1, 2, 1, 2]
    }, index=counts.columns)
    
    # Validate
    result = validate_analysis_inputs(counts, metadata, 'condition')
    
    print(f"Valid: {result.valid}")
    print(f"Errors: {result.errors}")
    print(f"Warnings: {len(result.warnings)}")
    print(f"Summary: {result.summary}")