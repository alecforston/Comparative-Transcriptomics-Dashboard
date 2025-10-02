"""Unit tests for validation module."""

import pytest
import pandas as pd
import numpy as np
from transcriptomics_dashboard.validation import (
    validate_count_matrix,
    validate_metadata,
    validate_analysis_inputs,
    read_count_matrix,
    ValidationError
)


@pytest.fixture
def valid_counts():
    """Create valid count matrix."""
    np.random.seed(42)
    return pd.DataFrame(
        np.random.poisson(100, (1000, 6)),
        index=[f"Gene_{i}" for i in range(1000)],
        columns=[f"Sample_{i}" for i in range(6)]
    )


@pytest.fixture
def valid_metadata():
    """Create valid metadata."""
    return pd.DataFrame({
        'condition': ['control'] * 3 + ['treatment'] * 3,
        'batch': [1, 2, 1, 2, 1, 2]
    }, index=[f"Sample_{i}" for i in range(6)])


class TestCountMatrixValidation:
    """Tests for count matrix validation."""
    
    def test_valid_count_matrix(self, valid_counts):
        """Test validation of valid count matrix."""
        result, schema = validate_count_matrix(valid_counts)
        
        assert result.valid
        assert len(result.errors) == 0
        assert schema.n_genes == 1000
        assert schema.n_samples == 6
        assert not schema.has_negative
        assert not schema.has_missing
    
    def test_negative_counts(self, valid_counts):
        """Test detection of negative counts."""
        invalid_counts = valid_counts.copy()
        invalid_counts.iloc[0, 0] = -5
        
        result, schema = validate_count_matrix(invalid_counts)
        
        assert not result.valid
        assert any('negative' in err.lower() for err in result.errors)
        assert schema.has_negative
    
    def test_missing_values(self, valid_counts):
        """Test detection of missing values."""
        invalid_counts = valid_counts.copy()
        invalid_counts.iloc[0, 0] = np.nan
        
        result, schema = validate_count_matrix(invalid_counts)
        
        assert not result.valid
        assert any('missing' in err.lower() for err in result.errors)
        assert schema.has_missing
    
    def test_non_integer_values(self, valid_counts):
        """Test detection of non-integer values."""
        invalid_counts = valid_counts.copy().astype(float)
        invalid_counts.iloc[0, 0] = 10.5
        
        result, schema = validate_count_matrix(invalid_counts)
        
        assert result.valid  # Should still be valid but with warning
        assert len(result.warnings) > 0
        assert any('non-integer' in w.message.lower() for w in result.warnings)
    
    def test_duplicate_gene_ids(self, valid_counts):
        """Test detection of duplicate gene IDs."""
        invalid_counts = valid_counts.copy()
        invalid_counts.index = ['Gene_0'] * len(invalid_counts)
        
        result, schema = validate_count_matrix(invalid_counts)
        
        assert not result.valid
        assert any('duplicate' in err.lower() for err in result.errors)
    
    def test_low_gene_count(self):
        """Test warning for low number of genes."""
        small_counts = pd.DataFrame(
            np.random.poisson(100, (100, 6)),
            index=[f"Gene_{i}" for i in range(100)],
            columns=[f"Sample_{i}" for i in range(6)]
        )
        
        result, schema = validate_count_matrix(small_counts)
        
        assert result.valid
        assert len(result.warnings) > 0
        assert any('low number of genes' in w.message.lower() for w in result.warnings)
    
    def test_library_size_warnings(self, valid_counts):
        """Test library size validation."""
        # Make one sample have very low counts
        low_counts = valid_counts.copy()
        low_counts.iloc[:, 0] = 1  # Very low library size
        
        result, schema = validate_count_matrix(low_counts)
        
        assert result.valid
        assert len(result.warnings) > 0
        assert any('low library size' in w.message.lower() for w in result.warnings)


class TestMetadataValidation:
    """Tests for metadata validation."""
    
    def test_valid_metadata(self, valid_metadata):
        """Test validation of valid metadata."""
        result, schema = validate_metadata(valid_metadata, condition_column='condition')
        
        assert result.valid
        assert len(result.errors) == 0
        assert schema.n_samples == 6
        assert schema.n_conditions == 2
        assert schema.replicates_per_condition['control'] == 3
        assert schema.replicates_per_condition['treatment'] == 3
    
    def test_missing_condition_column(self, valid_metadata):
        """Test error when condition column doesn't exist."""
        result, schema = validate_metadata(valid_metadata, condition_column='nonexistent')
        
        assert not result.valid
        assert any('not found' in err.lower() for err in result.errors)
    
    def test_insufficient_replicates(self):
        """Test error for insufficient replicates."""
        metadata = pd.DataFrame({
            'condition': ['control', 'treatment']
        }, index=['Sample_0', 'Sample_1'])
        
        result, schema = validate_metadata(metadata, condition_column='condition')
        
        assert not result.valid
        assert any('replicate' in err.lower() for err in result.errors)
    
    def test_sample_matching(self, valid_metadata, valid_counts):
        """Test sample name matching between counts and metadata."""
        # Metadata with different samples
        wrong_metadata = valid_metadata.copy()
        wrong_metadata.index = ['Wrong_' + str(i) for i in range(6)]
        
        result, schema = validate_metadata(
            wrong_metadata,
            count_samples=valid_counts.columns.tolist(),
            condition_column='condition'
        )
        
        assert not result.valid
        assert any('not in metadata' in err.lower() for err in result.errors)
    
    def test_missing_values_in_condition(self, valid_metadata):
        """Test error for missing values in condition column."""
        invalid_metadata = valid_metadata.copy()
        invalid_metadata.loc['Sample_0', 'condition'] = np.nan
        
        result, schema = validate_metadata(invalid_metadata, condition_column='condition')
        
        assert not result.valid
        assert any('missing values' in err.lower() for err in result.errors)


class TestAnalysisInputValidation:
    """Tests for complete analysis input validation."""
    
    def test_valid_complete_inputs(self, valid_counts, valid_metadata):
        """Test validation of complete valid inputs."""
        result = validate_analysis_inputs(valid_counts, valid_metadata, 'condition')
        
        assert result.valid
        assert len(result.errors) == 0
        assert 'counts' in result.summary
        assert 'metadata' in result.summary
    
    def test_combined_errors(self, valid_counts, valid_metadata):
        """Test that errors from both validations are combined."""
        # Create invalid counts
        invalid_counts = valid_counts.copy()
        invalid_counts.iloc[0, 0] = -5
        
        # Create invalid metadata
        invalid_metadata = valid_metadata.copy()
        invalid_metadata.index = ['Wrong_' + str(i) for i in range(6)]
        
        result = validate_analysis_inputs(invalid_counts, invalid_metadata, 'condition')
        
        assert not result.valid
        assert len(result.errors) >= 2  # At least one from each


class TestFileReading:
    """Tests for file reading functions."""
    
    def test_read_csv(self, tmp_path, valid_counts):
        """Test reading CSV file."""
        filepath = tmp_path / "counts.csv"
        valid_counts.to_csv(filepath)
        
        df = read_count_matrix(filepath)
        
        assert df.shape == valid_counts.shape
        assert list(df.columns) == list(valid_counts.columns)
    
    def test_read_tsv(self, tmp_path, valid_counts):
        """Test reading TSV file."""
        filepath = tmp_path / "counts.tsv"
        valid_counts.to_csv(filepath, sep='\t')
        
        df = read_count_matrix(filepath)
        
        assert df.shape == valid_counts.shape
    
    def test_read_excel(self, tmp_path, valid_counts):
        """Test reading Excel file."""
        filepath = tmp_path / "counts.xlsx"
        valid_counts.to_excel(filepath)
        
        df = read_count_matrix(filepath)
        
        assert df.shape == valid_counts.shape


if __name__ == "__main__":
    pytest.main([__file__, "-v"])