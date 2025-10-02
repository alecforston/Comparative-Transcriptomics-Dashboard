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