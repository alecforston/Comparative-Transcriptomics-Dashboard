"""Configuration management for transcriptomics dashboard."""

import os
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field
from pydantic_settings import BaseSettings
import yaml


class AnalysisDefaults(BaseModel):
    """Default analysis parameters."""
    
    fdr_threshold: float = Field(default=0.05, ge=0.0, le=1.0)
    log2fc_threshold: float = Field(default=1.0, ge=0.0)
    min_count: int = Field(default=10, ge=0)
    threads: int = Field(default=-1)  # -1 means use all available
    
    
class PathConfig(BaseModel):
    """Path configurations."""
    
    user_home: Path = Field(default_factory=lambda: Path.home() / ".transcriptomics_dashboard")
    data_dir: Optional[Path] = None
    reference_dir: Optional[Path] = None
    sessions_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None
    
    def __init__(self, **data):
        super().__init__(**data)
        # Set derived paths if not provided
        if self.data_dir is None:
            self.data_dir = self.user_home / "data"
        if self.reference_dir is None:
            self.reference_dir = self.user_home / "data" / "reference"
        if self.sessions_dir is None:
            self.sessions_dir = self.user_home / "sessions"
        if self.cache_dir is None:
            self.cache_dir = self.user_home / "cache"
            
    def create_directories(self):
        """Create all necessary directories."""
        for path in [self.user_home, self.data_dir, self.reference_dir, 
                     self.sessions_dir, self.cache_dir]:
            path.mkdir(parents=True, exist_ok=True)


class PerformanceConfig(BaseModel):
    """Performance-related settings."""
    
    max_cache_size_gb: float = Field(default=10.0, ge=0.1)
    chunk_size: int = Field(default=10000, ge=1000)
    parallel_backend: str = Field(default="loky")
    memory_warning_threshold: float = Field(default=0.8, ge=0.1, le=1.0)


class ReferenceDatabaseConfig(BaseModel):
    """Reference database settings."""
    
    go_obo_url: str = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    update_interval_days: int = Field(default=180, ge=1)
    default_organism: str = "human"
    supported_organisms: list[str] = ["human", "mouse", "rat"]


class Config(BaseSettings):
    """Main configuration class."""
    
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    paths: PathConfig = Field(default_factory=PathConfig)
    performance: PerformanceConfig = Field(default_factory=PerformanceConfig)
    reference: ReferenceDatabaseConfig = Field(default_factory=ReferenceDatabaseConfig)
    
    # App settings
    app_title: str = "Transcriptomics Dashboard"
    app_version: str = "0.1.0"
    debug: bool = False
    host: str = "127.0.0.1"
    port: int = 8050
    
    class Config:
        env_prefix = "TRANS_"
        env_nested_delimiter = "__"
    
    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        """Load configuration from YAML file."""
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    def to_yaml(self, path: Path):
        """Save configuration to YAML file."""
        data = self.model_dump()
        # Convert Path objects to strings
        def convert_paths(obj):
            if isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            elif isinstance(obj, Path):
                return str(obj)
            return obj
        
        data = convert_paths(data)
        
        with open(path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    
    def initialize(self):
        """Initialize the configuration (create directories, etc.)."""
        self.paths.create_directories()
        
        # Create default config file if it doesn't exist
        config_file = self.paths.user_home / "config.yaml"
        if not config_file.exists():
            self.to_yaml(config_file)


# Global configuration instance
_config: Optional[Config] = None


def get_config() -> Config:
    """Get the global configuration instance."""
    global _config
    if _config is None:
        # Try to load from default location
        default_config_path = Path.home() / ".transcriptomics_dashboard" / "config.yaml"
        if default_config_path.exists():
            _config = Config.from_yaml(default_config_path)
        else:
            _config = Config()
            _config.initialize()
    return _config


def set_config(config: Config):
    """Set the global configuration instance."""
    global _config
    _config = config


# Example config.yaml template
CONFIG_TEMPLATE = """
# Transcriptomics Dashboard Configuration
# This file controls default settings for analysis and application behavior

defaults:
  fdr_threshold: 0.05        # False discovery rate threshold
  log2fc_threshold: 1.0      # Log2 fold change threshold
  min_count: 10              # Minimum read count filter
  threads: -1                # Number of threads (-1 = use all available)

paths:
  user_home: ~/.transcriptomics_dashboard
  # Other paths are auto-generated under user_home if not specified
  # data_dir: ~/.transcriptomics_dashboard/data
  # reference_dir: ~/.transcriptomics_dashboard/data/reference
  # sessions_dir: ~/.transcriptomics_dashboard/sessions
  # cache_dir: ~/.transcriptomics_dashboard/cache

performance:
  max_cache_size_gb: 10.0              # Maximum cache size in GB
  chunk_size: 10000                    # Chunk size for large file operations
  parallel_backend: loky               # Backend for parallel processing
  memory_warning_threshold: 0.8        # Warn when memory usage exceeds this fraction

reference:
  go_obo_url: http://purl.obolibrary.org/obo/go/go-basic.obo
  update_interval_days: 180            # Days before warning about old reference data
  default_organism: human
  supported_organisms:
    - human
    - mouse
    - rat

# Application settings
app_title: Transcriptomics Dashboard
app_version: 0.1.0
debug: false
host: 127.0.0.1
port: 8050
"""


if __name__ == "__main__":
    # Example usage
    config = get_config()
    print(f"User home: {config.paths.user_home}")
    print(f"FDR threshold: {config.defaults.fdr_threshold}")
    print(f"Threads: {config.defaults.threads}")