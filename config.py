from __future__ import annotations

import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


@dataclass(slots=True)
class RetryConfig:
    """Retry settings applied to outbound HTTP calls."""

    attempts: int = 3
    backoff_factor: float = 0.5
    status_forcelist: tuple[int, ...] = (500, 502, 503, 504)


@dataclass(slots=True)
class ServiceConfig:
    """Configuration structure for external orientation services."""

    domains: Optional[list[str]] = None
    base_url: Optional[str] = None
    case_sensitive: bool = False
    timeout: tuple[int, int] = (5, 15)
    retries: int = 3
    priority: int = 0
    prefer_ready_pdb: bool = False
    enabled: bool = True


@dataclass(slots=True)
class PipelineConfig:
    """Derived paths and concurrency settings for the pipeline execution."""

    base_dir: Path
    input_template: Path
    output_excel: Path
    parquet_dir: Path
    exp_dir: Path
    pred_dir: Path
    aligned_dir: Path
    oriented_dir: Path
    log_file: Path
    window: int = 5
    max_workers: int = max(1, multiprocessing.cpu_count() or 1)
    api_retry: RetryConfig = field(default_factory=RetryConfig)
    prefer_native_sources: bool = True
    require_native_reference: bool = False
    native_source_priority: tuple[str, ...] = ("OPM", "PDBTM")
    # Feature flags
    # Use XlsxWriter for fast Excel export (default: True)
    use_xlsxwriter: bool = True
    # Enable PDBFixer-based repair/protonation (default: False)
    use_repair: bool = True
    # Attempt to install repair dependencies automatically
    auto_install_repair: bool = True
    # Backwards-compatible alias for structure repair (default: False)
    use_pdb_fixer: bool = True
    # Enable orientation fallback cascade (default: False)
    use_orientation_fallback: bool = True
    # Activate enhanced multi-source orientation
    enable_new_orientation_services: bool = True
    prefer_pdbtm_over_memprotmd: bool = False  # Legacy ordering preference toggle
    # Path to local PPM executable for orientation fallback
    ppm_path: Optional[str] = None
    # Use HTTP/2 for API/networking if available (default: False)
    use_http2: bool = False

    def __post_init__(self) -> None:
        """Keep legacy repair flags in sync."""
        if self.use_repair and not self.use_pdb_fixer:
            self.use_pdb_fixer = True
        elif self.use_pdb_fixer and not self.use_repair:
            self.use_repair = True
        if self.auto_install_repair and not self.use_repair:
            self.auto_install_repair = False

    @classmethod
    def from_args(
        cls,
        base: Optional[str] = None,
        input_template: Optional[str] = None,
        output_excel: Optional[str] = None,
        window: int = 5,
    ) -> "PipelineConfig":
        """Build a configuration instance from CLI-style inputs.

        Args:
            base (Optional[str]): Base working directory for all artefacts.
            input_template (Optional[str]): Path to the UniProt import template.
            output_excel (Optional[str]): Destination Excel filename.
            window (int): Sliding window size for alignment operations.

        Returns:
            PipelineConfig: Fully initialised configuration object.
        """
        base_path = Path(base).expanduser().resolve(
        ) if base else Path.cwd().resolve()
        input_path = (
            Path(input_template).expanduser().resolve()
            if input_template
            else base_path / "uniprot_template.xlsx"
        )
        output_path = (
            Path(output_excel).expanduser().resolve()
            if output_excel
            else base_path / "uniprot_results.xlsx"
        )
        cache_dir = base_path / "cache"
        parquet_dir = cache_dir / "parquet"
        return cls(
            base_dir=base_path,
            input_template=input_path,
            output_excel=output_path,
            parquet_dir=parquet_dir,
            exp_dir=base_path / "pdb_structures_exp",
            pred_dir=base_path / "pdb_structures_pred",
            aligned_dir=base_path / "pdb_structures_af_aligned",
            oriented_dir=base_path / "pdb_structures_oriented",
            log_file=base_path / "pipeline.log",
            window=window,
        )

    def ensure_directories(self) -> None:
        """Create required output directories if they are missing."""
        self.output_excel.parent.mkdir(parents=True, exist_ok=True)
        for directory in (
            self.base_dir,
            self.parquet_dir,
            self.exp_dir,
            self.pred_dir,
            self.aligned_dir,
            self.oriented_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)


SERVICE_CONFIG: Dict[str, ServiceConfig] = {
    "pdbtm": ServiceConfig(priority=1, enabled=True),
    "ppm": ServiceConfig(priority=2, enabled=True),
    "memembed": ServiceConfig(priority=3, enabled=True),
    "tmdet": ServiceConfig(priority=4, enabled=True),
    "axes": ServiceConfig(priority=5, enabled=True),
}

__all__ = [
    "PipelineConfig",
    "RetryConfig",
    "ServiceConfig",
    "SERVICE_CONFIG",
]
