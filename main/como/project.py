from __future__ import annotations

from pathlib import Path
from typing import ClassVar

from loguru import logger


class SingletonMeta(type):
    _instances: ClassVar[dict] = {}

    def __call__(cls, *args, **kwargs):
        """Validate that changes to the `__init__` argument do not affect the returned instance."""
        if cls not in cls._instances:
            instance = super().__call__(*args, **kwargs)
            cls._instances[cls] = instance
        return cls._instances[cls]


class Config(metaclass=SingletonMeta):
    def __init__(
        self,
        data_dir: Path | None = None,
        config_dir: Path | None = None,
        result_dir: Path | None = None,
    ) -> None:
        """Initialize the Config object."""
        current_dir = Path.cwd()

        self.data_dir = Path(data_dir) if data_dir else None
        if self.data_dir is None:
            logger.warning(f"'data_dir' not provided to Config, using {Path.cwd() / 'data'}")
            self.data_dir = current_dir / "data"
            self.data_dir.mkdir(parents=True, exist_ok=True)

        self.config_dir = Path(config_dir) if config_dir else None
        if self.config_dir is None:
            logger.warning(f"'config_dir' not provided to Config, using {self.data_dir / 'config_sheets'}")
            self.config_dir = self.data_dir / "config_sheets"
            self.config_dir.mkdir(parents=True, exist_ok=True)

        self.result_dir = Path(result_dir) if result_dir else None
        if self.result_dir is None:
            logger.warning(f"'results_dir' not provided to Config, using {self.data_dir / 'results'}")
            self.result_dir = self.data_dir / "results"
            self.result_dir.mkdir(parents=True, exist_ok=True)

        # Additional directories
        self.code_dir = current_dir / "main" / "como"
        self.log_dir = self.data_dir / "logs"
        self.matrix_dir = self.data_dir / "data_matrices"
        self.figures_dir = self.result_dir / "figures"

        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.matrix_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

    def update(self, **kwargs):
        """Update a key in the config object."""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, Path(value) if value else getattr(self, key))
            else:
                logger.warning(f"{key} is not a valid attribute of Config")

    def get_context_path(self, context_name: str, create: bool = True) -> Path:
        """Get path for a specific context, optionally creating it."""
        path = self.result_dir / context_name
        if create:
            path.mkdir(parents=True, exist_ok=True)
        return path

    def get_r_path(self, path: Path) -> str:
        """Convert a Path object to an R-compatible path string."""
        return path.as_posix()

    def get_matrix_path(self, context_name: str, filename: str) -> Path:
        """Get path for a matrix file in a specific context."""
        path = self.matrix_dir / context_name
        path.mkdir(parents=True, exist_ok=True)
        return path / filename
