from pathlib import Path

from loguru import logger


class SingletonMeta(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        """
        Possible changes to the value of the `__init__` argument do not affect
        the returned instance.
        """
        if cls not in cls._instances:
            instance = super().__call__(*args, **kwargs)
            cls._instances[cls] = instance
        return cls._instances[cls]


class Config(metaclass=SingletonMeta):
    def __init__(
        self,
        data_dir: Path = None,
        config_dir: Path = None,
        result_dir: Path = None,
    ) -> None:
        current_dir = Path.cwd()

        self.data_dir = data_dir
        if self.data_dir is None:
            logger.warning(f"'data_dir' not provided to Config, using {Path.cwd() / 'data'}")
            self.data_dir = current_dir / "data"
            self.data_dir.mkdir(parents=True, exist_ok=True)

        self.config_dir = config_dir
        if self.config_dir is None:
            logger.warning(f"'config_dir' not provided to Config, using {self.data_dir / 'config_sheets'}")
            self.config_dir = self.data_dir / "config_sheets"
            self.config_dir.mkdir(parents=True, exist_ok=True)

        self.result_dir = result_dir
        if self.result_dir is None:
            logger.warning(f"'results_dir' not provided to Config, using {self.data_dir / 'results'}")
            self.result_dir = self.data_dir / "results"
            self.result_dir.mkdir(parents=True, exist_ok=True)
