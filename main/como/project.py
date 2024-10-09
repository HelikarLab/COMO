from pathlib import Path


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
    def __init__(self, data_dir: Path = None, config_dir: Path = None, result_dir: Path = None, code_dir: Path = None) -> None:
        self.data_dir = data_dir
        if self.data_dir is None:
            self.data_dir = Path.cwd() / "data"
            self.data_dir.mkdir(exist_ok=True)

        self.config_dir = config_dir
        if self.config_dir is None:
            self.config_dir = Path.cwd() / "data" / "config_sheets"
            self.config_dir.mkdir(exist_ok=True)

        self.result_dir = result_dir
        if self.result_dir is None:
            self.result_dir = Path.cwd() / "data" / "results"
            self.result_dir.mkdir(exist_ok=True)

        self.code_dir = code_dir
        if self.code_dir is None:
            self.code_dir = Path.cwd() / "como"
