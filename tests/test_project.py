from pathlib import Path

from como import project


def test_config():
    configs = project.Config()
    root_dir = Path(__file__).parent.parent

    assert configs.data_dir == root_dir / "data"
    assert configs.config_dir == root_dir / "data" / "config_sheets"
    assert configs.result_dir == root_dir / "data" / "results"
