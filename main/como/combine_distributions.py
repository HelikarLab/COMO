from pathlib import Path
from typing import Literal, Union


def parse_contexts_zfpkm(wd: Union[str, Path], contexts: list[str], prep: Literal["mrna", "total"]):
    wd: Path = Path(wd)

    batches = []
    for context in contexts:
        dir_name = Path(wd, context, prep)
        files = dir_name.glob(f"zFPKM_Matrix_{prep}_*.csv")
        batches += [Path(file).stem for file in files]

    return batches


if __name__ == "__main__":
    result = parse_contexts_zfpkm(wd="/Users/joshl/PycharmProjects/COMO/main/data/results", contexts=["naiveB"], prep="total")
    print(result)
