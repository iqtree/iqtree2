"""
Set global variables to point to important paths
"""
from pathlib import Path

current_dir = Path(__file__).parent
repo_root = current_dir
while repo_root != Path("/") and not (repo_root / ".git").is_dir():
    repo_root = repo_root.parent

iqtree2_dir = repo_root / "build"
tests_root = repo_root / "tests"
tests_data = tests_root / "data"
iqtree1_dir = tests_root / "iqtree1"
