from pathlib import Path
import yaml

with open("config.yaml") as f:
    cfg = yaml.safe_load(f)

# config.yaml root
DATA_DIR = Path(cfg["paths"]["data_dir"])
CACHE_DIR = Path(cfg["paths"]["cache_dir"])
