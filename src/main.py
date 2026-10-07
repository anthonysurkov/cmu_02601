import pandas as pd
from pathlib import Path

from .paths import DATA_DIR


# Intended use from project directory: python -m src.main
def main():
    infile = Path(DATA_DIR / "bglb.csv")
    df = pd.read_csv(infile)
    print(df)

if __name__ == "__main__":
    main()
