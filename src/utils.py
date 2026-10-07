import pandas as pd
import json
from uuid import uuid5, NAMESPACE_URL

from .paths import DATA_DIR


def generate_uid(strings: pd.Series) -> str:
    if strings.empty:
        raise ValueError("strings cannot be empty")
    name = json.dumps(strings.tolist())
    return str(uuid5(NAMESPACE_URL, name))
