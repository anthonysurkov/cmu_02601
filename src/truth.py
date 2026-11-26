"""Author: Anthony Surkov
Nov 25 2026
Source-of-truth modules for machine learning pipe.
"""

from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import torch
import numpy as np
from sklearn.model_selection import KFold, train_test_split

from .features import get_protbert_dict
from .paths import REF_DIR
from .utils import generate_uid, mutation_to_seq


class BglBSubset(torch.utils.data.Dataset):
    """
    Stores BglB data as a Torch dataset, appropriate for use in Torch workflows.
    Yielded by BglBRepository on request.
    torch.utils.data.Dataset mandates an __init__, __len__, and __getitem__
    method.
    """
    def __init__(
        self,
        df: pd.DataFrame,
        protbert_dict: dict[str, torch.Tensor]
    ) -> None:
        self.df = df.copy()
        self.protbert_dict = protbert_dict

    def __len__(self) -> int:
        """Returns length of BglBSubset dataframe."""
        return len(self.df)

    def __getitem__(self, index) -> dict[str, torch.Tensor]:
        """Yields one dictionary of tensors for use with Torch."""
        row = self.df.iloc[index]
        return {
            "km": torch.tensor(
                row["km"],
                dtype=torch.float32
            ),
            "kcat": torch.tensor(
                row["kcat"],
                dtype=torch.float32
            ),
            "tm": torch.tensor(
                row["tm"],
                dtype=torch.float32
            ),
            "protbert": self.protbert_dict[row["seq"]]
        }

    def to_numpy(self, target: str = "tm") -> tuple[np.ndarray, np.ndarray]:
        """Written by a GPT model.
        Return mean-pooled protein features and the requested target.
        """
        X = np.stack([
            self.protbert_dict[seq]
            .detach()
            .mean(dim=0)
            .cpu()
            .numpy()
            for seq in self.df["seq"]
        ])
        y = self.df[target].to_numpy(dtype=np.float32)
        return X, y


class BglBRepository:
    """
    Owns canonicalized BglB data and features. Offers methods for normalization,
    data caching and loading, and train/val splitting (including K-folds).
    """
    def __init__(
        self,
        df: pd.DataFrame,
        seed: int = 42,
        persistent: bool = True,
        cache_path: Path = REF_DIR
    ) -> None:
        self.seed = seed

        self.cache_dir = cache_path
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        df = self._normalize(df)
        self.df = df.reset_index(drop=True)

        cached = self.load() if persistent else None
        if cached is None:
            self.protbert_dict = self.download_features()
            if persistent:
                self.save()
        else:
            self.protbert_dict: dict[str, torch.Tensor] = cached

    @staticmethod
    def _normalize(df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalizes BglB data: standardizes headers, handles NaN values,
        averages replicates, and returns relevant data.
        """
        assert "Variant" in df.keys()
        assert "KM (mM)" in df.keys() # original column headers from BglB DB
        assert "kcat (1/min)" in df.keys()
        assert "Tm (C)" in df.keys()

        df = df.copy()
        df["seq"] = df["Variant"].apply(mutation_to_seq)
        df = df.rename(columns={
            "KM (mM)": "km",
            "kcat (1/min)": "kcat",
            "Tm (C)": "tm"
        })
        df["tm"] = df.groupby("seq", sort=False)["tm"].transform("mean")
        df = df.dropna(subset=["km", "kcat", "tm"])
        df = df.groupby(
            "seq", as_index=False, sort=False
        )[["km", "kcat", "tm"]].mean()

        return df[["seq", "km", "kcat", "tm"]]

    def download_features(self) -> dict[str, torch.Tensor]:
        """Prompts features.py to download ProtBERT embeddings from Rostlab."""
        return get_protbert_dict(self.df["seq"])

    def save(self) -> None:
        """Saves downloaded ProtBERT embeddings to cache."""
        strings = pd.concat([
            self.df["seq"],
            self.df["km"],
            self.df["kcat"],
            self.df["tm"]
        ])
        uid = generate_uid(strings)
        torch.save(self.protbert_dict, self.cache_dir / f"{uid}_protbert.pt")

    def load(self) -> dict[str, torch.Tensor] | None:
        """Loads downloaded ProtBERT embeddings from cache."""
        strings = pd.concat([
            self.df["seq"],
            self.df["km"],
            self.df["kcat"],
            self.df["tm"]
        ])
        uid = generate_uid(strings)

        if not (self.cache_dir / f"{uid}_protbert.pt").exists():
            return None

        protbert_dict = torch.load(
            self.cache_dir / f"{uid}_protbert.pt",
            weights_only=True
        )
        return protbert_dict

    def get_train_val(self, size: float = 0.2) -> dict[str, BglBSubset]:
        """Provides one (1-x):x train:val split of the data as BglBSubsets."""
        train, val = train_test_split(
            self.df,
            test_size=size,
            random_state=self.seed
        )
        return {
            "train": BglBSubset(train, self.protbert_dict),
            "val": BglBSubset(val, self.protbert_dict)
        }

    def yield_train_val_folds(
        self,
        test_size: float = 0.20,
        n_splits: int = 5
    ) -> Iterator[dict[str, BglBSubset]]:
        """Provides n_split train-val folds of the data as BglBSubsets."""
        df, test = train_test_split(
            self.df,
            test_size=test_size,
            random_state=self.seed
        )
        folds = KFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=self.seed
        )
        for train_idx, val_idx in folds.split(df):
            yield {
                "train": BglBSubset(
                    df.iloc[train_idx],
                    self.protbert_dict
                ),
                "val": BglBSubset(
                    df.iloc[val_idx],
                    self.protbert_dict
                ),
                "test": BglBSubset(
                    test,
                    self.protbert_dict
                )
            }
