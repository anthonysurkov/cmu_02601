"""
(8) XGBoost training and evaluation pipeline

Takes:
- iterable of train/val/test BglBSubset dictionaries from pre-split k-fold CV
Outputs:
    - tuple containing the list of fold models, mean validation RMSE,
    and test RMSE
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
from sklearn.metrics import root_mean_squared_error
from xgboost import XGBRegressor

if TYPE_CHECKING:
    from .truth import BglBSubset


RANDOM_STATE = 42


def _train_single_fold(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
) -> tuple[XGBRegressor, float]:
    """
    Trains an XGBRegressor on one fold and computes its validation RMSE.

    This helper function fits a default XGBoost regressor on the training
    split and predicts on the validation split. Then finally returns the model 
    and RMSE.
    """
    model = XGBRegressor(random_state=RANDOM_STATE)
    model.fit(X_train, y_train)

    val_predictions = model.predict(X_val)
    val_rmse = float(root_mean_squared_error(y_val, val_predictions))

    return model, val_rmse

def infer_tm50(
    models: list[XGBRegressor],
    X: np.ndarray
) -> np.ndarray:
    """Predict one value per sample via fold ensemble average."""
    # some raise ValueErrors would be helpful in this function
    X = np.asarray(X)
    return np.stack(
        [model.predict(X) for model in models], axis=0
    ).mean(axis=0)

def evaluate_kfold(
    folds: Iterable[dict[str, BglBSubset]],
) -> tuple[list[XGBRegressor], float, float]:
    """
    Runs XGBoost training across pre-split k-fold train/val sets.
    """
    models: list[XGBRegressor] = []
    fold_rmses: list[float] = []
    X_test: np.ndarray | None = None
    y_test: np.ndarray | None = None

    # Training and validating across each pre-split fold:
    for fold in folds:
        if X_test is None:
            X_test, y_test = fold["test"].to_numpy()
        X_train, y_train = fold["train"].to_numpy()
        X_val, y_val = fold["val"].to_numpy()

        model, val_rmse = _train_single_fold(X_train, y_train, X_val, y_val)

        models.append(model)
        fold_rmses.append(val_rmse)

    if not models:
        raise ValueError("At least one fold is required.")

    mean_val_rmse = float(np.mean(fold_rmses))

    # Evaluate the ensemble once on the shared held-out test set:
    test_preds = infer_tm50(models, X_test)
    test_rmse = float(root_mean_squared_error(y_test, test_preds))

    return models, mean_val_rmse, test_rmse


if __name__ == "__main__":
    import pandas as pd
    import torch

    from .truth import BglBSubset

    def dummy_subset(size: int) -> BglBSubset:
        seqs = [str(i) for i in range(size)]
        df = pd.DataFrame({"seq": seqs, "tm": np.random.normal(size=size)})
        embeddings = {seq: torch.randn(1, 16) for seq in seqs}
        return BglBSubset(df, embeddings)

    # Testing
    dummy_folds = []
    dummy_test = dummy_subset(20)
    for fold in range(5):
        dummy_folds.append({
            "train": dummy_subset(80),
            "val": dummy_subset(20),
            "test": dummy_test,
        })

    models, val_rmse, test_rmse = evaluate_kfold(
        dummy_folds
    )
    print(f"Mean 5-Fold Validation RMSE: {val_rmse:.4f}")
    print(f"Test RMSE: {test_rmse:.4f}")
