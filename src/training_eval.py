"""
(8) XGBoost training and evaluation pipeline

Takes:
- list of (X_train, y_train, X_val, y_val) tuples from pre-split k-fold CV
- optional (X_test, y_test) arrays for final held-out evaluation
Outputs:
    - tuple containing the list of fold models, mean validation RMSE,
    and test RMSE
"""


import numpy as np
from sklearn.metrics import root_mean_squared_error
from xgboost import XGBRegressor


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


def evaluate_kfold(
    folds: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    X_test: np.ndarray | None = None,
    y_test: np.ndarray | None = None,
) -> tuple[list[XGBRegressor], float, float | None]:
    """
    Runs XGBoost training across pre-split k-fold train/val sets.

    This function iterates through each fold provided by the data loader,
    collects the trained models and validation RMSEs, and averages fold
    predictions if a held-out test set is provided.
    """
    models: list[XGBRegressor] = []
    fold_rmses: list[float] = []

    # Training and validating across each pre-split fold
    for X_train, y_train, X_val, y_val in folds:
        model, val_rmse = _train_single_fold(X_train, y_train, X_val, y_val)
        models.append(model)
        fold_rmses.append(val_rmse)

    mean_val_rmse = float(np.mean(fold_rmses))

    # Evalutation of average on held-out test set
    test_rmse: float | None = None
    if X_test is not None and y_test is not None:
        test_preds = np.mean([m.predict(X_test) for m in models], axis=0)
        test_rmse = float(root_mean_squared_error(y_test, test_preds))

    return models, mean_val_rmse, test_rmse


if __name__ == "__main__":
    # Testing
    dummy_folds = []
    for fold in range(5):
        X_train = np.random.normal(size=(80, 16))
        y_train = np.random.normal(size=(80,))
        X_val = np.random.normal(size=(20, 16))
        y_val = np.random.normal(size=(20,))
        dummy_folds.append((X_train, y_train, X_val, y_val))

    dummy_X_test = np.random.normal(size=(20, 16))
    dummy_y_test = np.random.normal(size=(20,))

    models, val_rmse, test_rmse = evaluate_kfold(
        dummy_folds, dummy_X_test, dummy_y_test
    )
    print(f"Mean 5-Fold Validation RMSE: {val_rmse:.4f}")
    print(f"Test RMSE: {test_rmse:.4f}")