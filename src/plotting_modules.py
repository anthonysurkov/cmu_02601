"""
This script creates plots detailing the quality of each of the models tested in
part 8.

Inputs:
1. CSVs with the true y values (y_val) and the predicted y values (y_pred) in
   the first and second columns, respectively, and no headers. There should be
   one CSV for each model evaluated. The CSVs are expected to be labeled
   according to the type of model being tested.
2. An output folder path where all the plots will be outputted.

Outputs:
1. A scatterplot for each model
2. A residual scatterplot for each model
3. One MSE bar graph with MSEs for all the models.
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from pathlib import Path

scatter_path = ""
residual_path = ""
mse_path = ""

def parse_csv(fn: str) -> (np.ndarray, np.ndarray):
    """
    Loads the CSVs into two numpy arrays do be used for downstream analysis
    """
    df = pd.read_csv(fn, header=None) # Will change header if needed later
    y_val = df.iloc[:, 0].to_numpy()
    y_pred = df.iloc[:, 1].to_numpy()
    return (y_val, y_pred)

def generate_scatterplot(y_val: np.ndarray, y_pred: np.ndarray) -> None:
    """
    Takes in the y_val and y_pred numpy arrays and generates a scatterplot. The
    scatterplot is saved to *_scatter.png in the output folder specified by the
    user.
    """
    r, p = pearsonr(y_val, y_pred)
    plt.figure(figsize=(8, 6))
    plt.scatter(y_val, y_pred)
    min_val = min(y_val.min(), y_pred.min())
    max_val = max(y_val.max(), y_pred.max())
    plt.plot([min_val, max_val], [min_val, max_val], "r--")
    plt.title(f"Predicted vs. Actual Values (r={r:.3f})")
    plt.xlabel("Actual")
    plt.ylabel("Predicted")
    plt.savefig(scatter_path, dpi=300)

def generate_residual_plot(y_val: np.ndarray, y_pred: np.ndarray) -> None:
    """
    Takes in the y_val and y_pred numpy arrays and generates a residual
    scatterplot. The scatterplot is saved to *_residual.png in the output folder
    specified by the user.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(y_pred, y_val - y_pred)
    plt.axhline(y=0, color="red", linestyle="--")
    plt.xlabel("Predicted")
    plt.ylabel("Residuals (actual - predicted)")
    plt.title("Residual Plot")
    plt.savefig(residual_path, dpi=300)

def generate_mse_plot(models: list[str], mse: list[float]) -> None:
    """
    Takes in the names of the models in a list and the corresponding list of
    MSEs for each of the models. A bar graph is generated the MSEs of each of
    the models. The plot is saved to mse.png in the output folder specified by
    the user.
    """
    plt.figure(figsize=(8, 6))
    plt.bar(models, mse)
    plt.xlabel("Model")
    plt.ylabel("Mean Squared Error")
    plt.title("MSE Comparison Between Models")
    plt.savefig(mse_path, dpi=300)

if __name__ == "__main__":

    '''
    Takes in one argument, the CSV name for each type of model evaluated.
    Last argument will be the directory which all the plots will be stored.
    Assuming CSV is shaped so that the first column is y_val and the second is
    y_pred (will change this later if needed).
    Also assuming the values for y_val and y_pred are numerical floats (will
    convert this later if needed)

    Run with test data:
    python src/plotting_modules.py data/plotting_testing/LR.csv 
    data/plotting_testing/random_forest.csv data/plotting_testing/svr.csv 
    data/plotting_testing/XGBoost.csv data/plotting_testing/results
    '''

    number_of_models = len(sys.argv) - 2
    directory = sys.argv[len(sys.argv) - 1] # assuming no / at end
    directory_path = Path(directory)
    directory_path.mkdir(parents=True, exist_ok=True)
    models = []
    mse = []

    for i in range(number_of_models):
        filename = sys.argv[i+1]
        file_prefix = filename.split('/')[-1].split('.')[0] # Grab name before extension
        y_val, y_pred = parse_csv(filename)

        # Actual vs. predicted scatter plot
        scatter_path = directory + "/" + file_prefix + "_scatter.png"
        scatter_path = Path(scatter_path)
        generate_scatterplot(y_val, y_pred)

        # Residual plot (y_val - y_pred)
        residual_path = directory + "/" + file_prefix + "_residual.png"
        residual_path = Path(residual_path)
        generate_residual_plot(y_val, y_pred)

        # MSE -> storing in array
        models.append(file_prefix)
        mse.append(mean_squared_error(y_val, y_pred))

    # MSE
    mse_path = directory + "/" + "mse.png"
    mse_path = Path(mse_path)
    generate_mse_plot(models, mse)
