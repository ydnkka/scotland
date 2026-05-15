from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.inspection import PartialDependenceDisplay
from sklearn.metrics import precision_recall_curve


def plot_outcome_counts(outcome_counts: pd.DataFrame, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    sns.barplot(
        data=outcome_counts,
        x="cluster_type",
        y="n_sequences",
        order=["singleton", "small_cluster", "large_cluster"],
        ax=ax,
        color="#4477AA",
    )
    ax.set_title("Sequence-level outcome counts")
    ax.set_xlabel("Cluster type")
    ax.set_ylabel("Sequences")
    return ax


def plot_confusion_matrix(cm: pd.DataFrame, title: str, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(cm, annot=True, fmt=",d", cmap="Blues", ax=ax)
    ax.set_title(title)
    return ax


def plot_precision_recall(y_true, y_score, title: str, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    ax.plot(recall, precision)
    ax.set_title(title)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax


def plot_importance(
    importance: pd.DataFrame,
    title: str,
    score_col: str | None = None,
    top_n: int = 25,
    color: str = "#4477AA",
    ax=None,
):
    if score_col is None:
        score_col = next(c for c in importance.columns if c.startswith("mean_drop_"))
    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, 0.35 * min(top_n, len(importance)))))
    sns.barplot(
        data=importance.head(top_n),
        y="feature",
        x=score_col,
        ax=ax,
        color=color,
    )
    ax.set_title(title)
    ax.set_xlabel(score_col.replace("_", " "))
    ax.set_ylabel("")
    return ax


def plot_binary_diagnostics(result, y_test, title_prefix: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    plot_confusion_matrix(result.confusion_matrix, f"{title_prefix}: confusion matrix", ax=axes[0])
    plot_precision_recall(y_test, result.probabilities, f"{title_prefix}: precision-recall", ax=axes[1])
    fig.tight_layout()
    return fig


def plot_partial_dependence_grid(
    model,
    X,
    features: list[str],
    target: int = 1,
    grid_resolution: int = 30,
):
    n_cols = 3
    n_rows = (len(features) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 3.6 * n_rows))
    axes = axes.ravel()
    for ax, feature in zip(axes, features):
        PartialDependenceDisplay.from_estimator(
            model,
            X,
            [feature],
            target=target,
            kind="average",
            grid_resolution=grid_resolution,
            ax=ax,
        )
        ax.set_title(feature)
    for ax in axes[len(features) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def partial_dependence_curves(
    model,
    X: pd.DataFrame,
    features: list[str],
    grid_resolution: int = 30,
    positive_class_index: int = 1,
) -> pd.DataFrame:
    """Compute PDP curves as mean predicted positive-class probability.

    This deliberately uses direct prediction on a sampled dataframe rather than
    sklearn's partial_dependence helper so the saved curves match the model
    pipeline exactly across sklearn versions.
    """
    rows = []
    for feature in features:
        values = X[feature].dropna()
        if values.empty:
            continue
        if values.nunique() <= grid_resolution:
            grid = np.sort(values.unique())
        else:
            grid = np.unique(np.nanquantile(values, np.linspace(0.05, 0.95, grid_resolution)))
        for value in grid:
            tmp = X.copy()
            tmp[feature] = value
            mean_pred = model.predict_proba(tmp)[:, positive_class_index].mean()
            rows.append(
                {
                    "feature": feature,
                    "grid_value": value,
                    "mean_pred_positive": mean_pred,
                    "n_rows": len(tmp),
                }
            )
    return pd.DataFrame(rows)


def plot_partial_dependence_curves(
    curves: pd.DataFrame,
    value_col: str = "grid_value",
    pred_col: str = "mean_pred_positive",
):
    features = list(curves["feature"].drop_duplicates())
    n_cols = 3
    n_rows = (len(features) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 3.6 * n_rows))
    axes = axes.ravel()
    for ax, feature in zip(axes, features):
        use = curves.loc[curves["feature"].eq(feature)]
        ax.plot(use[value_col], use[pred_col])
        ax.set_title(feature)
        ax.set_xlabel(feature)
        ax.set_ylabel("Mean predicted positive probability")
    for ax in axes[len(features) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def accumulated_local_effects_curves(
    model,
    X: pd.DataFrame,
    features: list[str],
    bins: int = 20,
    positive_class_index: int = 1,
    min_unique: int = 3,
) -> pd.DataFrame:
    """Compute first-order 1D accumulated local effects for numeric features.

    The returned ALE values are centred to have weighted mean zero, so the
    y-axis is the local contribution to predicted positive-class probability
    relative to the feature's average contribution in the sampled data.
    """
    rows = []
    for feature in features:
        if feature not in X.columns:
            continue
        values = pd.to_numeric(X[feature], errors="coerce")
        valid_mask = values.notna()
        if valid_mask.sum() == 0 or values.loc[valid_mask].nunique() < min_unique:
            continue

        X_valid = X.loc[valid_mask].copy()
        x = values.loc[valid_mask].astype(float)

        if x.nunique() <= bins + 1:
            edges = np.sort(x.unique())
        else:
            edges = np.unique(np.nanquantile(x, np.linspace(0, 1, bins + 1)))

        if len(edges) < 2:
            continue

        n_intervals = len(edges) - 1
        interval_ids = np.searchsorted(edges, x.to_numpy(), side="right") - 1
        interval_ids = np.clip(interval_ids, 0, n_intervals - 1)

        effects = []
        counts = []
        x_mid = []
        x_left = []
        x_right = []

        for interval in range(n_intervals):
            in_bin = interval_ids == interval
            n_bin = int(in_bin.sum())
            left = float(edges[interval])
            right = float(edges[interval + 1])
            if n_bin == 0 or left == right:
                effects.append(np.nan)
                counts.append(n_bin)
                x_mid.append((left + right) / 2)
                x_left.append(left)
                x_right.append(right)
                continue

            X_bin_low = X_valid.loc[in_bin].copy()
            X_bin_high = X_valid.loc[in_bin].copy()
            X_bin_low[feature] = left
            X_bin_high[feature] = right

            pred_low = model.predict_proba(X_bin_low)[:, positive_class_index]
            pred_high = model.predict_proba(X_bin_high)[:, positive_class_index]
            effects.append(float(np.mean(pred_high - pred_low)))
            counts.append(n_bin)
            x_mid.append((left + right) / 2)
            x_left.append(left)
            x_right.append(right)

        effects_array = np.asarray(effects, dtype=float)
        counts_array = np.asarray(counts, dtype=float)
        effects_array = np.nan_to_num(effects_array, nan=0.0)
        accumulated = np.cumsum(effects_array)
        if counts_array.sum() > 0:
            accumulated = accumulated - np.average(accumulated, weights=counts_array)
        else:
            accumulated = accumulated - accumulated.mean()

        for interval in range(n_intervals):
            rows.append(
                {
                    "feature": feature,
                    "bin": interval + 1,
                    "x_left": x_left[interval],
                    "x_right": x_right[interval],
                    "x_mid": x_mid[interval],
                    "n_rows": int(counts[interval]),
                    "local_effect": effects_array[interval],
                    "ale": float(accumulated[interval]),
                }
            )

    return pd.DataFrame(rows)


def plot_ale_curves(
    curves: pd.DataFrame,
    x_col: str = "x_mid",
    ale_col: str = "ale",
):
    features = list(curves["feature"].drop_duplicates())
    n_cols = 3
    n_rows = (len(features) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 3.6 * n_rows))
    axes = axes.ravel()
    for ax, feature in zip(axes, features):
        use = curves.loc[curves["feature"].eq(feature)]
        ax.plot(use[x_col], use[ale_col], marker="o", markersize=3)
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.5)
        ax.set_title(feature)
        ax.set_xlabel(feature)
        ax.set_ylabel("ALE on predicted positive probability")
    for ax in axes[len(features) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig
