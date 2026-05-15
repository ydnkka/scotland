from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .config import AnalysisConfig


@dataclass(slots=True)
class SplitData:
    X_train: pd.DataFrame
    X_test: pd.DataFrame
    y_train: pd.Series
    y_test: pd.Series
    groups_train: pd.Series
    groups_test: pd.Series
    train_idx: np.ndarray
    test_idx: np.ndarray
    splitter_name: str


@dataclass(slots=True)
class BinaryResult:
    name: str
    model: Pipeline
    metrics: dict
    confusion_matrix: pd.DataFrame
    classification_report: pd.DataFrame
    probabilities: np.ndarray
    predictions: np.ndarray
    importance: pd.DataFrame | None = None


@dataclass(slots=True)
class MulticlassResult:
    name: str
    model: Pipeline
    metrics: dict
    confusion_matrix: pd.DataFrame
    classification_report: pd.DataFrame
    probabilities: np.ndarray
    predictions: np.ndarray


def make_one_hot_encoder() -> OneHotEncoder:
    """Create a version-compatible one-hot encoder."""
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=True, min_frequency=20)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=True, min_frequency=20)


def make_rf_pipeline(
    config: AnalysisConfig,
    numeric_features: list[str],
    categorical_features: list[str],
) -> Pipeline:
    numeric_pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
        ]
    )
    categorical_pipe = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="constant", fill_value="Missing")),
            ("onehot", make_one_hot_encoder()),
        ]
    )
    preprocessor = ColumnTransformer(
        transformers=[
            ("num", numeric_pipe, numeric_features),
            ("cat", categorical_pipe, categorical_features),
        ],
        remainder="drop",
        sparse_threshold=0.3,
        verbose_feature_names_out=False,
    )
    rf = RandomForestClassifier(
        n_estimators=config.n_estimators,
        max_features="sqrt",
        min_samples_leaf=config.min_samples_leaf,
        class_weight="balanced_subsample",
        n_jobs=-1,
        random_state=config.random_state,
        oob_score=False,
    )
    return Pipeline(
        [
            ("preprocessor", preprocessor),
            ("rf", rf),
        ]
    )


def grouped_holdout_indices(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    config: AnalysisConfig,
    test_size_approx: float = 0.2,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Build a grouped holdout split, stratified when supported."""
    n_splits = round(1 / test_size_approx)
    try:
        splitter = StratifiedGroupKFold(
            n_splits=n_splits, shuffle=True, random_state=config.random_state
        )
        train_idx, test_idx = next(splitter.split(X, y, groups))
        split_name = f"StratifiedGroupKFold(n_splits={n_splits})"
    except Exception:
        splitter = GroupShuffleSplit(
            n_splits=1, test_size=test_size_approx, random_state=config.random_state
        )
        train_idx, test_idx = next(splitter.split(X, y, groups))
        split_name = "GroupShuffleSplit"
    return train_idx, test_idx, split_name


def make_split_data(
    X: pd.DataFrame,
    y: pd.Series,
    groups: pd.Series,
    config: AnalysisConfig,
    test_size_approx: float = 0.2,
) -> SplitData:
    """Split into grouped train/test sets and optionally downsample training clusters."""
    train_idx, test_idx, splitter_name = grouped_holdout_indices(
        X, y, groups, config, test_size_approx=test_size_approx
    )
    X_train = X.iloc[train_idx].copy()
    X_test = X.iloc[test_idx].copy()
    y_train = y.iloc[train_idx].copy()
    y_test = y.iloc[test_idx].copy()
    groups_train = groups.iloc[train_idx].copy()
    groups_test = groups.iloc[test_idx].copy()

    if config.max_train_rows is not None and len(X_train) > config.max_train_rows:
        rng = np.random.default_rng(config.random_state)
        sampled_clusters: list[str] = []
        n_rows = 0
        for cluster_id in rng.permutation(groups_train.unique()):
            mask_count = int((groups_train == cluster_id).sum())
            sampled_clusters.append(cluster_id)
            n_rows += mask_count
            if n_rows >= config.max_train_rows:
                break
        keep = groups_train.isin(sampled_clusters).to_numpy()
        X_train = X_train.loc[keep].copy()
        y_train = y_train.loc[keep].copy()
        groups_train = groups_train.loc[keep].copy()

    if not set(groups_train).isdisjoint(set(groups_test)):
        raise RuntimeError("Cluster leakage: at least one cluster appears in train and test.")

    return SplitData(
        X_train=X_train,
        X_test=X_test,
        y_train=y_train,
        y_test=y_test,
        groups_train=groups_train,
        groups_test=groups_test,
        train_idx=train_idx,
        test_idx=test_idx,
        splitter_name=splitter_name,
    )


def evaluate_binary_model(
    model: Pipeline,
    X_eval: pd.DataFrame,
    y_eval: pd.Series,
    label: str = "holdout",
    negative_name: str = "not_large",
    positive_name: str = "large",
) -> tuple[dict, pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    pred = model.predict(X_eval)
    proba = model.predict_proba(X_eval)[:, 1]
    metrics = {
        "label": label,
        "n_rows": len(y_eval),
        "balanced_accuracy": balanced_accuracy_score(y_eval, pred),
        "macro_f1": f1_score(y_eval, pred, average="macro"),
        "roc_auc": roc_auc_score(y_eval, proba),
        "average_precision": average_precision_score(y_eval, proba),
    }
    cm = pd.DataFrame(
        confusion_matrix(y_eval, pred, labels=[0, 1]),
        index=[f"true_{negative_name}", f"true_{positive_name}"],
        columns=[f"pred_{negative_name}", f"pred_{positive_name}"],
    )
    report = pd.DataFrame(
        classification_report(
            y_eval,
            pred,
            target_names=[negative_name, positive_name],
            digits=3,
            output_dict=True,
        )
    ).T
    return metrics, cm, report, proba, pred


def fit_binary_model(
    name: str,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    config: AnalysisConfig,
    numeric_features: list[str],
    categorical_features: list[str],
    negative_name: str,
    positive_name: str,
) -> BinaryResult:
    model = make_rf_pipeline(config, numeric_features, categorical_features)
    model.fit(X_train, y_train)
    metrics, cm, report, proba, pred = evaluate_binary_model(
        model,
        X_test,
        y_test,
        label=f"{name}_holdout",
        negative_name=negative_name,
        positive_name=positive_name,
    )
    return BinaryResult(
        name=name,
        model=model,
        metrics=metrics,
        confusion_matrix=cm,
        classification_report=report,
        probabilities=proba,
        predictions=pred,
    )


def evaluate_multiclass_model(
    model: Pipeline,
    X_eval: pd.DataFrame,
    y_eval: pd.Series,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
    pred = model.predict(X_eval)
    proba = model.predict_proba(X_eval)
    classes = list(model.named_steps["rf"].classes_)
    metrics = {
        "n_rows": len(y_eval),
        "balanced_accuracy": balanced_accuracy_score(y_eval, pred),
        "macro_f1": f1_score(y_eval, pred, average="macro"),
        "weighted_f1": f1_score(y_eval, pred, average="weighted"),
        "macro_ovr_auc": roc_auc_score(y_eval, proba, multi_class="ovr", average="macro"),
    }
    cm = pd.DataFrame(
        confusion_matrix(y_eval, pred, labels=classes),
        index=[f"true_{c}" for c in classes],
        columns=[f"pred_{c}" for c in classes],
    )
    report = pd.DataFrame(classification_report(y_eval, pred, output_dict=True, digits=3)).T
    return metrics, cm, report, proba, pred


def fit_multiclass_model(
    name: str,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    config: AnalysisConfig,
    numeric_features: list[str],
    categorical_features: list[str],
) -> MulticlassResult:
    model = make_rf_pipeline(config, numeric_features, categorical_features)
    model.fit(X_train, y_train)
    metrics, cm, report, proba, pred = evaluate_multiclass_model(model, X_test, y_test)
    return MulticlassResult(
        name=name,
        model=model,
        metrics=metrics,
        confusion_matrix=cm,
        classification_report=report,
        probabilities=proba,
        predictions=pred,
    )


def sample_for_importance(
    X_eval: pd.DataFrame,
    y_eval: pd.Series,
    config: AnalysisConfig,
) -> tuple[pd.DataFrame, pd.Series]:
    if config.max_perm_rows is None or len(X_eval) <= config.max_perm_rows:
        return X_eval, y_eval
    rng = np.random.default_rng(config.random_state)
    idx = rng.choice(len(X_eval), size=config.max_perm_rows, replace=False)
    return X_eval.iloc[idx].copy(), y_eval.iloc[idx].copy()


def permutation_importance_table(
    model: Pipeline,
    X_eval: pd.DataFrame,
    y_eval: pd.Series,
    features: list[str],
    config: AnalysisConfig,
    scoring: str = "average_precision",
) -> pd.DataFrame:
    X_perm, y_perm = sample_for_importance(X_eval, y_eval, config)
    perm = permutation_importance(
        model,
        X_perm,
        y_perm,
        scoring=scoring,
        n_repeats=config.perm_repeats,
        random_state=config.random_state,
        n_jobs=1,
    )
    return (
        pd.DataFrame(
            {
                "feature": features,
                f"mean_drop_{scoring}": perm.importances_mean,
                f"std_drop_{scoring}": perm.importances_std,
            }
        )
        .sort_values(f"mean_drop_{scoring}", ascending=False)
        .reset_index(drop=True)
    )
