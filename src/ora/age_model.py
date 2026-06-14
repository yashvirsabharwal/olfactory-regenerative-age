"""Donor-level ORA model training."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


TECHNICAL_PREFIXES = ("total_", "log10_total_")
TECHNICAL_EXACT = {
    "age",
    "sex",
    "race_ethnicity",
    "disease",
    "disease_group",
    "chemistry",
    "collection_method",
    "site",
    "sample_id",
    "has_age",
    "is_healthy",
    "is_ndd",
    "is_training_donor",
    "usable_for_ora_training",
    "lineage_cells",
    "mature_neurons",
}

MODEL_ORDER = ["null_model", "ridge", "lasso", "elastic_net", "random_forest"]


@dataclass
class ModelResult:
    performance: pd.DataFrame
    predictions: pd.DataFrame
    feature_importance: pd.DataFrame


@dataclass
class ProjectionResult:
    predictions: pd.DataFrame
    summary: pd.DataFrame


@dataclass
class RepeatedModelResult:
    repeat_performance: pd.DataFrame
    performance_summary: pd.DataFrame
    predictions: pd.DataFrame
    feature_stability: pd.DataFrame


def biological_feature_columns(features: pd.DataFrame, model_config: dict[str, Any] | None = None) -> list[str]:
    """Return biological numeric features, excluding technical/yield columns."""

    excluded = set(TECHNICAL_EXACT)
    for item in (model_config or {}).get("exclude_from_biological_features", []):
        excluded.add(str(item))
    cols: list[str] = []
    for col in features.columns:
        if col in {"donor_id"} or col in excluded:
            continue
        if any(col.startswith(prefix) for prefix in TECHNICAL_PREFIXES):
            continue
        if pd.api.types.is_numeric_dtype(features[col]):
            cols.append(col)
    return cols


def train_ora_models(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
) -> ModelResult:
    """Train composition-MVP ORA models with donor-level folds."""

    model_config = model_config or {}
    donor_meta = (
        manifest.sort_values(["donor_id", "sample_id"])
        .drop_duplicates("donor_id")
        .set_index("donor_id")
    )
    data = features.set_index("donor_id").join(donor_meta, how="inner", rsuffix="_meta").reset_index()
    train_mask = _boolean_series(data["usable_for_ora_training"]) & data["age"].notna()
    train = data[train_mask].copy()
    train = train.reset_index(drop=True)
    feature_cols = biological_feature_columns(train, model_config)
    if not feature_cols:
        raise ValueError("No biological numeric feature columns available for ORA modeling.")

    max_missing = float(model_config.get("missingness_max_fraction", 0.30))
    missing_fraction = train[feature_cols].isna().mean()
    feature_cols = [col for col in feature_cols if missing_fraction[col] <= max_missing]
    if not feature_cols:
        raise ValueError("All feature columns exceeded the missingness threshold.")

    y = train["age"].astype(float).to_numpy()
    folds = donor_cv_folds(train, model_config)
    predictions = []
    importances = []
    performance_rows = []

    for model_name in MODEL_ORDER:
        pred = np.full(train.shape[0], np.nan, dtype=float)
        fold_importances = []
        for fold_id, (train_idx, test_idx) in enumerate(folds):
            prep = fit_preprocessor(train.iloc[train_idx][feature_cols])
            x_train = transform_preprocessor(train.iloc[train_idx][feature_cols], prep)
            x_test = transform_preprocessor(train.iloc[test_idx][feature_cols], prep)
            y_train = y[train_idx]
            if model_name == "null_model":
                fold_pred = np.full(test_idx.size, float(np.mean(y_train)))
                coefs = np.zeros(len(feature_cols))
            elif model_name == "ridge":
                fold_pred, coefs = _fit_ridge_or_linear(x_train, y_train, x_test, model_config)
            elif model_name == "lasso":
                fold_pred, coefs = _fit_lasso_or_linear(x_train, y_train, x_test, model_config)
            elif model_name == "elastic_net":
                fold_pred, coefs = _fit_elastic_or_linear(x_train, y_train, x_test, model_config)
            else:
                fold_pred, coefs = _fit_random_forest_or_linear(x_train, y_train, x_test, model_config)
            pred[test_idx] = fold_pred
            if coefs is not None:
                fold_importances.append(pd.DataFrame({"feature": feature_cols, "importance": coefs, "fold": fold_id}))
        performance_rows.append(_performance_row(model_name, y, pred))
        predictions.append(
            pd.DataFrame(
                {
                    "donor_id": train["donor_id"].to_numpy(),
                    "model": model_name,
                    "chronological_age": y,
                    "ora": pred,
                }
            )
        )
        if fold_importances:
            imp = pd.concat(fold_importances, ignore_index=True)
            summary = imp.groupby("feature", as_index=False).agg(
                importance=("importance", "mean"),
                stability=("importance", lambda s: float(np.mean(np.sign(s) == np.sign(np.mean(s))))),
            )
            summary.insert(0, "model", model_name)
            importances.append(summary)

    pred_table = pd.concat(predictions, ignore_index=True)
    pred_table = add_oraa(pred_table, train)
    perf = pd.DataFrame(performance_rows)
    importance = pd.concat(importances, ignore_index=True) if importances else pd.DataFrame()
    return ModelResult(performance=perf, predictions=pred_table, feature_importance=importance)


def train_ora_models_repeated(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
    *,
    repeats: int | None = None,
) -> RepeatedModelResult:
    """Run repeated donor-level CV and summarize model uncertainty."""

    model_config = dict(model_config or {})
    repeats = int(repeats or model_config.get("outer_cv_repeats", 1))
    repeats = max(1, repeats)
    base_seed = int(model_config.get("random_seed", 42))
    performance_rows = []
    prediction_rows = []
    importance_rows = []
    for repeat in range(repeats):
        repeat_config = dict(model_config)
        repeat_config["random_seed"] = base_seed + repeat
        result = train_ora_models(features, manifest, repeat_config)
        perf = result.performance.copy()
        perf.insert(0, "repeat", repeat)
        performance_rows.append(perf)
        preds = result.predictions.copy()
        preds.insert(0, "repeat", repeat)
        prediction_rows.append(preds)
        if not result.feature_importance.empty:
            imp = result.feature_importance.copy()
            imp.insert(0, "repeat", repeat)
            importance_rows.append(imp)

    repeat_performance = pd.concat(performance_rows, ignore_index=True)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    feature_importance = pd.concat(importance_rows, ignore_index=True) if importance_rows else pd.DataFrame()
    return RepeatedModelResult(
        repeat_performance=repeat_performance,
        performance_summary=summarize_repeated_performance(repeat_performance),
        predictions=predictions,
        feature_stability=summarize_feature_stability(feature_importance),
    )


def summarize_repeated_performance(performance: pd.DataFrame) -> pd.DataFrame:
    """Summarize repeated-CV metrics with empirical 95% intervals."""

    if performance.empty:
        return pd.DataFrame()
    metrics = ["mae", "rmse", "r2", "spearman_r"]
    rows = []
    for model, frame in performance.groupby("model", observed=True, sort=False):
        row: dict[str, object] = {
            "model": model,
            "repeats": int(frame["repeat"].nunique()) if "repeat" in frame else 1,
            "n": int(frame["n"].median()) if "n" in frame else np.nan,
        }
        for metric in metrics:
            values = pd.to_numeric(frame.get(metric), errors="coerce").dropna().to_numpy(dtype=float)
            if values.size == 0:
                row[f"{metric}_mean"] = np.nan
                row[f"{metric}_sd"] = np.nan
                row[f"{metric}_ci_low"] = np.nan
                row[f"{metric}_ci_high"] = np.nan
                continue
            row[f"{metric}_mean"] = float(np.mean(values))
            row[f"{metric}_sd"] = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
            row[f"{metric}_ci_low"] = float(np.quantile(values, 0.025))
            row[f"{metric}_ci_high"] = float(np.quantile(values, 0.975))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_feature_stability(feature_importance: pd.DataFrame) -> pd.DataFrame:
    """Summarize feature importance stability across repeated CV runs."""

    if feature_importance.empty:
        return pd.DataFrame()
    frame = feature_importance.copy()
    frame["importance"] = pd.to_numeric(frame["importance"], errors="coerce").fillna(0.0)
    frame["selected"] = frame["importance"].abs().gt(1e-12)
    summary = (
        frame.groupby(["model", "feature"], observed=True)
        .agg(
            mean_importance=("importance", "mean"),
            sd_importance=("importance", "std"),
            selection_fraction=("selected", "mean"),
            repeats=("repeat", "nunique"),
        )
        .reset_index()
    )
    summary["abs_mean_importance"] = summary["mean_importance"].abs()
    return summary.sort_values(["model", "selection_fraction", "abs_mean_importance"], ascending=[True, False, False]).reset_index(drop=True)


def project_ora_models(
    features: pd.DataFrame,
    manifest: pd.DataFrame,
    model_config: dict[str, Any] | None = None,
) -> ProjectionResult:
    """Train frozen ORA models on healthy donors and project all donors with features."""

    model_config = model_config or {}
    donor_meta = (
        manifest.sort_values(["donor_id", "sample_id"])
        .drop_duplicates("donor_id")
        .set_index("donor_id")
    )
    data = features.set_index("donor_id").join(donor_meta, how="inner", rsuffix="_meta").reset_index()
    data["is_training_donor"] = _boolean_series(data["usable_for_ora_training"]) & data["age"].notna()
    train = data[data["is_training_donor"]].copy().reset_index(drop=True)
    project = data.copy().reset_index(drop=True)
    if train.shape[0] < 2:
        raise ValueError("At least two healthy age-known donors are required for frozen ORA projection.")

    feature_cols = biological_feature_columns(train, model_config)
    if not feature_cols:
        raise ValueError("No biological numeric feature columns available for ORA projection.")
    max_missing = float(model_config.get("missingness_max_fraction", 0.30))
    missing_fraction = train[feature_cols].isna().mean()
    feature_cols = [col for col in feature_cols if missing_fraction[col] <= max_missing]
    if not feature_cols:
        raise ValueError("All feature columns exceeded the missingness threshold.")

    prep = fit_preprocessor(train[feature_cols])
    x_train = transform_preprocessor(train[feature_cols], prep)
    x_project = transform_preprocessor(project[feature_cols], prep)
    y_train = train["age"].astype(float).to_numpy()
    rows = []
    for model_name in MODEL_ORDER:
        if model_name == "null_model":
            pred = np.full(project.shape[0], float(np.mean(y_train)))
        elif model_name == "ridge":
            pred, _ = _fit_ridge_or_linear(x_train, y_train, x_project, model_config)
        elif model_name == "lasso":
            pred, _ = _fit_lasso_or_linear(x_train, y_train, x_project, model_config)
        elif model_name == "elastic_net":
            pred, _ = _fit_elastic_or_linear(x_train, y_train, x_project, model_config)
        else:
            pred, _ = _fit_random_forest_or_linear(x_train, y_train, x_project, model_config)
        frame = _projection_frame(project, model_name, pred, train.shape[0], len(feature_cols))
        rows.append(frame)

    predictions = pd.concat(rows, ignore_index=True)
    predictions = add_frozen_oraa(predictions)
    summary = summarize_projection(predictions)
    return ProjectionResult(predictions=predictions, summary=summary)


def donor_cv_folds(data: pd.DataFrame, model_config: dict[str, Any] | None = None) -> list[tuple[np.ndarray, np.ndarray]]:
    """Create deterministic donor-level folds stratified approximately by age bins."""

    model_config = model_config or {}
    n_splits = int(model_config.get("outer_cv_folds", 5))
    random_seed = int(model_config.get("random_seed", 42))
    n = data.shape[0]
    if n < 2:
        raise ValueError("At least two donors are required for cross-validation.")
    n_splits = max(2, min(n_splits, n))
    rng = np.random.default_rng(random_seed)
    sortable = data[["donor_id", "age"]].copy()
    sortable["_bin"] = _age_bins(sortable["age"], model_config)
    fold_indices: list[list[int]] = [[] for _ in range(n_splits)]
    for _, group in sortable.groupby("_bin", dropna=False):
        indices = group.index.to_numpy()
        rng.shuffle(indices)
        for pos, idx in enumerate(indices):
            fold_indices[pos % n_splits].append(int(idx))
    folds = []
    all_idx = np.arange(n)
    for test in fold_indices:
        test_idx = np.array(sorted(test), dtype=int)
        train_idx = np.setdiff1d(all_idx, test_idx)
        if test_idx.size and train_idx.size:
            folds.append((train_idx, test_idx))
    return folds


def fit_preprocessor(frame: pd.DataFrame) -> dict[str, pd.Series]:
    medians = frame.median(numeric_only=True).fillna(0)
    filled = frame.fillna(medians)
    means = filled.mean()
    stds = filled.std(ddof=0).replace(0, 1).fillna(1)
    return {"medians": medians, "means": means, "stds": stds}


def transform_preprocessor(frame: pd.DataFrame, prep: dict[str, pd.Series]) -> np.ndarray:
    filled = frame.fillna(prep["medians"])
    scaled = (filled - prep["means"]) / prep["stds"]
    return scaled.to_numpy(dtype=float)


def add_oraa(predictions: pd.DataFrame, train_meta: pd.DataFrame) -> pd.DataFrame:
    output = predictions.copy()
    meta = train_meta[["donor_id", "age", "sex", "chemistry", "collection_method", "site"]].copy()
    output = output.merge(meta, on="donor_id", how="left")
    output["oraa"] = np.nan
    for model, frame in output.groupby("model"):
        idx = frame.index.to_numpy()
        y = frame["ora"].to_numpy(dtype=float)
        x_parts = [np.ones(frame.shape[0]), frame["age"].to_numpy(dtype=float)]
        for cov in ["sex", "chemistry", "collection_method", "site"]:
            if cov in frame and frame[cov].notna().nunique() > 1:
                dummies = pd.get_dummies(frame[cov].astype(str), drop_first=True, dtype=float)
                x_parts.extend(dummies[col].to_numpy(dtype=float) for col in dummies.columns)
        x = np.vstack(x_parts).T
        valid = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if valid.sum() <= x.shape[1]:
            output.loc[idx, "oraa"] = y - frame["age"].to_numpy(dtype=float)
            continue
        coef, *_ = np.linalg.lstsq(x[valid], y[valid], rcond=None)
        expected = x @ coef
        output.loc[idx, "oraa"] = y - expected
    return output


def add_frozen_oraa(predictions: pd.DataFrame) -> pd.DataFrame:
    """Residualize projected ORA against age/covariates using healthy training donors only."""

    output = predictions.copy()
    output["oraa"] = np.nan
    for model, frame in output.groupby("model", sort=False):
        idx = frame.index.to_numpy()
        age = pd.to_numeric(frame.get("chronological_age"), errors="coerce")
        ora = pd.to_numeric(frame.get("ora"), errors="coerce")
        train_mask = _boolean_series(frame.get("is_training_donor", pd.Series(False, index=frame.index)))
        valid_ref = train_mask & age.notna() & ora.notna()
        if valid_ref.sum() < 2:
            valid = age.notna() & ora.notna()
            output.loc[idx[valid.to_numpy()], "oraa"] = ora[valid] - age[valid]
            continue
        ref = frame.loc[valid_ref].copy()
        x_ref, columns = _expected_ora_design(ref)
        y_ref = pd.to_numeric(ref["ora"], errors="coerce").to_numpy(dtype=float)
        valid_design = np.isfinite(x_ref).all(axis=1) & np.isfinite(y_ref)
        if valid_design.sum() <= x_ref.shape[1]:
            valid = age.notna() & ora.notna()
            output.loc[idx[valid.to_numpy()], "oraa"] = ora[valid] - age[valid]
            continue
        coef, *_ = np.linalg.lstsq(x_ref[valid_design], y_ref[valid_design], rcond=None)
        x_all, _ = _expected_ora_design(frame, columns=columns)
        valid_all = age.notna().to_numpy() & ora.notna().to_numpy() & np.isfinite(x_all).all(axis=1)
        output.loc[idx[valid_all], "oraa"] = ora.to_numpy(dtype=float)[valid_all] - x_all[valid_all] @ coef
    return output


def summarize_projection(predictions: pd.DataFrame) -> pd.DataFrame:
    """Summarize projected ORA and ORAA by model and disease group."""

    if predictions.empty:
        return pd.DataFrame()
    frame = predictions.copy()
    for col in ["chronological_age", "ora", "oraa"]:
        frame[col] = pd.to_numeric(frame.get(col), errors="coerce")
    group_cols = ["model", "disease_group"]
    summary = (
        frame.groupby(group_cols, observed=True, dropna=False)
        .agg(
            donors=("donor_id", "nunique"),
            training_donors=("is_training_donor", lambda s: int(_boolean_series(s).sum())),
            ndd_donors=("is_ndd", lambda s: int(_boolean_series(s).sum())),
            mean_age=("chronological_age", "mean"),
            mean_ora=("ora", "mean"),
            mean_oraa=("oraa", "mean"),
            sd_oraa=("oraa", lambda s: float(pd.to_numeric(s, errors="coerce").std(ddof=0))),
        )
        .reset_index()
        .sort_values(["model", "disease_group"])
    )
    return summary


def _fit_elastic_or_linear(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    model_config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from sklearn.linear_model import ElasticNetCV  # type: ignore

        params = model_config.get("models", {}).get("elastic_net", {})
        model = ElasticNetCV(
            alphas=params.get("alphas", [0.01, 0.1, 1.0]),
            l1_ratio=params.get("l1_ratios", [0.1, 0.5, 0.9]),
            cv=min(5, max(2, len(y_train) // 4)),
            random_state=int(model_config.get("random_seed", 42)),
            max_iter=int(params.get("max_iter", 50000)),
            tol=float(params.get("tol", 1e-2)),
            n_jobs=int(params.get("n_jobs", -1)),
        )
        model.fit(x_train, y_train)
        return model.predict(x_test), np.asarray(model.coef_, dtype=float)
    except ModuleNotFoundError:
        coef, intercept = _ridge_closed_form(x_train, y_train)
        return x_test @ coef + intercept, coef


def _fit_ridge_or_linear(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    model_config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from sklearn.linear_model import RidgeCV  # type: ignore

        params = model_config.get("models", {}).get("ridge", {})
        model = RidgeCV(alphas=params.get("alphas", [0.1, 1.0, 10.0, 100.0]))
        model.fit(x_train, y_train)
        return model.predict(x_test), np.asarray(model.coef_, dtype=float)
    except ModuleNotFoundError:
        coef, intercept = _ridge_closed_form(x_train, y_train, alpha=10.0)
        return x_test @ coef + intercept, coef


def _fit_lasso_or_linear(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    model_config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from sklearn.linear_model import LassoCV  # type: ignore

        params = model_config.get("models", {}).get("lasso", {})
        model = LassoCV(
            alphas=params.get("alphas", [0.001, 0.01, 0.1, 1.0, 10.0]),
            cv=min(5, max(2, len(y_train) // 4)),
            random_state=int(model_config.get("random_seed", 42)),
            max_iter=int(params.get("max_iter", 50000)),
            tol=float(params.get("tol", 1e-3)),
            n_jobs=int(params.get("n_jobs", -1)),
        )
        model.fit(x_train, y_train)
        return model.predict(x_test), np.asarray(model.coef_, dtype=float)
    except ModuleNotFoundError:
        coef, intercept = _ridge_closed_form(x_train, y_train, alpha=25.0)
        return x_test @ coef + intercept, coef


def _fit_random_forest_or_linear(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    model_config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from sklearn.ensemble import RandomForestRegressor  # type: ignore

        params = model_config.get("models", {}).get("random_forest", {})
        model = RandomForestRegressor(
            n_estimators=int(params.get("n_estimators", 500)),
            max_depth=params.get("max_depth", 5),
            random_state=int(model_config.get("random_seed", 42)),
            n_jobs=-1,
        )
        model.fit(x_train, y_train)
        return model.predict(x_test), np.asarray(model.feature_importances_, dtype=float)
    except ModuleNotFoundError:
        coef, intercept = _ridge_closed_form(x_train, y_train, alpha=10.0)
        return x_test @ coef + intercept, np.abs(coef)


def _ridge_closed_form(x: np.ndarray, y: np.ndarray, alpha: float = 1.0) -> tuple[np.ndarray, float]:
    y_mean = float(np.mean(y))
    centered_y = y - y_mean
    penalty = alpha * np.eye(x.shape[1])
    coef = np.linalg.pinv(x.T @ x + penalty) @ x.T @ centered_y
    return coef, y_mean


def _performance_row(model_name: str, y: np.ndarray, pred: np.ndarray) -> dict[str, float | str]:
    valid = np.isfinite(y) & np.isfinite(pred)
    yv = y[valid]
    pv = pred[valid]
    resid = yv - pv
    mae = float(np.mean(np.abs(resid)))
    rmse = float(np.sqrt(np.mean(resid**2)))
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    corr = _spearman(yv, pv) if yv.size > 1 else np.nan
    return {
        "model": model_name,
        "n": int(yv.size),
        "mae": mae,
        "rmse": rmse,
        "r2": r2,
        "spearman_r": corr,
    }


def _age_bins(age: pd.Series, model_config: dict[str, Any]) -> pd.Series:
    bins_cfg = model_config.get("age_bins", {})
    output = pd.Series("unknown", index=age.index, dtype=object)
    for label, bounds in bins_cfg.items():
        if len(bounds) != 2:
            continue
        low, high = bounds
        output[(age >= low) & (age <= high)] = label
    return output


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ar = pd.Series(a).rank(method="average").to_numpy(dtype=float)
    br = pd.Series(b).rank(method="average").to_numpy(dtype=float)
    if np.std(ar) == 0 or np.std(br) == 0:
        return np.nan
    return float(np.corrcoef(ar, br)[0, 1])


def _projection_frame(
    project: pd.DataFrame,
    model_name: str,
    pred: np.ndarray,
    training_n: int,
    n_features: int,
) -> pd.DataFrame:
    meta_cols = [
        "donor_id",
        "age",
        "disease",
        "disease_group",
        "sex",
        "race_ethnicity",
        "chemistry",
        "collection_method",
        "site",
        "total_cells",
        "is_healthy",
        "is_ndd",
        "usable_for_ora_training",
        "is_training_donor",
    ]
    available = [col for col in meta_cols if col in project.columns]
    frame = project[available].copy()
    frame.insert(1, "model", model_name)
    frame["chronological_age"] = pd.to_numeric(frame.get("age"), errors="coerce")
    frame["ora"] = pred
    frame["training_n"] = training_n
    frame["n_features"] = n_features
    if "disease_group" not in frame:
        frame["disease_group"] = "unknown"
    if "is_ndd" not in frame:
        frame["is_ndd"] = frame["disease_group"].astype(str).isin(["ad", "pd"])
    return frame


def _expected_ora_design(
    frame: pd.DataFrame,
    columns: list[str] | None = None,
) -> tuple[np.ndarray, list[str]]:
    pieces = [pd.Series(1.0, index=frame.index, name="intercept")]
    pieces.append(pd.to_numeric(frame.get("chronological_age"), errors="coerce").rename("chronological_age"))
    if columns is None:
        dummy_pieces = []
        for cov in ["sex", "chemistry", "collection_method", "site"]:
            if cov in frame and frame[cov].notna().nunique() > 1:
                dummy_pieces.append(pd.get_dummies(frame[cov].fillna("unknown").astype(str), prefix=cov, drop_first=True, dtype=float))
        design = pd.concat([*pieces, *dummy_pieces], axis=1)
        return design.to_numpy(dtype=float), list(design.columns)
    design = pd.concat(pieces, axis=1)
    for col in columns:
        if col in design.columns:
            continue
        design[col] = 0.0
        for cov in ["sex", "chemistry", "collection_method", "site"]:
            prefix = f"{cov}_"
            if col.startswith(prefix) and cov in frame:
                value = col[len(prefix):]
                design[col] = (frame[cov].fillna("unknown").astype(str) == value).astype(float)
                break
    design = design.reindex(columns=columns, fill_value=0.0)
    return design.to_numpy(dtype=float), columns


def _boolean_series(values: Any) -> pd.Series:
    series = values if isinstance(values, pd.Series) else pd.Series(values)
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    return normalized.isin({"1", "true", "t", "yes", "y"})
