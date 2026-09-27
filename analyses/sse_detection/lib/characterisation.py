"""Descriptive Chapter 6 sensitivity analyses of saved detector results.

No scores, entropy references or Bayesian models are re-estimated here.
All exported tables are aggregates; sequence identifiers are used only to
validate the unit of observation in memory.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .sse.io import HIGH_PRIORITY_CANDIDATE_TIERS

SIZE_BANDS = ("6-9", "10-19", "20-49", "50-99", "100+")
ATTRIBUTES = (
    ("sex", "sex", "Sex"),
    ("age", "age_group", "Age group"),
    ("simd", "dz_simd_quintile", "SIMD quintile"),
    ("urban_rural", "dz_urban_rural_class", "Settlement type"),
    ("health_board", "dz_health_board", "Health board"),
)
ENTROPY_ATTRIBUTES = (*ATTRIBUTES, ("local_authority", "dz_local_authority", "Local authority"))
SCENARIOS = ("primary", "p_0.010", "p_0.025", "p_0.100", "size_10", "size_20", "route")
ROUTE_LABELS = {
    "high_priority_burst": "Burst only",
    "high_priority_burden": "Burden only",
    "high_priority_both_axes": "Both axes",
}


def candidate_mask(nodes: pd.DataFrame, alpha: float) -> pd.Series:
    """Reapply a cutoff to saved operational p values at the fixed testing floor."""
    return nodes["cluster_size"].ge(6) & nodes["sse_tested"].fillna(False) & (
        nodes["burst_score_upper_p"].le(alpha)
        | (nodes["burden_eligible"].fillna(False) & nodes["burden_score_upper_p"].le(alpha))
    )


def select_groups(nodes: pd.DataFrame, scenario: str) -> dict[str, pd.DataFrame]:
    eligible = nodes.loc[nodes.cluster_size.ge(6)].copy()
    primary = eligible.candidate_tier.isin(HIGH_PRIORITY_CANDIDATE_TIERS)
    if scenario == "route":
        return {
            "Background": eligible.loc[~primary],
            "Burden-eligible background": eligible.loc[~primary & eligible.burden_eligible.fillna(False)],
            **{label: eligible.loc[eligible.candidate_tier.eq(tier)] for tier, label in ROUTE_LABELS.items()},
        }
    if scenario.startswith("size_"):
        eligible = eligible.loc[eligible.cluster_size.ge(int(scenario.split("_")[1]))]
        selected = primary.reindex(eligible.index)
    elif scenario.startswith("p_"):
        selected = candidate_mask(eligible, float(scenario[2:]))
    elif scenario == "primary":
        selected = primary
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    return {"Background": eligible.loc[~selected], "CSC": eligible.loc[selected]}


def category_label(value: object) -> str:
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def build_summaries(nodes: pd.DataFrame, records: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Aggregate size, composition and entropy using explicit non-missing denominators."""
    if nodes.cluster_id.duplicated().any():
        raise ValueError("Detector cluster identifiers must be unique")
    eligible = nodes.loc[nodes.cluster_size.ge(6)].copy()
    if not candidate_mask(eligible, .05).equals(eligible.candidate_tier.isin(HIGH_PRIORITY_CANDIDATE_TIERS)):
        raise ValueError("Saved primary labels disagree with operational p <= 0.05")
    records = records.loc[records.cluster_id.isin(eligible.cluster_id)].copy()
    if records.duplicated(["window_id", "sequence_id"]).any():
        raise ValueError("Expected one record per sequence and retained window")
    sizes = records.groupby("cluster_id", observed=True).size().reindex(eligible.cluster_id, fill_value=0)
    if not np.array_equal(sizes.to_numpy(), eligible.cluster_size.to_numpy()):
        raise ValueError("Loaded sequence records do not reproduce detector cluster sizes")
    eligible["size_band"] = pd.cut(eligible.cluster_size, [5, 9, 19, 49, 99, np.inf], labels=SIZE_BANDS)
    matrices = {}
    for key, col, _ in ATTRIBUTES:
        observed = records.loc[records[col].notna(), ["cluster_id", col]].copy()
        observed[col] = observed[col].map(category_label)
        matrices[key] = pd.crosstab(observed.cluster_id, observed[col]).reindex(eligible.cluster_id, fill_value=0)

    size_rows, composition_rows, entropy_rows, cumulative_rows = [], [], [], []
    for scenario in SCENARIOS:
        for group, block in select_groups(eligible, scenario).items():
            if scenario == "primary":
                counts = block.groupby("cluster_size").size().sort_index()
                cumulative = counts.cumsum()
                mass = (counts * counts.index).cumsum()
                for size in counts.index:
                    cumulative_rows.append(dict(group=group, cluster_size=int(size),
                        cumulative_clusters=int(cumulative.loc[size]),
                        cumulative_records=int(mass.loc[size]),
                        cluster_percent=100*cumulative.loc[size]/len(block),
                        record_percent=100*mass.loc[size]/block.cluster_size.sum()))
            segments = [("All", block)]
            if scenario == "primary":
                segments += [(band, block.loc[block.size_band.eq(band)]) for band in SIZE_BANDS]
            for band, subset in segments:
                n, mass = len(subset), int(subset.cluster_size.sum())
                base = dict(scenario=scenario, group=group, size_band=band)
                q = subset.cluster_size.quantile([.25, .5, .75])
                size_rows.append(base | dict(n_clusters=n, n_records=mass,
                    size_q25=q.loc[.25], size_median=q.loc[.5], size_q75=q.loc[.75],
                    n_burden_eligible=int(subset.burden_eligible.fillna(False).sum()),
                    clusters_below20=int(subset.cluster_size.lt(20).sum()),
                    records_below20=int(subset.loc[subset.cluster_size.lt(20), "cluster_size"].sum())))
                for key, col, label in ATTRIBUTES:
                    count = matrices[key].reindex(subset.cluster_id)
                    valid = count.sum(axis=1)
                    proportions = count.div(valid.replace(0, np.nan), axis=0)
                    denom = int(valid.sum())
                    for category in count.columns:
                        numerator = int(count[category].sum())
                        composition_rows.append(base | dict(attribute=key, attribute_label=label,
                            category=category, n_clusters=n, n_records=mass,
                            n_nonmissing_records=denom, n_missing_records=mass-denom,
                            n_nonmissing_clusters=int(valid.gt(0).sum()),
                            n_missing_clusters=int(valid.eq(0).sum()), category_records=numerator,
                            record_proportion=numerator/denom if denom else np.nan,
                            cluster_mean_proportion=proportions[category].mean()))
                for key, _, label in ENTROPY_ATTRIBUTES:
                    for scale in ("obs", "z"):
                        values = subset[f"{key}_entropy_{scale}"].replace([np.inf, -np.inf], np.nan)
                        quantiles = values.quantile([.25, .5, .75])
                        entropy_rows.append(base | dict(attribute=key, attribute_label=label, scale=scale,
                            n_clusters=n, n_nonmissing=int(values.notna().sum()), n_missing=int(values.isna().sum()),
                            median=quantiles.loc[.5], q25=quantiles.loc[.25], q75=quantiles.loc[.75],
                            mean=values.mean()))
    tables = {"size_summary": pd.DataFrame(size_rows), "size_cumulative": pd.DataFrame(cumulative_rows),
              "composition": pd.DataFrame(composition_rows), "entropy": pd.DataFrame(entropy_rows)}
    contrasts = []
    for scenario in SCENARIOS:
        pairs = [("CSC", "Background")] if scenario != "route" else [
            ("Burst only", "Background"), ("Burden only", "Burden-eligible background"),
            ("Both axes", "Burden-eligible background")]
        df = tables["composition"].loc[lambda t: t.scenario.eq(scenario)]
        for candidate, background in pairs:
            keys = ["size_band", "attribute", "attribute_label", "category"]
            joined = df.loc[df.group.eq(candidate)].merge(df.loc[df.group.eq(background)],
                on=keys, suffixes=("_candidate", "_background"), validate="one_to_one")
            for weighting, source in [("records", "record_proportion"), ("clusters", "cluster_mean_proportion")]:
                out = joined[keys].copy()
                out["scenario"], out["candidate_group"], out["background_group"] = scenario, candidate, background
                out["weighting"] = weighting
                out["candidate_percent"] = 100*joined[source+"_candidate"]
                out["background_percent"] = 100*joined[source+"_background"]
                out["difference_pp"] = out.candidate_percent-out.background_percent
                for side in ("candidate", "background"):
                    for c in ("n_clusters", "n_records", "n_nonmissing_records", "n_missing_records", "n_nonmissing_clusters"):
                        out[f"{c}_{side}"] = joined[f"{c}_{side}"]
                contrasts.append(out)
    tables["composition_contrasts"] = pd.concat(contrasts, ignore_index=True)
    return tables


def validate_summaries(tables: dict[str, pd.DataFrame], *, scotland_baseline: bool = False) -> None:
    """Reconciliation checks for weighting, cutoffs and descriptive denominators."""
    sizes = tables["size_summary"]
    primary = sizes.loc[sizes.scenario.eq("primary")]
    for group, block in primary.groupby("group"):
        total = block.loc[block.size_band.eq("All")].iloc[0]
        bands = block.loc[block.size_band.ne("All")]
        for col in ("n_clusters", "n_records", "n_burden_eligible"):
            assert bands[col].sum() == total[col], (group, col)
        curve = tables["size_cumulative"].loc[lambda t: t.group.eq(group)]
        for col in ("cumulative_clusters", "cumulative_records", "cluster_percent", "record_percent"):
            assert curve[col].is_monotonic_increasing, col
        assert curve.iloc[-1].cumulative_clusters == total.n_clusters
        assert curve.iloc[-1].cumulative_records == total.n_records
        assert np.allclose(curve.iloc[-1][["cluster_percent", "record_percent"]].astype(float), 100)
    comp = tables["composition"]
    assert (comp.n_nonmissing_records + comp.n_missing_records).eq(comp.n_records).all()
    assert (comp.n_nonmissing_clusters + comp.n_missing_clusters).eq(comp.n_clusters).all()
    for _, block in comp.groupby(["scenario", "group", "size_band", "attribute"], observed=True):
        if block.n_nonmissing_records.iloc[0]:
            assert np.isclose(block.record_proportion.sum(), 1)
            assert np.isclose(block.cluster_mean_proportion.sum(), 1)
            assert block.category_records.sum() == block.n_nonmissing_records.iloc[0]
    ent = tables["entropy"]
    assert (ent.n_nonmissing + ent.n_missing).eq(ent.n_clusters).all()
    if scotland_baseline:
        total = primary.loc[primary.size_band.eq("All")]
        assert total.n_clusters.sum() == 8962
        assert total.n_records.sum() == 152177
        assert total.n_burden_eligible.sum() == 4216
        assert total.set_index("group").loc["CSC", "n_clusters"] == 631
