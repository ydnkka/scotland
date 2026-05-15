# Analysis plan: cluster-size category as a function of demographics, vaccination, and local context

**Question.** Conditional on local genomic surveillance intensity and local epidemic incidence, are particular sociodemographic groups, vaccination states, or reinfection profiles over-represented among sequences that fall into *large* transmission clusters — the category most plausibly enriched for superspreading-driven expansion?

**Outcome.** A three-level categorical label assigned per sequence at a single primary clustering resolution.

**Model.** Random forest classifier (multinomial), with permutation-based importance and SHAP for direction-of-effect. Grouped cross-validation by cluster.

---

## 1. Design choices fixed

| Choice | Setting |
|---|---|
| Primary clustering resolution | `resolution == 0.3` |
| Cluster size categories | `singleton` (size = 1), `small` (2–9), `large` (≥10) |
| Outcome | Multi-class: `{singleton, small, large}` |
| Window de-duplication | One row per sequence: the window where `|collection_date − wn_mid_date|` is smallest |

After filtering, the modelling table is **one row per sequence**, with a single `cluster_size_category` label.

---

## 2. Outcome construction

```
df0 = read parquet
df0 = df0[df0.resolution == 0.3]
df0 = df0[df0.nextclade_qc.isin(["good", "mediocre"])]   # drop bad QC

# pick the "most central" window per sequence
df0["dist_to_mid"] = (df0.collection_date - df0.wn_mid_date).abs()
df = df0.sort_values("dist_to_mid").drop_duplicates("sequence_id", keep="first")

# size category
def cat(n):
    if n == 1: return "singleton"
    if n <= 9: return "small"
    return "large"
df["cluster_size_category"] = df["cluster_size"].map(cat)
```

Sanity-check: tabulate `cluster_size_category` by `who_voc` and by year/quarter to confirm large clusters are not all concentrated in one variant. Class imbalance is expected — large is the smallest class.

---

## 3. Variables to include

Grouped by epistemic role. *Focal* variables answer the research question; *controls* adjust for surveillance and incidence (your stated confounders); *contextual* variables absorb time/variant/geography heterogeneity that would otherwise be soaked up by focal variables.

### 3.1 Focal — sociodemographic (individual)

| Variable | Encoding | Note |
|---|---|---|
| `age_midpoint` | continuous | Prefer over `age_band` factor — RF handles it cleanly and ordinal information is preserved. |
| `is_female` | binary | |
| `dz_simd_quintile` | ordinal 1–5 | Prefer over the 1–6,976 rank: collinear with other ranks, less interpretable. Use the quintile as the headline deprivation variable. |
| `dz_simd_access_rank` | continuous | Captures geographic remoteness — partly orthogonal to overall SIMD. Optional. |
| `dz_simd_housing_rank` | continuous | Housing-domain deprivation is independently associated with crowding and transmission. Optional. |
| `dz_urban_rural_class` | categorical (6) | One-hot encode. Strong potential modifier of cluster size — urban DZs structurally permit larger clusters. |
| `dz_population_density` | continuous, log1p | Heavy-tailed; log-transform. |

Do **not** simultaneously include `dz_simd_rank`, `dz_simd_quintile`, and the individual domain ranks — they are highly collinear and will fragment importance. The recommended minimal SIMD set is `dz_simd_quintile` plus, optionally, `access` and `housing` as theory-motivated domain ranks. The other domain ranks (`income`, `employment`, `education`, `health`, `crime`) should be excluded unless you have a specific hypothesis about them.

### 3.2 Focal — vaccination and infection history

| Variable | Encoding | Note |
|---|---|---|
| `vacc_dose_number` | ordinal (0–4+) | More informative than the binary `is_vaccinated`; subsumes it. |
| `days_since_vaccination` | continuous | NaN for unvaccinated — see §4 on missingness. |
| `vacc_booster` | binary | NaN → 0 for unvaccinated. |
| `vacc_product_name` | categorical | One-hot top-N (Pfizer, Moderna, AstraZeneca, Janssen, Novavax) + "Other/Unknown". |
| `is_reinfection` | binary | Patient previously positive ≥ 90 days earlier. |

Exclude `is_vaccinated` and `vacc_date_prior` — both are functions of variables already included.

### 3.3 Controls — local surveillance and incidence (your stated adjustments)

These are the variables you specifically wanted to adjust for. They must be in the model so that focal effects are interpretable as "over and above" local surveillance and incidence.

| Variable | Encoding | Note |
|---|---|---|
| `dz_cum_prop_sequenced` | continuous | Surveillance intensity in the DZ up to collection date. |
| `dz_cum_incidence_per_capita` | continuous, log1p | Local cumulative epidemic burden. |
| `dz_7d_test_positivity` | continuous | Local *recent* positivity — captures incidence pressure at the time of sampling. |
| `wn_prop_sequenced` | continuous | Window-level surveillance proxy. |
| `wn_no_sequences` | continuous, log1p | Window-level sampling size; larger windows produce larger-cluster opportunities. |

### 3.4 Contextual — variant, calendar time, region

These would otherwise be confounders for almost every focal variable (vaccination programme rolled out over time, variant fitness varied, etc.). Including them prevents focal-variable importances from being inflated.

| Variable | Encoding | Note |
|---|---|---|
| `who_voc` | categorical | One-hot. Use this rather than `pango_lineage` (too many levels for RF; lineage is captured by within-lineage clustering already). Add an "Other/Non-VOC" level for missing. |
| `collection_date` | continuous (days since 2020-01-01) | Calendar-time trend. RF will learn arbitrary functions. |
| `dz_health_board` | categorical (14) | Regional geography. Local authority (32 levels) is finer but invites overfitting; HB is the better default. |
| `hb_hospital_occupancy` | continuous, log1p | Regional epidemic pressure on the matched HB report date. |
| `hb_reinfection_rate` | continuous | Regional reinfection backdrop. |

### 3.5 Test-record covariates (keep, but interpret carefully)

| Variable | Encoding | Note |
|---|---|---|
| `test_reason` | categorical | `COMMUNITY`, `HOSPITAL`, `CARE_HOME`, other — strongly associated with large-cluster outcomes for structural reasons (care homes ⇒ outbreak sampling). Including it is necessary; otherwise care-home demographics will be mis-attributed to age/deprivation. |
| `test_type` | categorical | PCR vs ANTIGEN. |

Drop `s_gene_status` — informative only within a narrow window (Alpha-vs-Delta and the early Omicron period); largely captured by `who_voc`.

---

## 4. Variables to *exclude* (and why)

| Variable | Reason |
|---|---|
| `cluster_size`, `cluster_n_datazones`, `cluster_duration_days`, `cluster_start_date`, `cluster_end_date` | **Outcome leakage** — these define or are computed from the outcome. |
| `cluster_id`, `sequence_id`, `patient_id`, `specimen_id` | Identifiers, not features. (`cluster_id` and `patient_id` are used for grouping in CV.) |
| `window_idx`, `window_id`, `wn_start_date`, `wn_mid_date`, `wn_end_date` | Captured by `collection_date`; window-membership leaks the (window × lineage) cluster pool. |
| `dz_simd_rank`, `dz_simd_decile`, `dz_simd_vigintile`, plus most SIMD domain ranks | Collinear with `dz_simd_quintile`. |
| `dz_xcoord`, `dz_ycoord` | Spatial signal is already captured by `dz_health_board`, `dz_urban_rural_class`, `dz_population_density`. Re-introduce only if you want explicit spatial modelling. |
| `dz_local_authority`, `dz_local_authority_code` | Redundant with `dz_health_board` at this granularity; 32 levels invite overfit. |
| `dz_population`, `dz_working_age_population`, `dz_area_km2` | Subsumed by `dz_population_density`. |
| `dz_total_tests`, `dz_positive_tests`, etc. (single-day counts) | Noisy at the daily DZ level; the 7-day positivity and the cumulative measures are better. |
| `dz_cum_vaccinated`, `dz_total_vaccinated`, `dz_cum_prop_vaccinated` | Optional. `dz_cum_prop_vaccinated` is a candidate addition if you want a local vaccination-coverage control; the others are redundant or noisy. Recommended: include `dz_cum_prop_vaccinated`, exclude the rest. |
| `pango_lineage`, `clade` | Subsumed by `who_voc`; too many levels for RF without target encoding. |
| `is_vaccinated`, `vacc_date_prior` | Implied by `vacc_dose_number` and `days_since_vaccination`. |
| `nextclade_qc` | Used as a filter, not a predictor. |
| `s_gene_status` | Substantially redundant with `who_voc`; many NaN. |
| `sex` | Encoded by `is_female`. |
| `age_band` | Encoded by `age_midpoint`. |

So the final candidate set is roughly **24–28 features**.

---

## 5. Missing data

`days_since_vaccination`, `vacc_booster`, `vacc_product_name` are NaN for unvaccinated patients — that's structural, not random. Two clean approaches:

1. **Sentinel + indicator** (recommended for RF). Impute `days_since_vaccination` with a sentinel value (e.g., −1), impute `vacc_booster` with 0, set `vacc_product_name` to `"None"`. Add an explicit `is_unvaccinated` indicator (= 1 − `is_vaccinated`). The RF will learn the sentinel boundary naturally.
2. **Histogram gradient boosting**: `sklearn.ensemble.HistGradientBoostingClassifier` handles NaN natively and behaves similarly to RF in performance. Worth running as a sensitivity model.

For genuine MAR cases (`dz_7d_test_positivity` early in the pandemic, `hb_*` for DZs missing an HB code), use median imputation with an explicit `*_missing` indicator.

---

## 6. Cross-validation, grouping, and class imbalance

**Grouping is critical.** Sequences in the same `cluster_id` share the outcome label by construction. A standard random split will let the model memorise cluster membership via lineage × time × geography. Use:

- **GroupKFold (5 folds) with `groups = cluster_id`**, so an entire cluster is in train or in test, never both.
- Additionally stratify the *outer* evaluation by `who_voc` × calendar quarter to ensure each variant/period is represented in both train and test.

Class imbalance: use `class_weight="balanced_subsample"` and report macro-averaged metrics (macro-F1, macro-AUC one-vs-rest). Accuracy is uninformative here.

---

## 7. Model configuration

Primary model:

```python
from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(
    n_estimators=1000,
    max_features="sqrt",
    min_samples_leaf=20,           # guard against overfitting on cluster-shared rows
    class_weight="balanced_subsample",
    n_jobs=-1,
    random_state=42,
)
```

Tune `min_samples_leaf ∈ {10, 20, 50, 100}` and `max_depth ∈ {None, 20, 40}` via 5-fold `GroupKFold` CV on macro-F1. Keep `n_estimators` at 1000 throughout.

Sensitivity model: `HistGradientBoostingClassifier(loss="log_loss", class_weight="balanced", max_iter=600)`.

Calibration: fit `CalibratedClassifierCV(method="isotonic", cv=GroupKFold(5))` on top of the chosen RF for probability outputs you'll trust for downstream rate comparisons.

---

## 8. Interpretation

Three complementary views, not just one:

1. **Permutation importance** (`sklearn.inspection.permutation_importance` on the held-out folds, with `n_repeats=20`) — robust to mixed-type features and the high-cardinality categorical bias of impurity importance.
2. **SHAP TreeExplainer** — per-class SHAP values to read the *sign* of each variable's contribution to the `large` class specifically. Aggregate over the test set; report mean |SHAP| and SHAP-summary beeswarms restricted to the focal variables.
3. **Partial-dependence / ICE plots** for the focal variables — `age_midpoint`, `dz_simd_quintile`, `vacc_dose_number`, `days_since_vaccination`, `is_reinfection` — with the `large` class as the response. These translate cleanly into the manuscript figures.

A useful summary table for the paper: for each focal variable, report mean |SHAP on large class|, direction (sign of average SHAP at high vs low values), and 95% bootstrap CI.

---

## 9. Sensitivity analyses (planned, not optional)

1. Re-fit at `resolution ∈ {0.1, 0.2, 0.4, 0.5}` and confirm the top focal variables remain in the top-10 permutation importances.
2. Re-fit with a stricter large threshold (≥ 20) and confirm direction of focal effects.
3. Re-fit restricted to a single dominant VOC (Delta, then Omicron BA.1, then BA.2) to check whether sociodemographic associations are variant-specific.
4. Patient-level re-aggregation: collapse to one row per `patient_id` (most recent sequence) and re-fit — confirms results aren't driven by repeat-sequenced patients.
5. Drop `test_reason` and re-fit — show how much of the deprivation/age signal is mediated through care-home/hospital testing.

---

## 10. Pre-modelling diagnostics worth running first

Before fitting, generate these — they often reshape the variable list:

- Outcome class counts overall and per VOC and per quarter.
- Pairwise Spearman correlation matrix among the SIMD domain ranks (confirm the multicollinearity argument for using only the quintile).
- Correlation of `dz_cum_prop_sequenced` with `dz_cum_incidence_per_capita` — if very high, drop one.
- Within-cluster homogeneity of focal variables: for each cluster of size ≥ 10, compute the dominant level of `dz_health_board`, `age_band`, `test_reason`. If clusters are nearly monolithic on a variable, that variable's importance will be driven by between-cluster differences (which is the right thing) — but worth quantifying.

---

## 11. What this design will and won't show

It **will** identify sociodemographic and vaccination strata over-represented among sequences in large clusters, adjusting for local surveillance and incidence, variant, calendar time, and region.

It **will not** establish that any of those strata *cause* superspreading. Large clusters are over-detected in settings with intense local sampling and in long-stay congregate settings; both are partly controlled here, but residual confounding (e.g., unmeasured contact-rate heterogeneity) remains. Frame the manuscript in terms of *association with detected large-cluster membership*, not causation.

A complementary cluster-level model — one row per cluster, predicting `cluster_size` from cluster-aggregated covariates — would strengthen the inference. Consider it as a parallel analysis.
