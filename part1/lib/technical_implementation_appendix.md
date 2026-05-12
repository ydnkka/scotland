# Supplementary Methods Appendix 1  
# Technical implementation of the regression analysis

This appendix contains two parts. **Part A** is a concise manuscript-ready Methods section that can be copied into the main manuscript. **Part B** gives the detailed implementation notes, including equations and code excerpts corresponding to the Python analysis pipeline.

---

# Part A. Concise manuscript Methods text

## Study population and source data

We analysed SARS-CoV-2 whole-genome sequences collected from Scottish residents between July 2020 and February 2023 through Public Health Scotland and the COVID-19 Genomics UK (COG-UK) surveillance infrastructure. Sequences were aligned to the SARS-CoV-2 reference genome and quality-controlled using Nextclade; only sequences passing the `good` overall quality filter were retained. Records were linked to sequence metadata, including collection date, anonymised patient identifier, datazone of residence, age band, sex, and Pango lineage, and to the Scottish Index of Multiple Deprivation 2020 version 2 (SIMD), including overall rank and seven domain ranks: income, employment, education, health, geographic access, crime, and housing. Datazone-level testing aggregates and health-board-level surveillance trends were used to derive local surveillance covariates.

## Cluster inference pipeline

Cluster inference followed the EpiLink-Leiden framework described elsewhere. Briefly, sequences were grouped into 3-week sliding windows advanced in 1-week steps and then partitioned by Pango lineage. Within each window-lineage group, pairwise genetic distances were computed under the TN93 substitution model and paired with specimen-collection time differences. Each pair was scored using EpiLink, a process-based compatibility model that evaluates whether the observed genetic and temporal separation is consistent with simulated recent-transmission scenarios under specified assumptions about mutation accumulation, disease progression, and infectiousness. Pairwise compatibility scores were used as weighted graph edges, and recent-transmission clusters were inferred using the Leiden community-detection algorithm at the primary resolution $\gamma = 0.3$. Sequences in groups too small for pairwise processing were assigned to singleton clusters.

## Cluster-level outcomes and covariates

The long-format cluster table was collapsed to one row per inferred cluster at the primary Leiden resolution. For cluster $c$, let $N_c$ denote the number of distinct sequences and $D_c$ the number of distinct residential datazones represented in the cluster. The primary outcomes were cluster size $N_c$ and geographic dispersion $D_c$. Cluster duration was calculated as the interval between the earliest and latest sampled case but retained only descriptively because the fixed 3-week window structure mechanically constrains observed duration.

Because both primary count outcomes are bounded below by 1, we modelled excess counts above the structural minimum:

$$
Y^{(N)}_c=N_c-1,
\qquad
Y^{(D)}_c=D_c-1.
$$

The associated hurdle indicators were:

$$
H^{(N)}_c=\mathbb{I}(N_c > 1),
\qquad
H^{(D)}_c=\mathbb{I}(D_c > 1).
$$

Mean cluster-level deprivation was defined as the mean SIMD rank across sampled cases in the cluster, multiplied by $-1$ so that higher values indicate greater deprivation. Local cumulative incidence was log-transformed after scaling per 1,000 population. Local cumulative sequencing fraction, window-level sequencing proportion, and local 7-day test positivity were logit-transformed after clipping probabilities away from 0 and 1. Continuous covariates were standardised to mean 0 and unit variance. Calendar time was represented using a cubic B-spline basis over the window index, and Pango lineages represented by fewer than 50 clusters were pooled into an "Other rare lineages" category.

## Within-cluster sociodemographic mixing

Within-cluster sociodemographic mixing was measured using observed-minus-expected pairwise discordance. For categorical variable $V$, let $n_{ck}$ be the number of valid sampled cases in cluster $c$ belonging to category $k$, and let $n_c=\sum_k n_{ck}$. Observed pairwise discordance was defined as:

$$
d^{\mathrm{obs}}_{c,V}=1-\sum_k\frac{n_{ck}(n_{ck}-1)}{n_c(n_c-1)}.
$$

Expected discordance was computed from all sampled cases in the same window-lineage stratum. If $m_{slk}$ is the number of cases in window $s$, lineage $l$, and category $k$, with $m_{sl=\sum_k m_{slk}}$, then:

$$
d^{\mathrm{exp}}_{sl,V}=1-\sum_k\frac{m_{slk}(m_{slk}-1)}{m_{sl}(m_{sl}-1)}.
$$

The cluster-level excess-mixing outcome was:

$$
e_{c,V} =d^{\mathrm{obs}}_{c,V}-d^{\mathrm{exp}}_{s(c)l(c),V}.
$$

Positive values indicate more cross-category mixing than expected for the cluster's window-lineage context, whereas negative values indicate more assortative composition. Excess discordance was computed for SIMD quintile, age band, sex, and a joint SIMD-age-sex profile.

## Regression analyses

The primary analysis examined deprivation as the exposure. For each count outcome, we fitted a two-part hurdle model. The first component modelled whether the outcome exceeded its structural minimum:

$$
H_c \sim \operatorname{Bernoulli}(\pi_c),
\qquad
\operatorname{logit}(\pi_c)=\mathbf{x}_c^\top\boldsymbol{\gamma}.
$$

The second component modelled the positive excess count using a zero-truncated negative-binomial model:

$$
Y_c \mid Y_c > 0
\sim
\operatorname{ZTNB}(\mu_c,\alpha),
\qquad
\log(\mu_c)=\mathbf{x}_c^\top\boldsymbol{\beta}.
$$

Here, $\mathbf{x}_c$ included standardised deprivation, local cumulative incidence, local cumulative sequencing fraction, window-level sequencing proportion, local test positivity, calendar spline terms, and lineage fixed effects. For geographic dispersion, an additional pre-specified size-adjusted model included log cluster size as a covariate.

The ZTNB component used the negative-binomial variance function:

$$
\operatorname{Var}(Y_c)=\mu_c + \alpha\mu_c^2,
$$

with $\alpha>0$ estimated on the log scale. If $f_{\mathrm{NB}}(y;\mu,\alpha)$ is the ordinary negative-binomial probability mass function and $p_0=f_{\mathrm{NB}}(0;\mu,\alpha)$, then the zero-truncated likelihood contribution for $y=1,2,\ldots$ is:

$$
P(Y=y \mid Y>0)=\frac{f_{\mathrm{NB}}(y;\mu,\alpha)}{1-p_0}.
$$

For the excess-mixing outcomes, we fitted linear models among non-singleton clusters:

$$
e_{c,V}=\mathbf{x}_c^\top \boldsymbol{\delta}+\varepsilon_c,
$$

where $\mathbf{x}_c$ included the same deprivation and surveillance covariates, calendar spline terms, lineage fixed effects, and log cluster size. Coefficients were reported as percentage-point differences in excess discordance per 1 standard deviation higher covariate value.

To address the complementary question of whether bridging across sociodemographic strata predicted cluster scale, we refitted the count models with four standardised excess-mixing predictors: SIMD quintile, age band, sex, and joint SIMD-age-sex profile. These models retained deprivation, surveillance covariates, calendar splines, and lineage fixed effects as adjustments. Because excess mixing is undefined for singleton clusters, the cluster-size hurdle component was not estimable in this analysis; the positive cluster-size ZTNB component and both geographic-dispersion components were fitted among non-singleton clusters.

## Inference

All models used cluster-robust standard errors, clustered by window in the primary analysis to account for dependence induced by overlapping windows. A pre-specified sensitivity analysis clustered standard errors by health board. ZTNB standard errors were obtained using a sandwich covariance estimator with a numerical-Hessian bread and score-based meat. Let $\hat{\boldsymbol{\theta}}$ denote the fitted parameter vector and $s_c(\hat{\boldsymbol{\theta}})$ the observation-level score. For robust cluster $g$, the cluster-level score was:

$$
S_g(\hat{\boldsymbol{\theta}})=\sum_{c\in g}s_c(\hat{\boldsymbol{\theta}}).
$$

The sandwich covariance estimator was:

$$
\widehat{\operatorname{Var}}(\hat{\boldsymbol{\theta}}) =
\widehat{B}^{-1}
\left[\sum_g S_g(\hat{\boldsymbol{\theta}})S_g(\hat{\boldsymbol{\theta}})^\top\right]
\widehat{B}^{-1},
$$

where $\widehat{B}$ is the observed information matrix. A finite-sample correction,

$$
\frac{G}{G-1}\frac{n-1}{n-p},
$$

was applied when $G>1$ and $n>p$. Hurdle coefficients are reported as adjusted odds ratios, ZTNB coefficients as adjusted count ratios, and linear mixing coefficients as adjusted percentage-point differences.

## Sensitivity and extension analyses

We conducted pre-specified sensitivity analyses using health-board-clustered standard errors, a log offset for the number of sequences in the analysis window in the positive cluster-size model, size adjustment for positive geographic dispersion, index-case SIMD deprivation, 99th-percentile winsorisation of positive counts, an approximately non-overlapping window subset retaining every third window, and log-linear single-component models as comparators. We also fitted wave-stratified models for B.1.177, Alpha, Delta, BA.1, BA.2, BA.4, BA.5, and BQ.1, excluding XBB from regression models because of small sample size. Finally, domain-specific models were fitted for all seven SIMD domains and extended to the mixing-predictor count-model framework.

---

# Part B. Technical implementation details

## 1. Overview of the Python workflow

The implementation is organised around three layers:

1. **Shared constants and model specifications**  
   Defines the primary covariates, mixing variables, count-model outcomes, SIMD domains, wave labels, and run defaults.

2. **Data preparation**  
   Loads sequence-level rows, constructs the cluster-level table, calculates observed and expected discordance, applies covariate transformations, standardises variables, pools sparse lineages, and attaches calendar spline terms.

3. **Model fitting**  
   Builds design matrices, fits binomial hurdle models, implements the ZTNB likelihood and gradient, calculates cluster-robust standard errors, and orchestrates the primary, sensitivity, domain-specific, wave-specific, and mixing-predictor analyses.

The high-level mapping is:

| Analysis question | Main functions |
|---|---|
| Deprivation as exposure for count outcomes | `fit_count_models` |
| Deprivation as exposure for mixing outcomes | `fit_mixing_models` |
| Excess mixing as predictor for count outcomes | `fit_mixing_predictor_count_models` |
| SIMD domain analyses | `fit_domain_count_models`, `fit_domain_*_mixing_models` |
| Wave-stratified analyses | `fit_wave_outcome_models` |
| Log-linear sensitivity analyses | `fit_loglinear_models` |

---

## 2. Shared constants and model specifications

The primary run defaults are:

```python
QC_DEFAULT = "good"
PRIMARY_RESOLUTION = 0.3
LINEAGE_MIN_CLUSTERS = 50
CALENDAR_SPLINE_DF = 8
```

The primary adjustment set is:

```python
PRIMARY_TERMS = [
    "deprivation_z",
    "local_incidence_z",
    "local_seq_fraction_z",
    "window_seq_fraction_z",
    "test_positivity_z",
]
```

The four main excess-mixing variables are:

```python
MIXING_VARIABLES = {
    "simd": {
        "column": "dz_simd_quintile",
        "label": "SIMD quintile mixing",
        "short_label": "SIMD",
    },
    "age": {
        "column": "age_band",
        "label": "Age-band mixing",
        "short_label": "Age",
    },
    "sex": {
        "column": "sex",
        "label": "Sex mixing",
        "short_label": "Sex",
    },
    "profile": {
        "column": "socio_demographic_profile",
        "label": "Joint SIMD-age-sex profile mixing",
        "short_label": "Joint profile",
    },
}
```

The mixing-predictor count models use the corresponding standardised predictors:

```python
MIXING_PREDICTOR_TERMS = [
    f"{prefix}_excess_mixing_z" for prefix in MIXING_VARIABLES
]
```

The count-model specifications define the raw count, binary hurdle outcome, positive excess-count outcome, and component label:

```python
COUNT_MODEL_SPECS = [
    CountModelSpec(
        name="cluster_size",
        label="Cluster size",
        raw_outcome="cluster_size",
        binary_col="cluster_size_gt1",
        positive_col="cluster_size_excess",
        positive_label="Additional sequences among non-singleton clusters",
    ),
    CountModelSpec(
        name="geographic_dispersion",
        label="Geographic dispersion",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
    ),
    CountModelSpec(
        name="geographic_dispersion_size_adjusted",
        label="Geographic dispersion, size-adjusted",
        raw_outcome="cluster_n_datazones",
        binary_col="datazones_gt1",
        positive_col="datazones_excess",
        positive_label="Additional datazones among multi-datazone clusters",
        include_size=True,
    ),
]
```

Duration is not included in `COUNT_MODEL_SPECS`, because the fixed 3-week clustering windows constrain observed cluster span by design.

---

## 3. Sequence-row loading and profile construction

The primary sequence-level input contains cluster identifiers, sequence identifiers, window identifiers, lineage, collection date, datazone, age band, sex, SIMD rank/quintile, testing covariates, and health-board code. The primary loading function converts relevant variables to categorical type and parses dates:

```python
seq = load_analysis_columns_pandas(
    columns=SEQUENCE_COLUMNS,
    resolution=primary_resolution,
    qc=qc,
)

categorical = [
    "cluster_id", "sequence_id", "window_id", "datazone", "pango_lineage",
    "nextclade_qc", "age_band", "sex", "dz_simd_quintile",
]
for col in categorical:
    seq[col] = seq[col].astype("category")

seq["collection_date"] = pd.to_datetime(seq["collection_date"])
seq["wn_mid_date"] = pd.to_datetime(seq["wn_mid_date"])
```

The joint SIMD-age-sex profile is constructed only for rows with complete values:

```python
complete_profile = seq[["dz_simd_quintile", "age_band", "sex"]].notna().all(axis=1)
seq["socio_demographic_profile"] = pd.NA
seq.loc[complete_profile, "socio_demographic_profile"] = (
    seq.loc[complete_profile, "dz_simd_quintile"].astype(str)
    + "|"
    + seq.loc[complete_profile, "age_band"].astype(str)
    + "|"
    + seq.loc[complete_profile, "sex"].astype(str)
)
seq["socio_demographic_profile"] = seq["socio_demographic_profile"].astype("category")
```

For domain and wave analyses, a related loader derives domain-specific SIMD quintiles and assigns broad epidemic wave labels from Pango lineage.

---

## 4. Cluster table construction

The cluster-level analysis table is built by grouping sequence rows by `cluster_id`:

```python
clusters = (
    seq.groupby("cluster_id", observed=True, sort=False)
    .agg(
        cluster_size=("sequence_id", "nunique"),
        cluster_n_datazones=("datazone", "nunique"),
        cluster_start_date=("collection_date", "min"),
        cluster_end_date=("collection_date", "max"),
        resolution=("resolution", "first"),
        window_id=("window_id", "first"),
        window_idx=("window_idx", "first"),
        wn_mid_date=("wn_mid_date", "first"),
        pango_lineage=("pango_lineage", "first"),
        mean_simd_rank=("dz_simd_rank", "mean"),
        mean_local_incidence_per_capita=("dz_cum_incidence_per_capita", "mean"),
        mean_local_seq_fraction=("dz_cum_prop_sequenced", "mean"),
        mean_window_seq_fraction=("wn_prop_sequenced", "mean"),
        mean_test_positivity=("dz_7d_test_positivity", "mean"),
        wn_no_sequences=("wn_no_sequences", "first"),
        health_board=("dz_health_board_code", "first"),
    )
    .reset_index()
)
```

The index-case SIMD exposure is obtained from the earliest sampled sequence in each cluster:

```python
_seq_sorted = seq.sort_values("collection_date")
_index_simd = (
    _seq_sorted.groupby("cluster_id", observed=True)["dz_simd_rank"]
    .first()
    .rename("index_simd_rank")
)
clusters = clusters.merge(_index_simd, on="cluster_id", how="left")
```

The structural hurdle and positive-count variables are then created:

```python
clusters["duration_days"] = (
    clusters["cluster_end_date"] - clusters["cluster_start_date"]
).dt.days.astype(int)

clusters["cluster_size_gt1"] = (clusters["cluster_size"] > 1).astype(int)
clusters["duration_gt0"] = (clusters["duration_days"] > 0).astype(int)
clusters["datazones_gt1"] = (clusters["cluster_n_datazones"] > 1).astype(int)

clusters["cluster_size_excess"] = clusters["cluster_size"] - 1
clusters["duration_positive_days"] = clusters["duration_days"]
clusters["datazones_excess"] = clusters["cluster_n_datazones"] - 1
```

Thus, the ZTNB positive count is not the raw count itself but the count above the structural minimum.

---

## 5. Pairwise discordance and excess mixing

### 5.1 Pairwise discordance from category counts

For a categorical variable $V$, observed or expected discordance is calculated from category counts. The implementation first computes the number of valid observations and the number of ordered same-category pairs:

```python
totals = counts.groupby(group_cols, observed=True)["n"].sum().rename("n_valid")
same_pairs = (
    counts.assign(same_pairs=counts["n"] * (counts["n"] - 1))
    .groupby(group_cols, observed=True)["same_pairs"]
    .sum()
)
```

The denominator is the number of ordered pairs:

```python
denom = out["n_valid"] * (out["n_valid"] - 1)
```

and discordance is:

```python
out["discordance"] = np.nan
mask = denom > 0
out.loc[mask, "discordance"] = (
    1 - out.loc[mask, "same_pairs"] / denom.loc[mask]
)
```

This corresponds to:

$$
d=1-\frac{\sum_k n_k(n_k-1)}{n(n-1)}.
$$

The use of $n(n-1)$, rather than $n(n-1)/2$, is equivalent because the same ordered-pair convention is used in both numerator and denominator.

### 5.2 Observed cluster discordance

Observed discordance is computed within clusters:

```python
def observed_cluster_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    counts = (
        seq.dropna(subset=[variable])
        .groupby(["cluster_id", variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = _pairwise_discordance_from_counts(counts, ["cluster_id"])
    return out.rename(
        columns={
            "n_valid": f"{prefix}_n_valid",
            "discordance": f"{prefix}_discordance",
        }
    )
```

### 5.3 Expected window-lineage discordance

Expected discordance is computed using the sampled composition of the corresponding window-lineage stratum:

```python
def expected_stratum_discordance(
    seq: pd.DataFrame,
    variable: str,
    prefix: str,
) -> pd.DataFrame:
    strata = ["window_id", "pango_lineage"]
    counts = (
        seq.dropna(subset=[variable])
        .groupby(strata + [variable], observed=True)
        .size()
        .rename("n")
        .reset_index()
    )
    out = _pairwise_discordance_from_counts(counts, strata)
    return out.rename(
        columns={
            "n_valid": f"{prefix}_stratum_n_valid",
            "discordance": f"{prefix}_expected_discordance",
        }
    )
```

### 5.4 Excess discordance

For each mixing variable, observed and expected discordance are merged into the cluster table, and the excess value is calculated as:

```python
for prefix, spec in MIXING_VARIABLES.items():
    obs = observed_cluster_discordance(seq, spec["column"], prefix)
    exp = expected_stratum_discordance(seq, spec["column"], prefix)
    clusters = clusters.merge(obs, on="cluster_id", how="left")
    clusters = clusters.merge(exp, on=["window_id", "pango_lineage"], how="left")
    clusters[f"{prefix}_excess_discordance"] = (
        clusters[f"{prefix}_discordance"]
        - clusters[f"{prefix}_expected_discordance"]
    )
```

This gives:

$$
e_{c,V}=d^{\mathrm{obs}}_{c,V}-d^{\mathrm{exp}}_{s(c)l(c),V}.
$$

---

## 6. Covariate transformations and standardisation

The primary deprivation exposure is constructed by negating mean SIMD rank:

```python
clusters["deprivation_raw"] = -clusters["mean_simd_rank"]
clusters["index_deprivation_raw"] = -clusters["index_simd_rank"]
```

This ensures that larger values correspond to greater deprivation.

The surveillance covariates are transformed as follows:

```python
clusters["local_incidence_log"] = np.log1p(
    clusters["mean_local_incidence_per_capita"].clip(lower=0) * 1000
)
clusters["local_seq_fraction_logit"] = logit_clipped(clusters["mean_local_seq_fraction"])
clusters["window_seq_fraction_logit"] = logit_clipped(clusters["mean_window_seq_fraction"])
clusters["test_positivity_logit"] = logit_clipped(clusters["mean_test_positivity"].fillna(0))
clusters["log_cluster_size"] = np.log(clusters["cluster_size"])
```

The logit transform clips values to avoid infinite logits:

```python
def logit_clipped(values: pd.Series, eps: float = 1e-5) -> pd.Series:
    clipped = values.clip(lower=eps, upper=1 - eps)
    return np.log(clipped / (1 - clipped))
```

Standardisation uses the population standard deviation:

```python
def zscore(values: pd.Series) -> tuple[pd.Series, float, float]:
    mean = float(values.mean())
    sd = float(values.std(ddof=0))
    if not math.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardise {values.name!r}: zero or invalid SD.")
    return (values - mean) / sd, mean, sd
```

The main transform dictionary is:

```python
transforms = {
    "deprivation_z": "deprivation_raw",
    "index_deprivation_z": "index_deprivation_raw",
    "local_incidence_z": "local_incidence_log",
    "local_seq_fraction_z": "local_seq_fraction_logit",
    "window_seq_fraction_z": "window_seq_fraction_logit",
    "test_positivity_z": "test_positivity_logit",
    "log_cluster_size_z": "log_cluster_size",
}

transforms.update(
    {
        f"{prefix}_excess_mixing_z": f"{prefix}_excess_discordance"
        for prefix in MIXING_VARIABLES
    }
)
```

The scaling metadata are stored for reproducibility:

```python
for z_col, raw_col in transforms.items():
    clusters[z_col], mean, sd = zscore(clusters[raw_col])
    scaling_rows.append(
        {
            "standardised_column": z_col,
            "source_column": raw_col,
            "source_mean": mean,
            "source_sd": sd,
        }
    )
```

---

## 7. Calendar splines and lineage pooling

Sparse lineages are pooled before modelling:

```python
lineage_counts = clusters["pango_lineage"].astype(str).value_counts()
common_lineages = set(lineage_counts[lineage_counts >= lineage_min_clusters].index)

clusters["lineage_model"] = np.where(
    clusters["pango_lineage"].astype(str).isin(common_lineages),
    clusters["pango_lineage"].astype(str),
    "Other rare lineages",
)
```

Calendar time is represented using a cubic B-spline over `window_idx`:

```python
calendar = dmatrix(
    f"bs(window_idx, df={calendar_spline_df}, degree=3, include_intercept=False) - 1",
    clusters,
    return_type="dataframe",
)
calendar.columns = [f"calendar_spline_{i + 1}" for i in range(calendar.shape[1])]
```

In the primary analysis, `calendar_spline_df = 8`.

---

## 8. Design-matrix construction

Design matrices are constructed in a fixed order:

1. intercept;
2. numeric covariates;
3. calendar spline columns;
4. lineage dummy variables.

```python
parts = [
    pd.DataFrame({"const": np.ones(len(df), dtype=float)}, index=df.index),
    df[numeric_terms].astype(float),
    df[calendar_cols].astype(float),
]
```

Lineage dummies are constructed from a categorical variable with a common set of levels:

```python
lineages = pd.Categorical(
    df["lineage_model"].astype(str),
    categories=all_lineage_levels,
    ordered=False,
)

lineage_dummies = pd.get_dummies(
    pd.Series(lineages, index=df.index, name="lineage_model"),
    prefix="lineage",
    drop_first=True,
    dtype=float,
)
parts.append(lineage_dummies)
```

Empty lineage columns are pruned:

```python
zero_columns = [
    col
    for col in x.columns
    if col not in {"const", *numeric_terms} and float(x[col].abs().sum()) == 0
]
if zero_columns:
    x = x.drop(columns=zero_columns)
```

For wave-stratified models, a further rank check removes redundant columns by sequential orthogonalisation. The column order preserves substantive covariates before calendar and lineage terms wherever possible.

---

## 9. Hurdle model implementation

For a count outcome with structural minimum $a=1$, the analysis uses:

$$
H_c = \mathbb{I}(Y^{\mathrm{raw}}_c > a),
\qquad
Y_c = Y^{\mathrm{raw}}_c - a.
$$

The binary component is fitted with a binomial GLM:

```python
model = sm.GLM(y, x, family=sm.families.Binomial())
result = model.fit(
    maxiter=maxiter,
    cov_type="cluster",
    cov_kwds={"groups": groups},
)
```

The component-specific design matrix is built from:

```python
terms = _model_terms(
    spec,
    primary_terms=primary_terms,
    extra_terms=extra_terms,
)
```

where `_model_terms` appends any extra predictors and, for size-adjusted geographic dispersion, appends `log_cluster_size_z`:

```python
def _model_terms(spec, primary_terms=PRIMARY_TERMS, extra_terms=None):
    terms = list(primary_terms)
    if extra_terms:
        terms.extend(extra_terms)
    if spec.include_size:
        terms.append("log_cluster_size_z")
    return terms
```

Exponentiated logit coefficients are reported as odds ratios:

$$
\operatorname{aOR}_j=\exp(\hat{\gamma}_j).
$$

---

## 10. ZTNB likelihood and optimisation

### 10.1 Parameterisation

The ZTNB parameter vector is:

$$
\boldsymbol{\theta} = (\boldsymbol{\beta}, \log\alpha).
$$

The code transforms:

```python
beta = params[:-1]
log_alpha = float(np.clip(params[-1], -10.0, 8.0))
alpha = float(np.exp(log_alpha))
r = 1.0 / alpha
```

The mean model is:

```python
eta = x @ beta
if offset is not None:
    eta = eta + offset

eta = np.clip(eta, -30.0, 30.0)
mu = np.exp(eta)
```

Thus:

$$
\log(\mu_i)=\eta_i=\mathbf{x}_i^\top\boldsymbol{\beta}+\mathrm{offset}_i.
$$

### 10.2 Negative-binomial log-PMF

With $r=1/\alpha$, the ordinary negative-binomial log-PMF is:

$$
\log f_{\mathrm{NB}}(y_i;\mu_i,\alpha) = 
\log\Gamma(y_i+r) - 
\log\Gamma(r) - 
\log\Gamma(y_i+1) +
r\{\log r-\log(r+\mu_i)\} +
y_i\{\log\mu_i-\log(r+\mu_i)\}.
$$

The implementation is:

```python
log_r = math.log(r)
log_r_mu = np.log(r + mu)

logpmf = (
    gammaln(y + r)
    - gammaln(r)
    - gammaln(y + 1)
    + r * (log_r - log_r_mu)
    + y * (eta - log_r_mu)
)
```

### 10.3 Zero-truncation adjustment

The probability of zero under the untruncated negative binomial is:

$$
p_{0i}=\left(\frac{r}{r+\mu_i}\right)^r.
$$

The code computes:

```python
log_p0 = r * (log_r - log_r_mu)
p0 = np.exp(np.clip(log_p0, -745.0, -1e-12))
one_minus_p0 = np.clip(1.0 - p0, 1e-12, 1.0)
log_one_minus_p0 = np.log(one_minus_p0)
```

The zero-truncated log-likelihood contribution is:

$$
\ell_i= \log f_{\mathrm{NB}}(y_i;\mu_i,\alpha) - \log(1-p_{0i}).
$$

The implementation is:

```python
ll_obs = logpmf - log_one_minus_p0
llf = float(np.sum(ll_obs))
```

### 10.4 Analytical gradient

The implementation returns both summed and observation-level score contributions:

```python
score_obs = np.column_stack([score_beta, score_log_alpha])
score = np.sum(score_obs, axis=0)
return llf, score, score_obs
```

For the regression coefficients:

```python
score_eta = (
    r * (y - mu) / (r + mu)
    - mu * r * p0_ratio / (r + mu)
)
score_beta = x * score_eta[:, None]
```

For the dispersion parameter, the implementation first obtains the score with respect to $r$, then applies the chain rule for $\log\alpha$:

$$
\frac{\partial r}{\partial\log\alpha}=-r.
$$

```python
score_r = dlogpmf_dr + p0_ratio * dlogp0_dr
score_log_alpha = -r * score_r
```

### 10.5 Optimisation

The optimiser minimises the negative log-likelihood:

```python
def _ztnb_objective(params, y, x, offset=None):
    llf, score, _ = ztnb_loglike_score(params, y, x, offset)
    if not np.isfinite(llf) or not np.all(np.isfinite(score)):
        return 1e100, np.zeros_like(params)
    return -llf, -score
```

Starting values use a Poisson GLM when possible and a method-of-moments initial value for $\alpha$:

```python
beta = np.zeros(x.shape[1], dtype=float)
beta[0] = float(np.log(np.mean(y)))

try:
    poisson = sm.GLM(y, x, family=sm.families.Poisson()).fit(maxiter=50, disp=0)
    beta = np.asarray(poisson.params, dtype=float)
except Exception:
    pass

mean = float(np.mean(y))
var = float(np.var(y, ddof=1))
alpha = max((var - mean) / (mean * mean), 0.05)

return np.r_[beta, math.log(alpha)]
```

The model is fitted using L-BFGS-B with bounds on $\log\alpha$:

```python
opt = minimize(
    _ztnb_objective,
    start,
    args=(y, x_array, offset),
    method="L-BFGS-B",
    jac=True,
    bounds=[(None, None)] * x_array.shape[1] + [(-10.0, 8.0)],
    options={"maxiter": maxiter, "ftol": 1e-8, "gtol": 1e-5, "maxls": 50},
)
```

---

## 11. Cluster-robust covariance estimation

### 11.1 Hessian bread

The observed Hessian is calculated by finite differences of the analytical gradient:

```python
_, grad0, _ = ztnb_loglike_score(params, y, x, offset)
H = np.zeros((n, n))

for i in range(n):
    p_fwd = params.copy()
    p_fwd[i] += eps
    _, grad_fwd, _ = ztnb_loglike_score(p_fwd, y, x, offset)
    H[i] = (grad_fwd - grad0) / eps

return (H + H.T) / 2.0
```

The observed information matrix is:

$$
\widehat{B}=-H.
$$

The implementation uses a Hermitian pseudo-inverse:

```python
info = -H
bread_inv = pinvh(info, rtol=1e-10)
```

### 11.2 Score meat

Observation-level scores are summed within each robust cluster:

```python
group_codes, inverse = np.unique(groups, return_inverse=True)
cluster_scores = np.zeros((len(group_codes), len(params)), dtype=float)

for group_idx in range(len(group_codes)):
    cluster_scores[group_idx, :] = score_obs[inverse == group_idx].sum(axis=0)
```

The meat is:

```python
meat = cluster_scores.T @ cluster_scores
```

corresponding to:

$$
\widehat{M}=\sum_g S_g S_g^\top.
$$

### 11.3 Sandwich covariance

The robust covariance is:

```python
cov = bread_inv @ meat @ bread_inv
```

A finite-sample correction is applied when there is more than one robust cluster and $n>p$:

```python
n, p = x_array.shape
if len(group_codes) > 1 and n > p:
    correction = (len(group_codes) / (len(group_codes) - 1)) * ((n - 1) / (n - p))
    cov *= correction
```

Standard errors and Wald p-values are then:

```python
bse = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
z_values = np.divide(params, bse, out=np.full_like(params, np.nan), where=bse > 0)
pvalues = 2 * norm.sf(np.abs(z_values))
```

---

## 12. Cluster-robust binomial fallback path

For most binomial hurdle models, `statsmodels` clustered covariance is used directly. For wave-stratified models, the standard covariance path can fail when the Hessian is singular or near-singular. The implementation therefore includes a direct sandwich estimator.

Predicted probabilities are clipped:

```python
mu = expit(x_array @ params)
mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
```

The bread is based on the logistic information matrix:

```python
weights = mu * (1.0 - mu)
information = (x_array * weights[:, None]).T @ x_array
bread_inv = pinvh(information, rtol=1e-10)
```

The score contributions are:

```python
score_obs = x_array * (y_array - mu)[:, None]
```

and the same group-summed sandwich estimator is used.

The stable binomial log-likelihood clips fitted probabilities further to avoid $\log(0)$:

```python
mu = expit(x_array @ params)
mu = np.clip(mu, 1e-12, 1.0 - 1e-12)

llf = float(np.sum(y_array * np.log(mu) + (1.0 - y_array) * np.log1p(-mu)))
aic = float(-2.0 * llf + 2.0 * len(params))
```

---

## 13. Primary count models

The primary count models are fitted by iterating over `COUNT_MODEL_SPECS`.

For each outcome except the size-adjusted geographic-dispersion specification, the binary hurdle component is fitted:

```python
if not spec.include_size:
    rows, diag = fit_binary_component(
        clusters, spec, lineage_levels_all, calendar_cols, maxiter,
        cluster_by=cluster_by, primary_terms=primary_terms,
    )
```

The positive ZTNB component is fitted for all count specifications:

```python
rows, diag = fit_positive_component(
    clusters, spec, lineage_levels_all, calendar_cols, maxiter,
    cluster_by=cluster_by,
    use_size_offset=use_size_offset,
    winsorise_quantile=winsorise_quantile,
    primary_terms=primary_terms,
)
```

The positive component restricts to $Y_c>0$:

```python
use = clusters.loc[clusters[spec.positive_col] > 0].dropna(
    subset=[spec.positive_col, *terms, *calendar_cols, "lineage_model"]
)
y = use[spec.positive_col].astype(int).to_numpy()
```

---

## 14. Mixing outcome models

The deprivation-as-exposure mixing models use OLS among non-singleton clusters:

```python
terms = list(primary_terms) + ["log_cluster_size_z"]

for prefix, spec in MIXING_VARIABLES.items():
    outcome = f"{prefix}_excess_discordance"
    use = clusters.loc[clusters["cluster_size"] >= 2].dropna(
        subset=[outcome, *terms, *calendar_cols, "lineage_model"]
    )

    y = use[outcome].astype(float)
    x = build_exog(use, terms, calendar_cols, lineage_levels_all)
    groups = use[cluster_by].astype(str).to_numpy()

    model = sm.OLS(y, x)
    result = model.fit(cov_type="cluster", cov_kwds={"groups": groups})
```

The reported effect scale is percentage points:

```python
"coefficient_percentage_points": coef * 100
"ci_low_percentage_points": (coef - 1.96 * stderr) * 100
"ci_high_percentage_points": (coef + 1.96 * stderr) * 100
```

---

## 15. Mixing-predictor count models

The mixing-predictor models extend the primary count models by adding the four standardised excess-mixing variables:

```python
extra_terms = MIXING_PREDICTOR_TERMS
predictor_set = "primary_plus_mixing"
```

The cluster-size hurdle is skipped because singleton clusters cannot have defined mixing predictors:

```python
if spec.name == "cluster_size":
    diagnostics.append(
        _skipped_mixing_predictor_diag(
            spec,
            "hurdle_binary",
            spec.binary_col,
            (
                "mixing predictors require at least two valid cases, "
                "so the cluster-size hurdle has no singleton comparison group"
            ),
        )
    )
```

For geographic dispersion, both the hurdle and positive-count components are fitted with the extra predictors:

```python
rows, diag = fit_binary_component(
    clusters,
    spec,
    lineage_levels_all,
    calendar_cols,
    maxiter,
    cluster_by=cluster_by,
    primary_terms=primary_terms,
    extra_terms=MIXING_PREDICTOR_TERMS,
    analysis_population_label=analysis_population,
)
```

and:

```python
rows, diag = fit_positive_component(
    clusters,
    spec,
    lineage_levels_all,
    calendar_cols,
    maxiter,
    cluster_by=cluster_by,
    use_size_offset=use_size_offset,
    winsorise_quantile=winsorise_quantile,
    primary_terms=primary_terms,
    extra_terms=MIXING_PREDICTOR_TERMS,
)
```

---

## 16. Sensitivity analyses

### 16.1 Health-board-clustered standard errors

The primary `cluster_by` argument is `"window_id"`. Sensitivity analyses can pass `"health_board"` instead:

```python
fit_count_models(
    clusters,
    lineage_levels_all,
    calendar_cols,
    maxiter=maxiter,
    cluster_by="health_board",
)
```

### 16.2 Window-pool offset for positive cluster size

The offset sensitivity applies only to the cluster-size positive model:

```python
offset = None
if use_size_offset and spec.name == "cluster_size":
    wn_seq = use["wn_no_sequences"].to_numpy(dtype=float)
    offset = np.log(np.clip(wn_seq, 1.0, None))
```

The model becomes:

$$
\log(\mu_c)=\log(W_c)+\mathbf{x}_c^\top\boldsymbol{\beta},
$$

where $W_c$ is the number of sequences in the corresponding analysis window.

### 16.3 Winsorised positive counts

Extreme positive counts can be capped at a pre-specified quantile:

```python
if winsorise_quantile > 0.0:
    winsorise_cap = int(np.quantile(y, winsorise_quantile))
    y = np.minimum(y, winsorise_cap)
    winsorised = True
```

The planned sensitivity uses the 99th percentile.

### 16.4 Size-adjusted geographic dispersion

The size-adjusted geographic-dispersion specification is represented by:

```python
CountModelSpec(
    name="geographic_dispersion_size_adjusted",
    ...
    include_size=True,
)
```

When `include_size=True`, `_model_terms` appends:

```python
"log_cluster_size_z"
```

### 16.5 Index-case SIMD

Index-case deprivation uses:

```python
"index_deprivation_z"
```

instead of:

```python
"deprivation_z"
```

The exposure is derived from the earliest sampled sequence in each cluster.

### 16.6 Approximately non-overlapping windows

The approximately non-overlapping window sensitivity retains windows satisfying:

```python
window_idx % 3 == 0
```

This reduces dependence from overlapping 3-week windows advanced in 1-week steps.

### 16.7 Log-linear comparator models

The log-linear sensitivity fits OLS models to log-transformed raw outcomes:

```python
y_raw = use[spec["source"]].astype(float) + float(spec["log_plus"])
y = np.log(y_raw)

result = sm.OLS(y, x).fit(cov_type="cluster", cov_kwds={"groups": groups})
```

Exponentiated coefficients are interpreted as geometric mean ratios:

```python
"geometric_mean_ratio": float(np.exp(coef))
```

---

## 17. Domain-specific and wave-stratified extensions

### 17.1 Domain-specific deprivation

The domain analysis repeats the deprivation and mixing analyses using the seven SIMD domains:

```python
DOMAINS = {
    "overall": {"label": "Overall", "rank_col": "dz_simd_rank"},
    "income": {"label": "Income", "rank_col": "dz_simd_income_rank"},
    "employment": {"label": "Employment", "rank_col": "dz_simd_employment_rank"},
    "education": {"label": "Education", "rank_col": "dz_simd_education_rank"},
    "health": {"label": "Health", "rank_col": "dz_simd_health_rank"},
    "access": {"label": "Access", "rank_col": "dz_simd_access_rank"},
    "crime": {"label": "Crime", "rank_col": "dz_simd_crime_rank"},
    "housing": {"label": "Housing", "rank_col": "dz_simd_housing_rank"},
}
```

For each domain, raw deprivation is defined by negating the mean domain rank:

```python
raw_values = -clusters[f"{domain}_mean_rank"]
standardised_cols[z_col], mean, sd = zscore(raw_values)
```

### 17.2 Wave assignment

Broad wave labels are assigned from Pango lineages:

```python
def assign_wave(lineage: str) -> str:
    if not isinstance(lineage, str):
        return "Other"
    if lineage.startswith("B.1.177"):
        return "B.1.177"
    if lineage == "B.1.1.7" or lineage.startswith("B.1.1.7."):
        return "Alpha"
    if lineage.startswith("AY.") or lineage == "B.1.617.2":
        return "Delta"
    if lineage.startswith("BA.1"):
        return "BA.1"
    if lineage.startswith("BA.2"):
        return "BA.2"
    if lineage.startswith("BA.4"):
        return "BA.4"
    if lineage.startswith("BA.5") or lineage.startswith("BE."):
        return "BA.5"
    if lineage.startswith("BQ."):
        return "BQ.1"
    if lineage.startswith("XBB"):
        return "XBB"
    return "Other"
```

The pre-specified wave order is:

```python
WAVE_ORDER = [
    "B.1.177",
    "Alpha",
    "Delta",
    "BA.1",
    "BA.2",
    "BA.4",
    "BA.5",
    "BQ.1",
    "XBB",
]
```

XBB is excluded from regression models where sample size is insufficient.

---

## 18. Model outputs and diagnostics

Each model returns two tidy tables:

1. **Results table**  
   One row per coefficient of interest, with:
   - term name;
   - term label;
   - coefficient;
   - clustered standard error;
   - z statistic;
   - p-value;
   - exponentiated ratio and confidence interval, where appropriate.

2. **Diagnostics table**  
   Model-level metadata, including:
   - outcome;
   - model component;
   - response variable;
   - number of observations;
   - number of events for hurdle models;
   - number of features;
   - number of windows;
   - number of lineage terms;
   - convergence status;
   - log-likelihood;
   - AIC;
   - estimated $\alpha$ for ZTNB models;
   - optimiser message;
   - warnings;
   - whether the model used offsets or winsorisation.

For ratio models, rows are created by exponentiating coefficients and Wald intervals:

```python
row.update(
    {
        "coefficient": coef,
        "std_error_clustered_by_window": stderr,
        "z": coef / stderr if stderr > 0 else np.nan,
        "p_value": float(pvalues[i]),
        "ratio": float(np.exp(coef)),
        "ratio_ci_low": float(np.exp(coef - 1.96 * stderr)),
        "ratio_ci_high": float(np.exp(coef + 1.96 * stderr)),
    }
)
```

---

## 19. Interpretation guide

### 19.1 Hurdle component

For a coefficient $\gamma_j$, the reported adjusted odds ratio is:

$$
\exp(\gamma_j).
$$

For cluster size, this is the association with being a non-singleton cluster. For geographic dispersion, it is the association with spanning more than one datazone.

### 19.2 Positive ZTNB component

For a coefficient $\beta_j$, the reported adjusted count ratio is:

$$
\exp(\beta_j).
$$

This is the multiplicative association with the expected positive excess count, conditional on exceeding the structural minimum.

For cluster size, the positive outcome is $N_c-1$, so the coefficient refers to additional sequences among non-singleton clusters. For geographic dispersion, the positive outcome is $D_c-1$, so the coefficient refers to additional datazones among multi-datazone clusters.

### 19.3 Linear excess-mixing models

For a coefficient $\delta_j$, the reported effect is:

$$
100\delta_j.
$$

percentage points of excess discordance per 1 standard deviation higher covariate value.

### 19.4 Mixing-predictor models

For the mixing-predictor count models, coefficients for standardised excess-mixing variables estimate whether clusters with greater-than-expected bridging across categories are larger or more geographically dispersed after adjustment for deprivation, surveillance covariates, calendar time, and lineage.

---

## 20. References

The references below support the main statistical and computational methods used in the manuscript and appendix. Replace bracket numbers with the final manuscript numbering during journal formatting.

1. COVID-19 Genomics UK (COG-UK) Consortium. An integrated national scale SARS-CoV-2 genomic surveillance network. *The Lancet Microbe*. 2020;1(3):e99--e100. doi:10.1016/S2666-5247(20)30054-9.

2. da Silva Filipe A, Shepherd JG, Williams T, et al. Genomic epidemiology reveals multiple introductions of SARS-CoV-2 from mainland Europe into Scotland. *Nature Microbiology*. 2021;6:112--122. doi:10.1038/s41564-020-00838-z.

3. Aksamentov I, Roemer C, Hodcroft EB, Neher RA. Nextclade: clade assignment, mutation calling and quality control for viral genomes. *Journal of Open Source Software*. 2021;6(67):3773. doi:10.21105/joss.03773.

4. Scottish Government. Scottish Index of Multiple Deprivation 2020: technical notes. 2020.

5. Tamura K, Nei M. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. *Molecular Biology and Evolution*. 1993;10(3):512--526. doi:10.1093/oxfordjournals.molbev.a040023.

6. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*. 2019;9:5233. doi:10.1038/s41598-019-41695-z.

7. Mullahy J. Specification and testing of some modified count data models. *Journal of Econometrics*. 1986;33(3):341--365. doi:10.1016/0304-4076(86)90002-3.

8. Cameron AC, Trivedi PK. *Regression Analysis of Count Data*. 2nd ed. Cambridge University Press; 2013.

9. Hilbe JM. *Negative Binomial Regression*. 2nd ed. Cambridge University Press; 2011.

10. White H. Maximum likelihood estimation of misspecified models. *Econometrica*. 1982;50(1):1--25. doi:10.2307/1912526.

11. Liang K-Y, Zeger SL. Longitudinal data analysis using generalized linear models. *Biometrika*. 1986;73(1):13--22. doi:10.1093/biomet/73.1.13.

12. Zeileis A. Econometric computing with HC and HAC covariance matrix estimators. *Journal of Statistical Software*. 2004;11(10):1--17. doi:10.18637/jss.v011.i10.

13. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*. 2020;17:261--272. doi:10.1038/s41592-019-0686-2.

14. Seabold S, Perktold J. statsmodels: Econometric and statistical modelling with Python. In: *Proceedings of the 9th Python in Science Conference*. 2010;92--96. doi:10.25080/Majora-92bf1922-011.

15. Csárdi G, Nepusz T. The igraph software package for complex network research. *InterJournal, Complex Systems*. 2006;1695.

16. Tange O. GNU Parallel 20240422. Zenodo. 2024. doi:10.5281/zenodo.11098743.

---

## 21. Notes for manuscript integration

- Keep Part A in the main manuscript or supplementary methods, depending on target journal word limits.
- Keep Part B as a reproducibility appendix.
- During final formatting, harmonise the reference numbers in Part A with the manuscript bibliography.
- The code excerpts in this appendix are explanatory and should be checked against the final committed repository version before submission.
