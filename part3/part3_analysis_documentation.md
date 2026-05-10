# Part 3 Analysis Documentation

## Policy phases, Alpha emergence, and SARS-CoV-2 genomic cluster structure in Scotland

This document records the current Part 3 analysis implemented in:

```text
part3/part3_analysis.py
```

and the supplementary meta-cluster extension implemented in:

```text
part3/notebooks/alpha_pre_l2_meta_cluster_network.ipynb
part3/notebooks/alpha_top6_meta_cluster_demographics_over_time.ipynb
part3/manuscript/make_alpha_meta_cluster_supp_figures.py
```

It is written as an implementation record: what data were used, how policy
periods were attached, what whole-epidemic descriptive summaries were computed,
what interrupted time-series models were fitted, how the Alpha emergence case
study was constructed, and what the pre-L2 Alpha meta-cluster extension found.

The Part 3 research questions are:

1. How did median genomic cluster size and geographic dispersion change across
   Scottish COVID-19 policy phases over the full epidemic period?
2. Did specific policy transitions — the pre-tier tightening (T1 onset), the
   lockdown-to-stay-local transition (L2 to SL), and the move to near-normal
   conditions (NN onset) — coincide with measurable changes in cluster
   structure and demographic mixing?
3. How did the Alpha variant (B.1.1.7) emerge and expand in Scotland during
   the F5 and L2 policy periods, and what does its growth advantage imply about
   the timing of non-pharmaceutical interventions?
4. What were the characteristics of the major pre-L2 Alpha meta-clusters, and
   what mutation signatures distinguished the dominant component (AM001)?

Part 3 is deliberately not an exhaustive policy catalogue. The analysis keeps
the full 16-period context table, then develops four policy moments where
policy, variant advantage, and cluster structure are especially informative.
All associations between policy periods and genomic outcomes are descriptive
and strongly confounded.

---

## 1. Source Data

### 1.1 Cluster-level data

Policy-period summaries and ITS analyses use the Part 1 cluster cache:

```text
part1/main/cache/main_cluster_table.parquet
```

This cache is treated as the primary cluster-level input because it already
contains the Part 1 cluster outcomes and excess-mixing metrics. The unit is
one inferred cluster/window row.

### 1.2 Sequence-level data

Alpha emergence and mutation analyses use:

```text
data/processed/scotland_clustering_analysis_dataset.parquet
```

with fallback to the raw Nextclade TSV from `config.yaml`:

```text
data/raw/cog_all_scotland_nextclade.tsv
```

The same QC filter and primary resolution used in Parts 1 and 2 are applied:

| Filter            | Value  |
|-------------------|-------:|
| Nextclade QC      | `good` |
| Leiden resolution |  `0.3` |

Alpha sequences are identified from Pango lineages beginning `B.1.1.7`.
Mutation trajectories are rebuilt from Nextclade amino-acid substitutions for
`S:N501Y`, `S:A222V`, and supporting Alpha markers.

### 1.3 Policy period definitions

Policy periods are attached with `utils.policy.attach_period_pandas` using
`wn_mid_date`. Period definitions reside in `utils/policy.py` and are not
duplicated in Part 3 scripts.

---

## 2. Why This Is Part 3

Part 1 modelled deprivation effects on cluster size, geographic spread, and
within-cluster demographic mixing using regression models with surveillance
adjustments. Part 2 added vaccination characterisation without causal
inference. Part 3 places both within a policy and variant-phase context for
three reasons:

1. Cluster size and geographic dispersion varied substantially across policy
   periods, and this variation provides a useful qualitative check on whether
   the genomic cluster structure tracked plausible epidemiological signals.
2. The Alpha emergence in late 2020 — during the Five-Tier Framework (F5) and
   the second lockdown (L2) — is the clearest example in the dataset where a
   specific policy context coincided with a well-characterised variant-growth
   event.
3. The pre-L2 Alpha meta-cluster extension strengthens the timing argument by
   showing that pre-L2 Alpha burden was concentrated in a small number of
   high-amplification components, especially AM001.

---

## 3. Policy Period Definitions

The full context table covers 16 policy periods. Part 3 develops four as
selected phases and retains the remainder for supplementary context.

| Period | Label                     | Start date | End date   | Intensity | Role           |
|--------|---------------------------|------------|------------|----------:|----------------|
| E0     | Emergence                 | 2020-03-01 | 2020-03-23 |        15 | supplement     |
| L1     | First lockdown            | 2020-03-24 | 2020-05-28 |       100 | supplement     |
| P1     | Route map phase 1         | 2020-05-29 | 2020-06-18 |        72 | supplement     |
| P2     | Route map phase 2         | 2020-06-19 | 2020-07-09 |        52 | supplement     |
| P3     | Route map phase 3         | 2020-07-10 | 2020-10-01 |        30 | **selected**   |
| T1     | Pre-tier tightening       | 2020-10-02 | 2020-11-01 |        55 | **selected**   |
| F5     | Five-tier framework       | 2020-11-02 | 2021-01-04 |        65 | **selected**   |
| L2     | Second lockdown           | 2021-01-05 | 2021-04-01 |        95 | **selected**   |
| SL     | Stay local — Level 3      | 2021-04-02 | 2021-04-25 |        65 | **selected**   |
| L3     | Level 3                   | 2021-04-26 | 2021-05-16 |        55 | supplement     |
| L21    | Level 2 / Level 1         | 2021-05-17 | 2021-07-18 |        38 | supplement     |
| L0     | Level 0                   | 2021-07-19 | 2021-08-08 |        20 | **selected**   |
| NN     | Near-normal               | 2021-08-09 | 2021-11-28 |        10 | **selected**   |
| OM     | Omicron wave              | 2021-11-29 | 2022-01-23 |        42 | context        |
| FE     | Final easing              | 2022-01-24 | 2022-04-17 |        15 | context        |
| PR     | Post-restriction          | 2022-04-18 | 2023-05-05 |         3 | context        |

Policy intensity is a pre-specified ordinal score assigned in `utils/policy.py`.
The early periods (E0 through P2) fall outside the analysis-window coverage and
therefore contribute no cluster rows.

---

## 4. Whole-Epidemic Policy Context

### 4.1 Period-level descriptive summaries

Current descriptive results from `period_descriptives.csv` for observed periods:

| Period | Label                | Clusters | Sequences | Singleton fraction | Median log size | Median log datazones | Mean SIMD excess | Mean age excess |
|--------|----------------------|---------:|----------:|--------------------|----------------:|---------------------:|-----------------:|----------------:|
| P3     | Route map phase 3    |   1,570  |    10,625 | 51.1%              |            0.00 |                 0.00 |           −0.111 |          −0.094 |
| T1     | Pre-tier tightening  |   1,379  |     6,266 | 57.8%              |            0.00 |                 0.00 |           −0.087 |          −0.053 |
| F5     | Five-tier framework  |   2,428  |     7,470 | 56.5%              |            0.00 |                 0.00 |           −0.120 |          −0.071 |
| L2     | Second lockdown      |   9,247  |    54,710 | 45.8%              |            0.69 |                 0.00 |           −0.194 |          −0.062 |
| SL     | Stay local           |   1,689  |    10,828 | 38.2%              |            0.69 |                 0.69 |           −0.256 |          −0.084 |
| L3     | Level 3              |   1,169  |     9,898 | 39.3%              |            0.69 |                 0.69 |           −0.237 |          −0.084 |
| L21    | Level 2/1            |  10,184  |    74,290 | 54.2%              |            0.00 |                 0.00 |           −0.157 |          −0.085 |
| L0     | Level 0              |   5,245  |    14,313 | 58.6%              |            0.00 |                 0.00 |           −0.160 |          −0.087 |
| NN     | Near-normal          |  48,407  | 183,190   | 53.4%              |            0.00 |                 0.00 |           −0.173 |          −0.079 |
| OM     | Omicron wave         |  21,394  | 112,158   | 58.7%              |            0.00 |                 0.00 |           −0.146 |          −0.063 |
| FE     | Final easing         |  38,699  | 178,704   | 54.7%              |            0.00 |                 0.00 |           −0.146 |          −0.061 |
| PR     | Post-restriction     |  51,701  | 126,790   | 62.9%              |            0.00 |                 0.00 |           −0.161 |          −0.095 |

The L2 period contains the highest non-singleton fraction (45.8% singletons)
and the largest median log cluster size (0.69) outside the SL/L3 transition.
SL has the lowest singleton fraction across all observed periods (38.2%) and
the most negative mean SIMD excess discordance (−0.256), consistent with a
period when Alpha, during a restriction phase, produced the most homogeneous
clusters relative to expectation.

### 4.2 Clustering Rate and Dispersion Parameter k

Because cluster size is highly right-skewed across all periods, the log median
is an insensitive summary of the full distribution. Two additional per-period
metrics are reported in `period_clustering_dispersion.csv`:

**Clustering rate** quantifies what fraction of sequenced isolates belong to
transmission chains rather than appearing as isolated detections:

```
clustering_rate = (n_c - c) / n_c
```

where `n_c` is the total number of sequenced isolates (sum of all cluster
sizes in the period) and `c` is the total number of clusters. When every
cluster is a singleton, `c = n_c` and the rate is 0; when all isolates belong
to a single large cluster, the rate approaches 1. This is equivalent to
1 − (clusters / sequences).

**Dispersion parameter k** is estimated by the method of moments (MME):

```
k_hat = x̄² / (s² − x̄)
```

where `x̄` is the mean cluster size and `s²` is the sample variance, computed
across all clusters (including singletons). Overdispersion requires `s² > x̄`,
which holds for all observed periods. A smaller `k` indicates greater dominance
of a few large clusters ("superspreading-like" heterogeneity); as `k → ∞` the
distribution approaches Poisson. A complementary estimate `k_non_singleton` is
computed from the non-singleton population only.

The same two metrics are computed in parallel for **geographic spread**
(`cluster_n_datazones`), giving `geo_clustering_rate` and `geo_k_all`:

```
geo_clustering_rate = (sum(cluster_n_datazones) − c) / sum(cluster_n_datazones)
```

This is the fraction of "datazone-slots" that are secondary (i.e., belong to a
cluster spanning more than one datazone). `geo_k_all` is fitted to the
distribution of `cluster_n_datazones` with the same MME formula.

All four metrics (clustering_rate, k_all, geo_clustering_rate, geo_k_all) are
also computed per analysis window and stored in `weekly_summaries.csv`, where
they drive Figures 1 and 2. ITS segmented-regression fitted values for each
metric are included in the three `its_weekly_*.csv` files.

Computed from `part3/cache/policy_cluster_table.parquet` and saved to
`part3/tables/period_clustering_dispersion.csv`.

| Period | Label                | Clustering rate | k (all) | k (non-singleton) | Mean size | Var size |
|--------|----------------------|----------------:|--------:|------------------:|----------:|---------:|
| P3     | Route map phase 3    |          0.8522 |  0.1083 |            0.2060 |      6.77 |    429.7 |
| T1     | Pre-tier tightening  |          0.7799 |  0.1715 |            0.3588 |      4.54 |    124.9 |
| F5     | Five-tier framework  |          0.6750 |  0.2646 |            0.4721 |      3.08 |     38.9 |
| L2     | Second lockdown      |          0.8310 |  0.1444 |            0.2472 |      5.92 |    248.3 |
| SL     | Stay local           |          0.8440 |  0.3623 |            0.6148 |      6.41 |    119.9 |
| L3     | Level 3              |          0.8819 |  0.0763 |            0.1188 |      8.47 |    948.4 |
| L21    | Level 2/1            |          0.8629 |  0.0171 |            0.0325 |      7.29 |  3,118.9 |
| L0     | Level 0              |          0.6335 |  0.2157 |            0.3595 |      2.73 |     37.3 |
| NN     | Near-normal          |          0.7358 |  0.1429 |            0.2469 |      3.78 |    104.0 |
| OM     | Omicron wave         |          0.8093 |  0.0102 |            0.0197 |      5.24 |  2,692.7 |
| FE     | Final easing         |          0.7834 |  0.0240 |            0.0419 |      4.62 |    893.4 |
| PR     | Post-restriction     |          0.5922 |  0.0817 |            0.1269 |      2.45 |     76.1 |

Key patterns:

- k is below 1.0 in all periods, confirming strong overdispersion throughout
  the epidemic. The distribution is never Poisson-like.
- The lowest k values occur in long multi-variant periods (L21 k = 0.017,
  OM k = 0.010), where extreme variance is driven by mixing of many small
  post-variant clusters with a few large outbreak clusters.
- Among selected phases, SL has the highest k (0.362) — relatively less
  extreme overdispersion — consistent with Alpha declining and producing a
  less concentrated size distribution. F5 has the next highest (0.265).
- Clustering rate is highest in the L3/SL/L2 periods (Alpha wave, 0.83–0.88),
  meaning 83–88% of sequences in those periods were secondary cluster members
  rather than singletons. It is lowest in PR (0.592) and L0 (0.634).
- The variance column explains why log median is insensitive: variance in L21
  (3,119) is 250× that in F5 (39), but the log median is 0.0 in both.

### 4.3 Policy intensity correlations

Spearman correlations between policy intensity and weekly outcomes across all
134 observed analysis windows (from `intensity_correlations.csv`). All
correlations are described as descriptive and confounded.

| Weekly outcome              | Spearman ρ | p-value | Note                          |
|-----------------------------|-----------:|--------:|-------------------------------|
| Total clusters              |     −0.256 |  0.0028 | Weak negative                 |
| Non-singleton clusters      |     −0.098 |  0.2622 | Not clearly different from 0  |
| Singleton fraction          |     −0.621 | < 0.001 | Strong negative               |
| Median log cluster size     |     +0.741 | < 0.001 | Strong positive               |
| Median log datazones        |     +0.576 | < 0.001 | Moderate positive             |
| Mean SIMD excess discordance|     +0.019 |  0.832  | Near null                     |
| Mean age excess discordance |     +0.587 | < 0.001 | Moderate positive             |
| Mean sex excess discordance |     −0.205 |  0.0176 | Weak negative                 |
| Mean profile excess disc.   |     +0.260 |  0.0024 | Weak positive                 |

Higher policy intensity is strongly correlated with larger median cluster size
(ρ = 0.741) and lower singleton fraction (ρ = −0.621). These are descriptive
patterns confounded by variant phase, population immunity, testing eligibility,
sequencing capacity, and calendar time.

---

## 5. Interrupted Time-Series Analysis

### 5.1 Model form

The ITS model is a linear segmented regression:

```text
y_t = beta0 + beta1*t + beta2*post + beta3*(post*t) + error
```

where `t` is week index centred on the transition date and `post` is an
indicator for observations after the transition. `beta2` is the immediate
level change and `beta3` is the change in slope after the transition.

The primary window is ±8 weeks around each transition. Sensitivity fits are
run at ±6, ±10, and ±12 weeks.

### 5.2 Outcomes

The main ITS outcomes are median log cluster size and median log datazones.
Mixing outcomes (mean SIMD excess discordance, mean age excess discordance)
are fitted in the supplementary figure.

### 5.3 Selected transitions

| Transition  | Label       | From → To    | Date       |
|-------------|-------------|--------------|------------|
| `t1_onset`  | T1 onset    | P3 → T1      | 2020-10-02 |
| `l2_to_sl`  | L2 to SL    | L2 → SL      | 2021-04-02 |
| `nn_onset`  | NN onset    | L0 → NN      | 2021-08-09 |

### 5.4 Primary ITS results (±8-week window)

Primary coefficients from `its_coefficients.csv`:

**T1 onset (P3 → T1, 2020-10-02)**

| Outcome                   | Level change (post) | Slope change (post × t) | adj-R² |
|---------------------------|--------------------:|------------------------:|-------:|
| Median log cluster size   |  −0.080 (p = 0.197) |  −0.086 (p < 0.001)     |  0.505 |
| Median log datazones      |  −0.027 (p = 0.866) |  −0.077 (p = 0.028)     |  0.314 |
| Mean SIMD excess disc.    |  +0.032 (p = 0.185) |  −0.027 (p < 0.001)     |  0.635 |
| Mean age excess disc.     |  +0.010 (p = 0.655) |  −0.009 (p = 0.075)     |  0.575 |

The T1 onset coincided with a significant decline in the growth trend of median
log cluster size and log datazones, but no significant immediate level change.
The SIMD mixing slope also declined significantly after T1 onset.

**L2 to SL (L2 → SL, 2021-04-02)**

| Outcome                   | Level change (post) | Slope change (post × t) | adj-R² |
|---------------------------|--------------------:|------------------------:|-------:|
| Median log cluster size   |  −0.171 (p = 0.102) |  −0.072 (p = 0.001)     |  0.677 |
| Median log datazones      |  −0.355 (p < 0.001) |  −0.014 (p = 0.344)     |  0.437 |
| Mean SIMD excess disc.    |  −0.007 (p = 0.558) |  +0.022 (p < 0.001)     |  0.754 |
| Mean age excess disc.     |  +0.007 (p = 0.264) |  +0.007 (p = 0.002)     |  0.680 |

The L2-to-SL transition coincided with an immediate and significant drop in
median log datazones (−0.355 log units, p < 0.001) at the level-change term,
consistent with geographic clustering contracting as Alpha wound down during
SL. The SIMD mixing slope turned positive after the transition, suggesting
clusters became progressively more mixed by deprivation as the restriction
phase lifted.

**NN onset (L0 → NN, 2021-08-09)**

| Outcome                   | Level change (post) | Slope change (post × t) | adj-R² |
|---------------------------|--------------------:|------------------------:|-------:|
| Median log cluster size   |  +0.116 (p = 0.016) |  +0.051 (p < 0.001)     |  0.725 |
| Median log datazones      |  +0.312 (p = 0.036) |  +0.048 (p = 0.038)     |  0.354 |
| Mean SIMD excess disc.    |  −0.007 (p = 0.749) |  +0.001 (p = 0.882)     | −0.151 |
| Mean age excess disc.     |  −0.006 (p = 0.701) |  −0.003 (p = 0.448)     | −0.107 |

The NN onset coincided with significant immediate increases in both median
log cluster size (+0.116 log units, p = 0.016) and median log datazones
(+0.312 log units, p = 0.036), and the growth slope for both outcomes also
increased significantly after the transition. Mixing outcomes showed no
significant change, consistent with the near-null SIMD-intensity correlation
in the whole-epidemic analysis.

### 5.5 Sensitivity window results

Coefficients are written to `its_coefficients.csv` for all four windows (±6,
±8, ±10, ±12 weeks). The direction of the primary slope-change terms is
consistent across the ±8, ±10, and ±12-week windows. The ±6-week window
has fewer observations and wider confidence intervals. Key results are
robust to the window choice for the cluster-size and datazones outcomes at
T1 onset and L2-to-SL. The NN-onset mixing outcomes remain non-significant
at all window widths.

---

## 6. Alpha Emergence Case Study: F5 and L2

### 6.1 Phase definitions

Alpha emergence is rebuilt from the primary sequence-level table. The three
reported phases are defined by observed window mid-dates:

| Phase                  | Windows  | Window IDs | Observed mid-date range         |
|------------------------|----------|------------|---------------------------------|
| `cryptic_early`        | W016–W021| W016–W021  | 2020-10-27 to 2020-12-01        |
| `multi_region_expansion`| W022–W024| W022–W024 | 2020-12-08 to 2020-12-22        |
| `f5_l2_bridge`         | W025     | W025       | 2020-12-29                      |

### 6.2 Phase descriptives

Current results from `alpha_phase_summary.csv`:

| Phase                    | Total sequences | Alpha sequences | Alpha fraction | Alpha clusters | Median cluster size | Max cluster size | Health boards |
|--------------------------|----------------:|----------------:|---------------:|---------------:|--------------------:|-----------------:|--------------:|
| Cryptic/early            |           2,484 |              51 |          2.1%  |             15 |                 5.0 |               33 |             7 |
| Multi-region expansion   |           1,387 |             291 |         21.0%  |            113 |                 1.0 |               99 |            11 |
| F5/L2 bridge (W025)      |           1,065 |             458 |         43.0%  |            102 |                 1.0 |              135 |            11 |

During the cryptic/early phase, Alpha was predominantly detected in Greater
Glasgow and Clyde (77 of 51 sequences — note: sequences may be counted across
windows). By the multi-region expansion phase, Alpha had spread to 11 health
boards and 28 local authorities. In the F5/L2 bridge window, Alpha comprised
43% of all sequences and was present across 29 local authorities, indicating
near-national spread before L2 was enacted.

The largest early Alpha cluster (max size 99 during expansion, 135 at the
bridge) reflects the concentrated AM001 component documented in the
meta-cluster extension (Section 7).

---

## 7. Pre-L2 Alpha Meta-Cluster Extension

### 7.1 Network construction

The supplementary Alpha meta-cluster analysis is implemented as exploratory
notebooks and converted to manuscript supplementary figures. The meta-cluster
network is restricted to unique Alpha sequences collected before the L2 start
date (2021-01-05). Network nodes are Alpha-containing cluster/window
assignments. Edges connect nodes in adjacent rolling windows when they share
at least one sequence. Connected components of this graph are treated as
Alpha meta-clusters.

This construction captures continuity across overlapping rolling-window
cluster calls. It should not be interpreted as proving that a component is a
single transmission chain or that a specific superspreading event occurred.

### 7.2 Network summary

| Quantity                               | Value |
|----------------------------------------|------:|
| Unique pre-L2 Alpha sequences          |   442 |
| Connected components (meta-clusters)   |    78 |
| Single-sequence meta-clusters          |    49 |
| Meta-clusters with ≥10 unique sequences|     6 |
| Ambiguous sequence assignments         |     0 |

### 7.3 Six largest meta-clusters

Results from `alpha_pre_l2_meta_cluster_summary.csv` and
`alpha_top6_meta_cluster_compact_summary.csv`:

| Meta-cluster | Unique sequences | First collection | Last collection  | Nodes | Edges | Main health board signal       |
|--------------|----------------:|------------------|------------------|------:|------:|-------------------------------|
| AM001        |             234 | 2020-11-04       | 2021-01-04       |    46 |    47 | Greater Glasgow and Clyde (137; 58.5%) |
| AM003        |              19 | 2020-12-07       | 2020-12-18       |     5 |     4 | Grampian (14; 73.7%)           |
| AM004        |              19 | 2020-12-24       | 2021-01-04       |     5 |     4 | GGC (6), Tayside (4), Grampian (4) |
| AM035        |              17 | 2020-12-27       | 2021-01-04       |     2 |     1 | Highland (14; 82.4%)           |
| AM024        |              13 | 2020-12-23       | 2021-01-04       |     3 |     2 | Borders (13; 100.0%)           |
| AM034        |              10 | 2020-12-27       | 2021-01-04       |     3 |     2 | Greater Glasgow and Clyde (7; 70.0%) |

Together these six meta-clusters contain 312 of 442 pre-L2 Alpha sequences
(70.6%). AM001 alone contains 234 of 442 sequences (52.9%).

### 7.4 AM001 demographics

From `alpha_top6_meta_cluster_compact_summary.csv`:

| Variable          | AM001 composition                                   |
|-------------------|-----------------------------------------------------|
| Age group         | 40–64 (85; 36.3%), 65+ (79; 33.8%), 18–39 (60; 25.6%) |
| Sex               | Female (144; 61.5%), Male (90; 38.5%)               |
| SIMD quintile     | Q1 most deprived (86; 36.8%), Q4 (42; 17.9%), Q2 (38; 16.2%) |
| Test reason       | Missing (122; 52.1%), Symptomatic (101; 43.2%)      |

The concentration of AM001 sequences in Q1 (most deprived, 36.8%) is notable,
though this reflects population distribution in Greater Glasgow and Clyde rather
than a targeted ascertainment of deprived areas.

### 7.5 Mutation-signature summaries

Raw Nextclade amino-acid substitutions are joined by `sequence_id`. For each
sufficiently large meta-cluster, recurrent substitutions are tested for
enrichment against the other pre-L2 Alpha sequences using Fisher exact tests
with Benjamini-Hochberg false-discovery correction.

The main AM001-associated marker is `ORF1a:L730F`, present in 199/234 AM001
sequences (85.0%) and 38/208 non-AM001 pre-L2 Alpha sequences (18.3%). It is
strongly enriched for AM001 but not private to AM001. Other candidate signatures
include `ORF3a:L15F`, `ORF1a:F2780L`, `ORF1b:K1383R`, `S:V1129I`, and
`ORF8:K68*` for smaller components.

Results are written to:

```text
part3/notebooks/tables/alpha_pre_l2_meta_cluster_signature_mutations.csv
part3/notebooks/tables/alpha_pre_l2_meta_cluster_signature_impact_summary.csv
```

---

## 8. Growth Models and Counterfactuals

### 8.1 Growth model specification

Binomial GLMs are fitted to weekly marker frequencies:

- Alpha marker `S:N501Y` during the F5 period
- Alpha marker `S:N501Y` during the L2 period
- B.1.177 comparator marker `S:A222V` during the L2 period

The primary model uses weekly positive tests as frequency weights. Sensitivity
models use sequence-count weights, equal weekly weights, and coverage-adjusted
weights.

### 8.2 Primary growth model results

From `alpha_growth_params.csv` (positive-test-weighted):

| Analysis               | Marker  | Period | Windows | Slope/day | Slope/week | OR/week      | Doubling/halving time | Pseudo-R² |
|------------------------|---------|--------|--------:|----------:|-----------:|:-------------|----------------------:|----------:|
| Alpha F5 S:N501Y       | S:N501Y | F5     |       9 | +0.0851   |  +0.596    | 1.815 (1.802–1.828) | 8.1 days doubling | 0.936 |
| Alpha L2 S:N501Y       | S:N501Y | L2     |      13 | +0.0618   |  +0.433    | 1.542 (1.534–1.549) | 11.2 days doubling | 0.974 |
| B.1.177 L2 S:A222V     | S:A222V | L2     |      13 | −0.0655   |  −0.459    | 0.632 (0.629–0.636) | 10.6 days halving | 0.971 |

Alpha grew faster during F5 (doubling time 8.1 days) than during L2 (11.2
days), consistent with L2 slowing but not halting Alpha expansion. The B.1.177
comparator was declining at a rate symmetric with Alpha's L2 growth, with a
halving time of 10.6 days. The primary model fits are robust across weighting
schemes (sensitivity results in `alpha_growth_model_sensitivity.csv`).

### 8.3 Counterfactual timing analysis

The counterfactual analysis projects `S:N501Y` frequency under earlier switches
to the fitted L2 growth rate, using the L2 binomial growth model with the
initial frequency at the nearest observed expansion window.

Results from `alpha_counterfactual_projections.csv`:

| Scenario               | Switch date | Nearest window | Starting N501Y frequency | Projected 50% date | Days vs actual |
|------------------------|-------------|----------------|-------------------------:|--------------------|---------------:|
| Actual L2 start        | 2021-01-05  | W026           |                    60.4% | 2020-12-29         |              0 |
| Earlier L2: 2020-12-08 | 2020-12-08  | W022           |                    17.7% | 2021-01-09         |            +11 |
| Earlier L2: 2020-12-02 | 2020-12-02  | W021           |                     3.2% | 2021-01-11         |            +13 |
| L2 from F5 start       | 2020-11-02  | W017           |                     0.7% | 2021-01-23         |            +25 |

Under the most aggressive scenario (L2 from the F5 start on 2020-11-02), the
projected 50% S:N501Y frequency date would be 2021-01-23, only 25 days later
than the actual trajectory. An earlier switch from 2020-12-08 delays projected
50% by only 11 days. These are descriptive fitted scenarios, not causal policy
estimates. They depend on the assumption that the L2 growth rate would have
applied unchanged from an earlier start date, which requires extrapolating the
growth model well beyond its observation window.

---

## 9. Outputs

### 9.1 Primary analysis tables (part3/tables/)

| File                                | Contents                                                                  |
|-------------------------------------|---------------------------------------------------------------------------|
| `period_descriptives.csv`           | Per-period cluster counts, singleton fractions, medians, mixing means     |
| `period_clustering_dispersion.csv`  | Per-period clustering rate and MME dispersion parameter k (all and non-singleton) |
| `weekly_summaries.csv`              | Weekly cluster outcomes, dominant policy period, and policy intensity     |
| `intensity_correlations.csv`        | Spearman correlations between policy intensity and weekly outcomes         |
| `its_coefficients.csv`             | ITS model coefficients for all three transitions and four window widths   |
| `its_weekly_t1_onset.csv`           | Primary fitted weekly data for the T1 onset transition                    |
| `its_weekly_l2_to_sl.csv`           | Primary fitted weekly data for the L2-to-SL transition                   |
| `its_weekly_nn_onset.csv`           | Primary fitted weekly data for the NN onset transition                    |
| `alpha_phase_summary.csv`           | Per-phase Alpha sequence counts, cluster counts, and geographic spread    |
| `alpha_cluster_emergence.csv`       | Alpha cluster-level emergence data                                        |
| `alpha_health_board_weekly.csv`     | Weekly Alpha sequence counts by health board                              |
| `alpha_local_authority_weekly.csv`  | Weekly Alpha sequence counts by local authority                           |
| `alpha_mutation_trajectories.csv`   | Weekly S:N501Y and S:A222V frequency trajectories                        |
| `alpha_growth_params.csv`           | Primary binomial GLM growth estimates for all three marker-period pairs   |
| `alpha_growth_model_sensitivity.csv`| Growth estimates under four weighting schemes                             |
| `alpha_counterfactual_projections.csv`| Counterfactual 50% dates under earlier L2 scenarios                   |
| `alpha_counterfactual_trajectories.csv`| Weekly projected S:N501Y trajectories under each counterfactual       |

### 9.2 Meta-cluster notebook tables (part3/notebooks/tables/)

| File                                             | Contents                                                      |
|--------------------------------------------------|---------------------------------------------------------------|
| `alpha_pre_l2_meta_cluster_summary.csv`          | Per-meta-cluster node count, sequence count, geography        |
| `alpha_pre_l2_meta_cluster_sequence_membership.csv` | Sequence-to-meta-cluster assignments                       |
| `alpha_pre_l2_meta_cluster_signature_mutations.csv` | Fisher-exact enrichment results for candidate signatures   |
| `alpha_pre_l2_meta_cluster_signature_impact_summary.csv` | Signature enrichment summary                          |
| `alpha_pre_l2_meta_cluster_nodes.csv`            | Network node attributes                                       |
| `alpha_pre_l2_meta_cluster_edges.csv`            | Network edge list                                             |
| `alpha_pre_l2_signature_mutation_trajectories.csv` | Calendar-time trajectories for signature mutations          |
| `alpha_top6_meta_cluster_sequence_metadata.csv`  | Sequence-level metadata for the six largest meta-clusters     |
| `alpha_top6_meta_cluster_weekly_counts.csv`      | Weekly sequence counts by meta-cluster                        |
| `alpha_top6_meta_cluster_overall_composition.csv`| Overall demographic and geographic composition per meta-cluster |
| `alpha_top6_meta_cluster_compact_summary.csv`    | Compact demographics/geography/testing summary per meta-cluster |
| `alpha_top6_meta_cluster_composition_chi_square.csv` | Chi-square tests for compositional differences            |
| `alpha_top6_meta_cluster_testing_composition.csv`| Test-reason composition per meta-cluster                      |
| `alpha_top6_meta_cluster_top_test_reasons.csv`   | Top test reasons per meta-cluster                             |
| `alpha_top6_meta_cluster_weekly_composition.csv` | Weekly demographic composition per meta-cluster               |

### 9.3 Manuscript figures (part3/manuscript/figures/)

Six figures generated by `part3/manuscript/make_figures.py` and
`part3/manuscript/make_alpha_meta_cluster_supp_figures.py`. Figures are
written as PDF, PNG, and LZW-compressed TIFF at 600 dpi.

| Figure stem                                   | Description                                                                          |
|-----------------------------------------------|--------------------------------------------------------------------------------------|
| `fig1_policy_timeline_cluster_structure`      | Three-panel weekly overview: (A) median cluster size + IQR, (B) clustering rate, (C) dispersion k̂ |
| `fig2_selected_policy_transitions`            | ±8-week ITS (main): log-median cluster size + log-median datazones, with IQR error bars |
| `fig3_alpha_emergence_f5_l2`                 | S:N501Y and S:A222V frequencies plus health-board expansion during F5/L2             |
| `fig4_alpha_counterfactual_timing`            | Fitted and counterfactual S:N501Y trajectories                                       |
| `supp_fig1_its_mixing_outcomes`               | ±8-week ITS plots for SIMD and age excess-discordance outcomes                       |
| `supp_fig2_alpha_meta_cluster_amplification`  | AM001 amplification and meta-cluster growth curves                                   |
| `supp_fig2a_its_clustering_rate`              | ±8-week ITS (supplementary A): clustering rate — cluster size (left) and geographic spread (right) |
| `supp_fig2b_its_dispersion`                   | ±8-week ITS (supplementary B): dispersion k̂ — cluster size (left) and geographic spread (right) |
| `supp_fig3_alpha_top6_meta_cluster_context`   | Demographics and geography of the six largest pre-L2 Alpha meta-clusters             |

---

## 10. Reproducibility

Run the full Part 3 pipeline from the repository root:

```bash
# Step 1 — run the primary Part 3 analysis (requires Part 1 cache)
conda run -n PhD python part3/part3_analysis.py

# Step 2 — generate main manuscript figures
conda run -n PhD python part3/manuscript/make_figures.py

# Step 3 — generate supplementary Alpha meta-cluster figures
#           (requires notebooks to have been executed first)
conda run -n PhD python part3/manuscript/make_alpha_meta_cluster_supp_figures.py
```

The meta-cluster notebooks must be executed separately:

```bash
conda run -n PhD jupyter nbconvert --to notebook --execute \
  part3/notebooks/alpha_pre_l2_meta_cluster_network.ipynb

conda run -n PhD jupyter nbconvert --to notebook --execute \
  part3/notebooks/alpha_top6_meta_cluster_demographics_over_time.ipynb
```

Important defaults in `part3_analysis.py`:

| Setting                    | Default                                       |
|----------------------------|-----------------------------------------------|
| Primary QC filter          | `good`                                        |
| Primary Leiden resolution  | `0.3`                                         |
| Alpha lineage prefix       | `B.1.1.7`                                     |
| ITS primary window         | ±8 weeks                                      |
| ITS sensitivity windows    | ±6, ±10, ±12 weeks                            |
| Growth model weight scheme | Positive-test-weighted (primary)              |
| Counterfactual scenarios   | 4 (actual, Dec-08, Dec-02, F5 start)          |

---

## 11. Current Main Interpretation

The current Part 3 results support four linked conclusions:

1. **Policy intensity correlates with cluster structure but is strongly
   confounded.** Higher policy intensity is associated with larger median
   cluster size (ρ = 0.741) and lower singleton fraction (ρ = −0.621) across
   134 analysis windows. These correlations run opposite to the naive
   expectation that restrictions reduce cluster size, because high-intensity
   restrictions coincided with the Alpha and Delta waves, which drove larger
   and more geographically dispersed clusters regardless of policy. Policy
   intensity shows near-null correlation with SIMD excess discordance (ρ = 0.019).
   Because cluster size is strongly right-skewed, per-period clustering rate and
   the MME dispersion parameter k complement the log-median summary. Clustering
   rate peaked at 0.88 in the Alpha-wave periods (L3, SL), meaning 88% of
   sequences were secondary cluster members rather than isolated detections.
   k remained below 1.0 in all periods (strong overdispersion throughout), with
   the lowest values in long multi-variant periods (OM k = 0.010; L21 k = 0.017)
   reflecting mixing of many small clusters with a few very large outbreak events,
   and relatively higher values in F5 (k = 0.265) and SL (k = 0.362).

2. **Specific policy transitions coincided with changes in cluster structure
   trends, but immediate level changes were inconsistent.** At T1 onset and
   L2-to-SL, the growth slope of median log cluster size declined significantly.
   The L2-to-SL transition showed a significant immediate drop in median log
   datazones (−0.355 log units) but no significant immediate cluster-size
   level change. The NN onset showed significant immediate increases in both
   outcomes, consistent with the Delta wave coinciding with eased restrictions.
   Mixing outcomes did not change significantly at the NN transition.

3. **Alpha expanded rapidly during F5 before L2 was enacted.** The growth
   advantage of Alpha during F5 implies a doubling time of 8.1 days. By the
   F5/L2 bridge window (W025), Alpha comprised 43% of sequenced cases and
   had spread to 29 local authorities. The counterfactual timing analysis
   suggests that even a switch to L2-equivalent conditions from 2020-12-08
   would have delayed the Alpha 50% frequency date by only 11 days, and
   switching from the F5 start (2020-11-02) by only 25 days. These projections
   are highly sensitive to model assumptions and should not be treated as
   causal policy estimates.

4. **Pre-L2 Alpha burden was concentrated in a small number of high-amplification
   meta-clusters, especially AM001.** AM001 accounted for 234 of 442 pre-L2
   Alpha sequences (52.9%), was dominated by Greater Glasgow and Clyde (58.5%),
   and carried a signature substitution `ORF1a:L730F` in 85% of its sequences.
   The five other meta-clusters with ≥10 sequences account for an additional
   17.6% of pre-L2 Alpha sequences, suggesting early Alpha was not yet
   randomly distributed across Scotland but was concentrated in a small number
   of expanding lineages.

---

## 12. Caveats

1. Policy periods are strongly entangled with variant replacement, changing
   population immunity, testing eligibility and behaviour, sequencing coverage,
   hospital burden, and seasonality. No Part 3 result should be interpreted
   as causal evidence that a policy caused a cluster outcome.
2. The analysis is associational and descriptive throughout. The ITS models
   estimate correlation between timing of a policy transition and a change in
   cluster-structure trends, not a policy effect.
3. Sequence data do not identify a specific superspreading event. Components
   with rapid growth and large sequence counts are described as
   high-amplification or superspreading-like meta-clusters only in the
   descriptive sense.
4. The counterfactual timing analysis assumes that the L2 growth rate would
   have applied unchanged from an earlier start date. This requires extrapolating
   the binomial growth model well outside its observation window and ignores
   the possibility that restrictions themselves alter the growth trajectory.
5. The meta-cluster extension uses rolling-window overlap as a proxy for
   within-chain linkage. Components with rapid growth and many shared sequences
   across adjacent windows will appear as large meta-clusters regardless of
   whether the underlying transmission is a single chain or multiple independent
   introductions.
6. Early policy periods (E0 through P2) fall outside the analysis-window
   coverage. Results for P3 and T1 are based on a relatively small number of
   clusters (1,570 and 1,379 respectively) compared with later periods.
7. An earlier restriction scenario would not necessarily have prevented Alpha
   establishment; the counterfactual trajectories show modest delays in Alpha
   reaching 50% frequency, not prevention of establishment.
8. The AM001 demographic composition reflects the population structure of
   Greater Glasgow and Clyde, not a selective transmission pattern. The
   concentration in SIMD Q1 (36.8%) is broadly consistent with GGC's deprivation
   distribution.
