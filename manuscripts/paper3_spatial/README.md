# Paper 3 — Spatial structure of SARS-CoV-2 transmission clusters in Scotland

## One-line scope

How the *geographic footprint* of SARS-CoV-2 transmission clusters — measured from data-zone population-weighted centroids — varies with deprivation, urban-rural status, and VOC epoch, and whether spatial signatures of onward transmission attenuated as the pandemic progressed.

## Primary hypotheses

1. Within-cluster geographic spread (mean pairwise distance between members) is smaller in clusters drawn from more deprived data zones, consistent with localised, high-density transmission.
2. The **distance–decay** of co-membership (probability of being in the same cluster vs. pairwise distance) is steepest during Alpha/Delta and flattens during Omicron.
3. Urban clusters are smaller but more spatially concentrated than rural clusters; rural clusters are smaller-N but spatially dispersed.

## Unit of analysis

Primary: one row per `(window_id, cluster_id)` summarising the cluster's spatial footprint. Distance-decay analyses use pairs of sequences within a window.

## Figure list

| # | File                                         | Description |
|---|----------------------------------------------|-------------|
| 1 | `figures/fig1_cluster_centroid_map.py`       | Map of Scotland, one point per cluster centroid, size ∝ n_sequences, colour = dominant VOC |
| 2 | `figures/fig2_within_cluster_spread.py`      | Distribution of within-cluster mean pairwise distance (km), by SIMD quintile and by VOC epoch |
| 3 | `figures/fig3_distance_decay.py`             | Probability that two sequences share a cluster vs. pairwise distance (log-log), one line per VOC epoch |
| 4 | `figures/fig4_urban_rural_footprint.py`      | Boxplot comparing urban vs. rural cluster sizes and footprints (rural = SIMD access domain decile ≥ 8 proxy) |
| 5 | `figures/fig5_spread_by_deprivation.py`      | Hexbin of (log cluster size, log footprint) with deprivation-gradient colour overlay |

## Statistical conventions (read once, then skim the figure notes)

- **Distance** is Euclidean between **British National Grid (EPSG:27700) data-zone population-weighted centroids**, reported in **kilometres** (metres ÷ 1000). This is a DZ-level proxy, so *all* within-DZ pairs collapse to 0 km and the smallest resolvable distance is the nearest-neighbour DZ centroid distance (tens to hundreds of metres in urban cores; several km in rural areas).
- **Mean pairwise km** per cluster — the mean of `{d(i,j) : i < j, i,j ∈ cluster}`. Singletons and co-DZ-only clusters have mean = 0 by construction. A cluster with mean ≈ 10 km is confined to one city; 50+ km crosses council boundaries; 100+ km is transnationally dispersed within Scotland.
- **Bounding-box diagonal** — the Euclidean diagonal of the axis-aligned BNG bounding box of the cluster's DZ centroids. Sensitive to outliers, good for answering "how spread out is this cluster at its extremes?".
- **Co-cluster probability (distance-decay)** — of all sampled within-window sequence pairs at pairwise distance `d`, the fraction that fall in the same cluster. Under panmictic mixing this would be ~constant in `d`; steep negative slope on a log-log plot = strong locality.
- **Wilson 95% CI** — used instead of normal-approximation on binomial proportions because some distance bins have small `k/n`. Formula-driven; no bootstrap needed.
- **Linear regression for log-spread** (`models.spatial_mixing`) — on `log(mean_pairwise_km + 1) ~ simd_quintile_mode + voc`. The coefficient on (e.g.) `simd_quintile_mode == 1` is interpretable as "expected log-km difference vs. the Q5 reference" — exponentiate and subtract 1 for a percentage-point interpretation of geographic spread.
- **Log-log axes in Fig. 3** make distance-decay power laws visible as straight lines; the slope is the decay exponent.
- **Urban / rural proxy** (Fig. 4) is computed from the SIMD **access domain rank** (low rank = access-deprived = rural): bottom 20% ⇒ Rural, top 20% ⇒ Urban, middle 60% ⇒ Mixed. This is a convenience proxy that correlates strongly but not perfectly with the Scottish Government 6-fold urban/rural classification — see caveat 2.

## How to read each figure

### Fig. 1 — Map of cluster centroids by VOC

Small multiples: one map panel per dominant VOC (Alpha / Delta / Omicron). Every point is one cluster placed at the mean BNG easting/northing of its members' data zones; **point size is `6 + 6·log1p(n_sequences)`** (so singletons are tiny and a cluster of 1000 is 4–5× bigger than a cluster of 10); **colour encodes the VOC** via `style.WHO_VOC_PALETTE`. The faint Scotland coastline is drawn if the DZ shapefile is present (see Running) — otherwise axes are BNG-tick-only. Read it for (a) **spatial coverage** of sampling (dense in Glasgow / Edinburgh / Aberdeen; sparse in the Highlands), (b) **VOC penetration** (any white regions in Alpha that fill in during Delta?), and (c) the per-panel cluster count in the subtitle — a sanity check on sample sizes.

### Fig. 2 — Within-cluster geographic spread by SIMD quintile × VOC epoch

Twin panels. **Panel A (ridgeline)** — one ridge per SIMD quintile of `log(1 + mean_pairwise_km)`, normalised to the same peak height; x-axis is log1p-km. A ridge concentrated at 0 = clusters confined to one or two DZs (localised); a right-shifted ridge = geographically dispersed clusters. **Hypothesis 1 (paper 3) predicts Q1 (most deprived) sits further left than Q5**. **Panel B (boxplot matrix)** — the same metric but stratified by (epoch × SIMD quintile). Columns are epochs, and within each column five narrow boxplots are the five quintiles. Compare within-epoch patterns (deprivation gradient holding) and across-epoch trends (localisation weakening under Omicron). Singletons are excluded because their `mean_pairwise_km = 0` by definition and would dominate the left edge.

### Fig. 3 — Distance-decay of cluster co-membership

One line per VOC epoch on a **log-log plot**. X-axis: pairwise DZ-centroid distance (km, log). Y-axis: probability that two sequences sampled from the same sliding window fall in the same inferred cluster (log). Each point is a log-spaced distance bin; the band is a **Wilson 95% CI** on the binomial proportion. **Read it as a spatial transmission-kernel:** a steeper downward slope means co-membership drops off faster with distance, i.e. transmission is more local. Flattening curves over epochs (e.g. Alpha steep → Omicron shallow) are consistent with hypothesis 2 (spatial structure erodes as mixing intensifies). The y-values are calibrated by sampling intensity, so absolute vertical position is less informative than *slope differences* between epochs. Companion table `tables/fig3_distance_decay.csv`: `epoch`, `distance_km_mid` (bin midpoint), `p_co_cluster` (k/n), `ci_low`/`ci_high` (Wilson), `n_pairs`.

### Fig. 4 — Urban vs. rural cluster footprint

Three side-by-side violins per locale (Urban / Mixed / Rural): (A) `log1p(n_sequences)`, (B) `log1p(mean_pairwise_km)`, (C) `log1p(bbox_diag_km)`. The black dot on each violin marks the median. **Expected pattern for hypothesis 3:** Urban violins skew taller in (A) but shorter in (B) and (C); Rural violins are shorter in (A) but taller in (B) / (C). "Mixed" acts as a reference middle. Because locales are defined by access-rank quantiles within *this* dataset, the Urban / Rural tails are roughly equal-N (20% each) — don't read density differences between violins as Scotland-wide population claims.

### Fig. 5 — Size × footprint plane, with deprivation overlay

Two panels, same axes — **x = log(1 + cluster size)**, **y = log(1 + mean pairwise km)**. **Left panel** is a greyscale hexbin density: every cluster contributes one hex count. **Right panel** re-plots only hexes with ≥ 10 clusters, colouring each hex by the **mean modal SIMD quintile** of clusters in that hex (RdBu diverging palette: red = Q1-dominated, blue = Q5-dominated); marker size scales with √n. **Read it for deprivation's position in the size-footprint plane:** if red hexes concentrate in the bottom-right (large, spatially compact) and blue hexes in the upper-left (small, dispersed), that's the Paper 3 deprivation signature — deprived clusters bigger but more geographically concentrated. A uniform blue-red mix means deprivation doesn't sort by size-vs-footprint; the figure is then a null result worth stating in the Discussion.

## Inputs

- `data/processed/scotland_clustering_analysis_dataset.parquet` — needs `dz_xcoord`, `dz_ycoord` (British National Grid, easting/northing in metres).
- Optionally `data/raw/datazone/sg_datazone_bdry_2011.shp` for the map basemap (not required; Fig. 1 falls back to a rough Scotland coastline overlay if shapefile unavailable).

## Running

```bash
python -m manuscripts.paper3_spatial.make_figures --figures manuscripts/paper3_spatial/figures
```

## Statistical models

`models/spatial_mixing.py`:

- **Within-cluster spread** — log mean pairwise distance per cluster vs. `simd_quintile_mode` + VOC, linear regression.
- **Distance-decay** — fit a pair-level logistic regression: `co_cluster ~ log(distance) + voc`.

## Known caveats

1. Coordinates are data-zone centroids, not patient-level. Within-DZ transmission is invisible; footprint of 0 km for singletons or co-DZ clusters.
2. Rural/urban via SIMD geographic-access decile is a proxy only; preferred is the Scottish Government 6-fold rural/urban classification if the TSV can be linked by `datazone`.
3. Distance-decay is sensitive to the sequencing effort and residential bias of COG-UK; results should be presented alongside sampling-intensity maps.
