# Paper 1 — socioeconomic deprivation. Run from repository root.

PY ?= python3
OUT := manuscripts/paper1_socioeconomic/output

.PHONY: all clean figs fig1 fig2 fig3 fig4 fig5 fig6

all: figs

figs:
	$(PY) -m manuscripts.paper1_socioeconomic.make_figures --output $(OUT)

fig1:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig1_sequences_by_simd_over_time --output $(OUT)
fig2:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig2_cluster_size_by_simd --output $(OUT)
fig3:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig3_simd_domain_forest --output $(OUT)
fig4:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig4_singleton_odds_by_epoch --output $(OUT)
fig5:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig5_deprivation_lineage_heatmap --output $(OUT)
fig6:
	$(PY) -m manuscripts.paper1_socioeconomic.figures.fig6_domain_decomposition --output $(OUT)

clean:
	rm -f $(OUT)/*.pdf $(OUT)/*.png
	rm -f manuscripts/paper1_socioeconomic/tables/*.csv
