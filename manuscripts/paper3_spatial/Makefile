# Paper 3 — spatial structure. Run from repository root.

PY ?= python3
OUT := manuscripts/paper3_spatial/output

.PHONY: all figs clean fig1 fig2 fig3 fig4 fig5

all: figs

figs:
	$(PY) -m manuscripts.paper3_spatial.make_figures --output $(OUT)

fig1:
	$(PY) -m manuscripts.paper3_spatial.figures.fig1_cluster_centroid_map --output $(OUT)
fig2:
	$(PY) -m manuscripts.paper3_spatial.figures.fig2_within_cluster_spread --output $(OUT)
fig3:
	$(PY) -m manuscripts.paper3_spatial.figures.fig3_distance_decay --output $(OUT)
fig4:
	$(PY) -m manuscripts.paper3_spatial.figures.fig4_urban_rural_footprint --output $(OUT)
fig5:
	$(PY) -m manuscripts.paper3_spatial.figures.fig5_spread_by_deprivation --output $(OUT)

clean:
	rm -f $(OUT)/*.pdf $(OUT)/*.png
	rm -f manuscripts/paper3_spatial/tables/*.csv
