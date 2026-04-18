# Paper 2 — demographic and vaccination correlates. Run from repository root.

PY ?= python3
OUT := manuscripts/paper2_demographic/output

.PHONY: all figs clean fig1 fig2 fig3 fig4 fig5 fig6

all: figs

figs:
	$(PY) -m manuscripts.paper2_demographic.make_figures --output $(OUT)

fig1:
	$(PY) -m manuscripts.paper2_demographic.figures.fig1_age_over_time_by_epoch --output $(OUT)
fig2:
	$(PY) -m manuscripts.paper2_demographic.figures.fig2_vaccination_vs_cluster_size --output $(OUT)
fig3:
	$(PY) -m manuscripts.paper2_demographic.figures.fig3_age_homogeneity --output $(OUT)
fig4:
	$(PY) -m manuscripts.paper2_demographic.figures.fig4_sex_composition --output $(OUT)
fig5:
	$(PY) -m manuscripts.paper2_demographic.figures.fig5_demographic_forest --output $(OUT)
fig6:
	$(PY) -m manuscripts.paper2_demographic.figures.fig6_voc_stratified_shifts --output $(OUT)

clean:
	rm -f $(OUT)/*.pdf $(OUT)/*.png
	rm -f manuscripts/paper2_demographic/tables/*.csv
