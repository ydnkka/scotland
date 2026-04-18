# Top-level Makefile for all three manuscripts. Run from repository root.

PY ?= python3

.PHONY: all paper1 paper2 paper3 clean

all: paper1 paper2 paper3

paper1:
	$(MAKE) -C manuscripts/paper1_socioeconomic figs

paper2:
	$(MAKE) -C manuscripts/paper2_demographic figs

paper3:
	$(MAKE) -C manuscripts/paper3_spatial figs

clean:
	$(MAKE) -C manuscripts/paper1_socioeconomic clean
	$(MAKE) -C manuscripts/paper2_demographic clean
	$(MAKE) -C manuscripts/paper3_spatial clean
