# Build the screening summary PDF.
#
#   make           generate numbers/tables/figures, then compile
#   make check     page count and leftover LaTeX warnings
#   make clean     remove build output, keep data/
#
# Needs data/per_molecule.csv, data/summary.json and data/redock_control.json
# (optionally data/ligands.json) pulled from the cluster. make_report.py
# prints the exact pull command if they are missing.

# Prefer a local .venv if one exists: this machine's system python3 has
# matplotlib but not pandas, so `make venv` is usually needed once.
PY      ?= $(shell [ -x .venv/bin/python ] && echo .venv/bin/python || echo python3)
DATADIR ?= data

.PHONY: all venv check clean

all: summary.pdf

# The generator must run before latexmk: numbers.tex and tables/ are inputs to
# the document, not products of it. A macro used but never generated shows up
# as an "Undefined control sequence" and -halt-on-error turns that into a
# failed build, which is the check that matters.
summary.pdf: summary.tex make_report.py
	$(PY) make_report.py --data-dir $(DATADIR)
	latexmk -pdf -interaction=nonstopmode -halt-on-error summary.tex

venv:
	python3 -m venv .venv
	.venv/bin/pip install --quiet --upgrade pip
	.venv/bin/pip install --quiet pandas numpy matplotlib
	@echo 'created .venv; make will now use it'

check: summary.pdf
	@echo "--- pages (expect 2-4) ---"
	@pdfinfo summary.pdf | grep -i '^Pages'
	@echo "--- overfull boxes / warnings ---"
	@grep -c -e Overfull -e Underfull summary.log || echo 0
	@grep -i 'warning' summary.log | grep -vi 'hyperref\|microtype' || echo "none"

clean:
	latexmk -C summary.tex 2>/dev/null || true
	rm -rf figures tables numbers.tex summary.pdf
