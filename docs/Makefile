SPHINXBUILD ?= python -m sphinx
.PHONY: html
html:
	$(SPHINXBUILD) -n -W --keep-going -b html . _build/html
	python ../scripts/check_docs.py _build/html
