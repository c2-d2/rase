#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

.PHONY: \
        all clean cleanall install \
        test \
        pylint flake8 yapf \
        inc pypi sha256 \
        docs readme wpypi wconda \
        deppip depconda \
        help

PIP=/usr/bin/env pip
PYTHON=/usr/bin/env python3

ROOT_DIR = $(shell pwd)

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p

.SECONDARY:

.SUFFIXES:


###############
# BASIC RULES #
###############


all: ## Run everything

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	make -C pipeline clean

cleanall: ## Clean all files (including bam files and logs)
cleanall: clean
	make -C pipeline cleanall


install: ## Install using PIP
install: hooks
	$(PIP) uninstall -y rase || true
	$(PIP) install .


###########
# TESTING #
###########

test: ## Run tests
	$(MAKE) -C tests clean
	$(MAKE) -C tests

pylint: ## Run PyLint
	$(PYTHON) -m pylint rase

flake8: ## Run Flake8
	flake8

yapf: ## Run YAPF (inline replacement)
	yapf -i --recursive rase setup.py tests scripts

#############
# RELEASING #
#############

inc: ## Increment version
inc: hooks
	./rase/increment_version.py

pypi: ## Upload to PyPI
pypi: hooks
	$(MAKE) clean
	$(PYTHON) setup.py sdist bdist_wheel upload

sha256: ## Compute sha256 for the PyPI package
sha256:
	s=$$(curl https://pypi.python.org/pypi/rase  2>/dev/null| perl -pe 's/#/\n/g' | grep -o 'https.*\.tar\.gz' | xargs curl -L 2>/dev/null | shasum -a 256 | awk '{print $$1;}'); echo $$s; echo $$s | pbcopy


#######################
# DOCUMENTATION & WEB #
#######################

readme: ## Generate readme.html
	markdown_py readme.md > readme.html

wconda: ## Open Bioconda webpage
	open https://bioconda.github.io/recipes/rase/README.html

wpypi: ## Open PyPI webpage
	open https://pypi.python.org/pypi/rase

########################
# INSTALL DEPENDENCIES #
########################

depconda: ## Install dependencies using Conda
	cat requirements.txt | perl -pe 's/==.*//g' | xargs conda install

deppip: ## Install dependencies using PIP
	cat requirements.txt | xargs $(PIP) install

