#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

.PHONY: all help clean cleanall
.SECONDARY:
.SUFFIXES:

SHELL=/usr/bin/env bash -eo pipefail

include ../conf.mk

all: $(foreach onek,$(k),.$(onek).complete)

.%.complete:
	cp ../resistance/res_cat.tsv $(name).k$*.tsv
	cp ../index/$(name).k$*.tar.gz $(name).k$*.tar.gz
	touch $@

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f .*.complete *.tar.gz *.tsv

cleanall: clean
