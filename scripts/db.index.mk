.PHONY: all help clean cleanall
.SECONDARY:
.SUFFIXES:

SHELL=/usr/bin/env bash -eo pipefail

include ../conf.mk

all: $(foreach onek,$(k),.compress.$(onek).complete)

.compress.%.complete: .index.%.complete
	prophyle compress $(name).k$* $(name).k$*.tar.gz
	touch $@

.index.%.complete:
	prophyle index -k $* -A -g ../isolates ../tree/tree.newick $(name).k$*
	touch $@

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -fr $$(ls -d */) *.tar.gz .*.complete

cleanall: clean ## Cleanall
