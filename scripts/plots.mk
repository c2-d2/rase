#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

.PHONY: all help clean cleanall merge
.SECONDARY:
.SUFFIXES:

SHELL=/usr/bin/env bash -eo pipefail

include ../conf.mk

all: .plot.complete .stats.complete .merge.mic.complete .merge.res_tree.complete

.plot.complete: $(foreach ant,$(antibiotics),.plot.$(ant).complete)

.stats.complete:
	../../scripts/compute_stats.py ../resistance/res_cat.tsv > $(name).phylogroups.tsv
	touch $@

.plot.%.complete:
	../../scripts/plot_mic.R ../resistance/res_cat.tsv -b $($*) $* res_cat.$*.basic.pdf
	../../scripts/reproducible_pdf.sh                              res_cat.$*.basic.pdf
	../../scripts/plot_mic.R -z ../resistance/res_cat.tsv -b $($*) $* res_cat.$*.zoomed.pdf
	../../scripts/reproducible_pdf.sh                                 res_cat.$*.zoomed.pdf
	touch $@

.merge.mic.complete: .plot.complete
	/usr/bin/env gs \
		-sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default \
		-dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 \
		-dAutoRotatePages=/None \
		-sOutputFile=$(name).mic.pdf \
		$(foreach ant,$(antibiotics),res_cat.$(ant).zoomed.pdf)

	../../scripts/reproducible_pdf.sh $(name).mic.pdf
	touch $@

.merge.res_tree.complete: $(wildcard ../resistance/_res_tree.*.pdf)
	/usr/bin/env gs \
		-sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default \
		-dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 \
		-dAutoRotatePages=/None \
		-sOutputFile=$(name).res_tree.pdf \
		$(foreach ant,$(antibiotics),../resistance/_res_tree.$(ant).pdf)

	../../scripts/reproducible_pdf.sh $(name).res_tree.pdf
	touch $@

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f res_cat.*.pdf .*.complete
	rm -f res_tree.*.pdf

cleanall: clean ## Clean
	rm -f *.pdf *.tsv
