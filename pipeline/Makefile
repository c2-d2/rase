#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

.PHONY: all help clean cleanall cluster export clock pacbio

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p

.SECONDARY:

.SUFFIXES:

all: ## Run everything
	$(SM)

clock: ## Run everything in the clock mode
	$(SM) --configfile config.clock-mode.yaml

pacbio: ## Run everything in the pacbio mode
	$(SM) --configfile config.pacbio.yaml

test: ## Run the smallest experiment only
	ls database/*.tsv || \
		(cd database \
		&& wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tar.gz \
		&& wget https://github.com/c2-d2/rase-db/releases/download/v01/spneumoniae_sparc.k18.tsv)
	ls reads/*.fq || \
		(cd reads \
		&& wget https://zenodo.org/record/1405173/files/sp10_norwich_P33.filtered.fq)
	$(SM) test

replot: ## Replot all figures
	rm plots/*.pdf
	$(SM)

export: ## Export all outputs
	f="rase_results.$$d.tar"; \
	rm -f $$f; \
	tar cvf $$f --files-from /dev/null; \
	for x in prediction/*.tsv plots/*.pdf; do \
		echo $$x; \
		tar --append --file=$$f $$x; \
	done; \
	for x in prediction/*/; do \
		y=$$(ls $${x}*.tsv | sort | tail -n 1) || continue; \
		echo $$y; \
		tar --append --file=$$f $$y; \
	done

cluster: ## Submit jobs to a cluster
	snakemake --cores 9999 -p --keep-going \
		--cluster-config cluster.json \
		--cluster 'sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem={cluster.memory}'

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -vf plots/*.pdf benchmarks/*.plot.log plots/.*.plot{,_timeline,_snapshots}.complete
	rm -vf benchmarks/*.predict.log prediction/*.predict.tsv prediction/.*.predict.complete
	find prediction -name '*.tsv' | xargs rm -vf
	rm -vf exported.*.tar

cleanall: ## Clean all files (including bam files and logs)
cleanall: clean
	rm -f matching/*.{fq,fq.complete} benchmarks/*.readprep.log
	rm -f matching/*.bam
	rm -f prediction/.*.complete
	rm -f benchmarks/*.classify.log
	rm -fr database/.*.complete $$(ls -d database/*/ 2>/dev/null || true) benchmarks/decompress.log

