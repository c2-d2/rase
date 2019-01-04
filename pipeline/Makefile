#
# Author:  Karel Brinda <kbrinda@hsph.harvard.edu>
#
# License: MIT
#

.PHONY: all help clean cleanall cluster export

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p

.SECONDARY:

.SUFFIXES:

all: ## Run everything
	$(SM)

test: ## Run the smallest experiment
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
	d=$$(date +"%Y-%m-%d_%H-%M"); \
	f="rase_exported.$$d.tar"; \
	rm -f $$f; \
	tar cvf $$f --files-from /dev/null; \
	for x in benchmarks/*.log prediction/*.tsv plots/*.pdf; do \
		tar --append --file=$$f $$x; \
	done

cluster: ## Submit jobs to a cluster
	snakemake --cores 9999 -p \
		--cluster-config cluster.json \
		--cluster 'sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem={cluster.memory} --job-name {cluster.name} -o {cluster.output} -e {cluster.error}'

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f plots/*.pdf benchmarks/*.plot.log
	rm -f prediction/.*.quantify.complete benchmarks/*.quantify.log
	rm -f prediction/*.predict.tsv
	rm -f prediction/*.predict.complete
	find prediction -name '*.tsv' | xargs rm -f
	rm -f benchmarks/*.predict.log
	rm -f exported.*.tar

cleanall: ## Clean all files (including bam files and logs)
cleanall: clean
	rm -f prediction/*.{fq,fq.complete} benchmarks/*.readprep.log
	rm -f prediction/*.bam
	rm -f prediction/.*.complete
	rm -f benchmarks/*.classify.log
	rm -fr database/.*.complete $$(ls -d database/*/ 2>/dev/null || true) benchmarks/decompress.log

