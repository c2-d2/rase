.PHONY: all help clean cleanall o2

SHELL=/usr/bin/env bash -eo pipefail

SM=snakemake -j -p

.SECONDARY:

.SUFFIXES:

all:
	$(SM)

replot:
	rm plots/*.pdf
	$(SM)

o2:
	snakemake --cores 9999 -p \
		--cluster-config cluster.o2.json \
		--cluster 'sbatch -p {cluster.queue} -c {cluster.n} -t {cluster.time} --mem={cluster.memory}'

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f plots/*.pdf benchmarks/*.plot.log
	rm -f prediction/*.quantify.complete benchmarks/*.quantify.log
	rm -f prediction/*.predict.{tsv,complete} prediction/*/*.tsv benchmarks/*.predict.log

cleanall: clean
	rm -f prediction/*.{fq,fq.complete} benchmarks/*.readprep.log
	rm -f prediction/*.{bam,bam.complete} benchmarks/*.classify.log
	rm -fr database/*.complete $$(ls -d database/*/ || true) benchmarks/decompress.log

