"""This Snakefile serves for submitting the index construction to cluster. This
can be useful especially for a parallel construction of the same RASE DB
with different k-mer lenghts. To be symlinked.

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

rule all:
    input:
        #". .complete"
        [".k={}.complete".format(k) for k in range(10, 33)]

rule index:
    output: touch(".{params}.complete")
    shell:
        """
            make {wildcards.params}
        """
