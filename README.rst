==================================================================
Calling *Platygyra daedalea* parental genotypes from mapping files
==================================================================

Preamble
--------
This repository contains the scripts written to parse the large ``*.bam`` files produced from mapping into intermediate ``*.tsv`` files, to produce a compiled table (``pdae.all_parent_scores.tsv.gz``).

The first script
----------------
``*.bam`` files produced from mapping were first converted to mpileup tables with the following command (``samtools mpileup`` with default settings)::

  for a in *.sorted.bam; do samtools mpileup ${a} | gzip > ${a/sorted.bam/sorted.tsv.gz}; done

(the bam files and the compressed mpileups are massive, hence not uploaded.)

``parse_mpileup_tsvs.py`` then reads the fifth column of the ``mpileup`` outputs, and parses the per-base composition at every position. The parsed tables are provided here, in the form ``???.basecomp.tsv.gz``. A short excerpt from an example file::

  zcat F73.basecomp.tsv.gz | head -2
  
  scaffold1|size5511861   6544    6       0       0       0
  scaffold1|size5511861   6545    0       0       6       0

For scaffold1 position 6544, there are 6 A, 0 C, 0 G, and 0 T; for scaffold1 position 6545, there are 0 A, 0 C, 6 G, and 0 T.

The second script
-----------------
Using the ``*.basecomp.tsv.gz`` files and the *Platygyra daedalea* genome, ``score_all_snps.py`` checks whether the per-base composition resembles that of the homozygous reference, heterozygous, or homozygous alt; numeric scores are given to facilitate downstream correlation/heatmap plotting. These calls are subject to coverage cutoffs, with position with < 5 coverage deemed NA.

Genotype calls across all files are then merged into a single table (provided in this repo, ``pdae.all_parent_scores.tsv.gz``).
