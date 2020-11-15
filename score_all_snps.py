#!/usr/bin/env python3

"""
> score_all_snps.py <

Compare all parental genotypes, see which samples are more closely related
to which others.
"""
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

import parse_fasta

base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def get_dist_from_ref_base(scaf, pos, ref_base, A, C, G, T):
    """
    Based on the counts of A, C, G and T at a specific genomic loci, check
    what the genomic reference is, and calculate distance from reference.
    """
    # exit function with distance 0 if ref base isn't A/C/G/T
    if ref_base not in ['A', 'C', 'G', 'T']: return 0
    
    basecomp = np.array([A, C, G, T], np.int64)
    
    # institute minimum coverage filter--if less than 5, return NaN (i.e.
    # we don't know the distance to ref base for that low-covered position)
    if np.sum(basecomp) < 5:
        return np.NaN
    
    # zero dictionary entry corresponding to the reference base, so that 
    # we can simply equate distance from reference by counting non-zeros
    basecomp[base_to_int[ref_base]] = 0
    
    # distance: ref/ref = 0; ref/alt = 1; alt/alt = 2
    # some positions have three base calls--prevent max distance from being
    # greater than 2
    dist_from_ref_base = min(np.count_nonzero(basecomp), 2)
    
    return dist_from_ref_base

# hardcode important raw data files
PDAE_GENOME_FILE = Path('../raw_data/pdae_genome.v1.fa.gz')
pdae_genome = parse_fasta.get_all_sequences(PDAE_GENOME_FILE, 'fasta')

# use a pandas dataframe to hold the compiled information
compiled_dists = pd.DataFrame()

# iterate through *.basecomp.tsv.gz files
parent_tsvs = Path('./').glob('*basecomp.tsv.gz')
for pt in parent_tsvs:
    parent_label = pt.name.split('.')[0]
    print (f'[{time.asctime()}] Processing file {pt.name}...', file=sys.stderr)
    
    # process per-parent basecomp file
    parent_basecomp = pd.read_table(pt, names=['scaf', 'pos', 'A', 'C', 'G', 'T'])
    # remember that `pos` values are 1-based, but arrays are 0-based
    parent_basecomp['ref_base'] = parent_basecomp.apply(
        lambda x: pdae_genome[x['scaf']][x['pos'] - 1], axis=1)
    parent_basecomp[parent_label] = parent_basecomp.apply(
        lambda x: get_dist_from_ref_base(x['scaf'], x['pos'], x['ref_base'], 
                                         x['A'], x['C'], x['G'], x['T']),
        axis=1)
    
    # merge calculated distances into the compiled dataframe
    parent_basecomp = parent_basecomp.drop(columns=['A', 'C', 'G', 'T'])
    if compiled_dists.empty:
        compiled_dists = parent_basecomp
    else:
        compiled_dists = compiled_dists.merge(
            parent_basecomp, how='outer', on=['scaf', 'pos', 'ref_base'])

# save table
compiled_dists.to_csv('pdae.all_parent_scores.tsv', sep='\t', na_rep='NA', index=False)
