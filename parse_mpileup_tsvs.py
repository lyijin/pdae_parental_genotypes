#!/usr/bin/env python3

"""
> parse_mpileup_tsvs.py <

Script parses mpileup tsvs, collates files belonging to same parent, and outputs
  1. scaffold name
  2. coord
  3. base composition at coordinate
    a. number of A
    b. number of C
    c. number of G
    d. number of T

To counter spurious SNPs arising from sequencing errors, base counts that are
< 5% of coverage would be zeroed out. e.g.
  (F77) scaffold9|size2634320  978066  1  0  1214  0  -->  0  0  1214  0
"""
import csv
import gzip
from pathlib import Path
import re
import sys
import time

import natural_sort

# hardcode folder path containing mpileup tsvs
MPILEUP_FOLDER = Path('../raw_data/mpileup_tsvs')

mpileup_tsvs = MPILEUP_FOLDER.glob('*.tsv.gz')
unique_parents = sorted(list(set([x.name[:4] for x in mpileup_tsvs])))

# iterate through all files per unique parent
for u in unique_parents:
    u_tsvs = MPILEUP_FOLDER.glob(f'{u}*.tsv.gz')
    compiled_coverage = {}
    lineno = 0
    
    for ut in u_tsvs:
        print (f'[{time.asctime()}] Parsing {str(ut)}...', file=sys.stderr)
        
        tsv_reader = csv.reader(gzip.open(str(ut), 'rt'), delimiter='\t')
        for row in tsv_reader:
            scaf, start, _, _ = row[0].split('_')
            offset = row[1]
            pos = int(start) + int(offset) - 1
            bases = row[4]
            
            # mpileup adds additional text to denote "start of read" (^.) and
            # "end of read" ($). for "start of read", the ascii character
            # following caret is the Phred+33 quality of the base.
            # strip these text to retain sequenced bases
            bases = bases.replace('$', '')
            bases = re.sub(r'\^.', '', bases)
            
            scaf_pos = f'{scaf}_{pos}'
            if scaf_pos in compiled_coverage:
                compiled_coverage[scaf_pos] += bases.upper()
            else:
                compiled_coverage[scaf_pos] = bases.upper()
            
            # indicate progress
            lineno += 1
            if lineno % 1_000_000 == 0:
                print (f'[{time.asctime()}] Reading line {lineno:,}...',
                       file=sys.stderr)
    
    # print compiled coverages out
    outfile = gzip.open(f'{u.replace("_", "")}.basecomp.tsv.gz', 'wt')
    for scaf_pos in natural_sort.natural_sort(compiled_coverage):
        n_a = compiled_coverage[scaf_pos].count('A')
        n_c = compiled_coverage[scaf_pos].count('C')
        n_g = compiled_coverage[scaf_pos].count('G')
        n_t = compiled_coverage[scaf_pos].count('T')
        n_acgt = n_a + n_c + n_g + n_t
        
        # zero out anything not > 5% * n_acgt
        n_a = n_a if n_a > 0.05 * n_acgt else 0
        n_c = n_c if n_c > 0.05 * n_acgt else 0
        n_g = n_g if n_g > 0.05 * n_acgt else 0
        n_t = n_t if n_t > 0.05 * n_acgt else 0
        
        print (*scaf_pos.split('_'), n_a, n_c, n_g, n_t, sep='\t', file=outfile)
    outfile.close()
