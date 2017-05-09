'''Parse and merge various sources with known AMR loci
'''

import pandas as pd

sources = {
    'sh_2010': '/nfs/users/nfs_b/bb9/workspace/rotation3/misc/homoplasies_st239_harris_et_al_2010_science_table1.annotated.csv',
    'ly_xls': '/nfs/users/nfs_b/bb9/workspace/rotation3/misc/amr_loci/Resistance mechanisms aureus_mh.xls',
    'lit_review': '/nfs/users/nfs_b/bb9/workspace/rotation3/misc/amr_loci/sa_amr_lit_review.xlsx'
}

amr_sh = pd.read_csv(sources['sh_2010'])
amr_ly = pd.read_excel(sources['ly_xls'])
amr_lit = pd.read_excel(sources['lit_review'])

#  TODO Make sure MGE .tab files we have are excluded.

# Custom filtering for each data source

amr_sh_genes = set(amr_sh.loc[(amr_sh['Antibiotic'] != '-') & (~pd.isnull(amr_sh['Antibiotic'])), 'gene'].map(lambda x: eval(x)[0][0]))

amr_ly_genes = list(g for g in amr_ly.loc[~pd.isnull(amr_ly['Gene']), 'Gene'] if not g in ('-', 'Not assigned'))
# Split on , and ; then add the individual entries
for i, g in enumerate(amr_ly_genes):
    if ';' in g:
        amr_ly_genes[i:i+1] = [x.strip() for x in g.split(';')]
    elif ',' in g:
        amr_ly_genes[i:i+1] = [x.strip() for x in g.split(',')]
# Take first component of the name
for i, g in enumerate(amr_ly_genes):
    amr_ly_genes[i] = g.split(' ')[0].strip()
amr_ly_genes = set(amr_ly_genes)

amr_lit_genes = list(amr_lit.loc[~pd.isnull(amr_lit['gene_symbol']), 'gene_symbol'])
for i, g in enumerate(amr_lit_genes):
    amr_lit_genes[i] = g.split(' ')[0].strip()
# Parse alternatives
for i, g in enumerate(amr_lit_genes):
    if '/' in g:
        if not g.count('/') == 1:
            raise ValueError(g)
        prefix, alt_suffix = g.split('/')
        amr_lit_genes[i:i+1] = [prefix, prefix[: -len(alt_suffix)] + alt_suffix]
amr_lit_genes = set(amr_lit_genes)

# Merge and write
all_genes = amr_sh_genes | amr_ly_genes | amr_lit_genes

with open('.output/known_amr_loci.txt', 'w') as outfile:
    outfile.write('\n'.join(all_genes))

