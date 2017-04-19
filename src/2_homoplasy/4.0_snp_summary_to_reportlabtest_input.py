
'''

Idea is to plot the tip states for homoplasic sites in and around genes against the tree.
'''


import collections
import itertools
import pandas as pd
import sys
import os
sys.path.append(os.path.abspath("/nfs/users/nfs_b/bb9/workspace/rotation3/src/scripts/"))
import parse_ft_tab_file

#  Read in snp state for each taxon
snp_summary_file = '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_snps/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised.out'
snp_summary = pd.read_table(snp_summary_file)
snp_summary.columns = map(str.strip, snp_summary.columns)
# Set ref bases
snp_summary = snp_summary.apply(lambda x: x.replace('.', x['Ref_base']), axis=1)

#  Reportlabtest also allows metadata to be added to trees from a csv file. The
#  csv file must be comma separated, with the first column containing taxon names
#  equivalent to those in the tree. The csv file should also contain a header row
#  containing headings for each column of data. The metadata file is provided with
#  the -m option.

# Create taxa by snps matrix
snps = snp_summary.iloc[:, 10:].transpose()
snps.index.name = 'taxon'
snps.columns = snp_summary['Position_in_CC8_USA300_FPR3757']
snps.columns.name = 'loc'

# Create metadata df
meta = snps
meta['taxon'] = meta.index
#
meta_long = pd.melt(meta, id_vars=['taxon'], value_name='snp')

# Find closest locus for each location
hp_file = '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies.csv'
hp = pd.read_csv(hp_file)
# Check that all sites in hp appear in the summary
assert(all(hp['loc'].isin(snps.columns)))
#  Check changes in snp_summary and hp match
# TODO
#
loc_to_locus_tag = dict(zip(hp['loc'], hp['locus_tag']))
loc_to_homoplasic = dict(zip(hp['loc'], hp['n_homoplasic_acctran'] > 0))
#
meta_long['locus_tag'] = meta_long['loc'].apply(lambda x: loc_to_locus_tag[x])
meta_long['homoplasic'] = meta_long['loc'].apply(lambda x: loc_to_homoplasic[x])

# Get all snps associated with a locus
locus_tag = 'SAUSA300_2565'
meta_long = meta_long[meta_long['locus_tag'].apply(lambda x: locus_tag in str(x))]

# Colour loc based on whether a loc is homoplasic
colours = {
    ('.', False): 'white', 
    ('N', False): 'white', 
    ('A', False): 'lightgreen', 
    ('C', False): 'lightblue', 
    ('G', False): 'gray', 
    ('T', False): 'lightred',
    ('.', True): 'white', 
    ('N', True): 'white', 
    ('A', True): 'green', 
    ('C', True): 'blue', 
    ('G', True): 'black', 
    ('T', True): 'red'
}
colours = {
    ('.', False): 13, 
    ('N', False): 13, 
    ('A', False): 13, 
    ('C', False): 13, 
    ('G', False): 13, 
    ('T', False): 13,
    ('.', True): 1, 
    ('N', True): 13, 
    ('A', True): 3, 
    ('C', True): 4, 
    ('G', True): 14, 
    ('T', True): 2 
}
meta_long.loc[:, 'colour'] = meta_long.copy().apply(lambda x: colours[(x['snp'], x['homoplasic'])], axis=1)

# Generate .tab file features
def meta_row_to_tab_ft(x):
    '''Convert metadata table row to .tab file feature.
    '''
    ft = parse_ft_tab_file.Feature(
            key='SNP',
            loc=x['loc'], 
            qualifs={'taxa': x['taxon'], 'colour': x['colour']}
    )
    return(ft)
#
ft_tab_file = parse_ft_tab_file.FtTabFile()
ft_tab_file.features = list(meta_long.apply(meta_row_to_tab_ft, axis=1))
ft_tab_file.write_tab('test.tab')

cmd = ' '.join((r'python2 /nfs/users/nfs_b/bb9/workspace/rotation3/src/2_homoplasy/reportlabtest_modified.py',
r'-t "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_acctran_steps.tre"',
r'-q taxa',
r'-b {start} -e {end}',
r'-l 2 -A 2 -E 15',
r'-a 2 -L left',
r'-o "test.pdf"',
r'"{tab}" "{embl}"')).format(start=min(2774340, min(meta_long['loc'])), end=max(2777039, max(meta_long['loc'])), tab='test.tab', embl='/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/CC8_USA300_FPR3757.embl')
#
with open('test.sh', 'w') as test_fhandle:
    test_fhandle.write(cmd)

#  FT   SNP             2084488
#  FT                   /node="0->12673_2#33"
#  FT                   /SNP="G->A"
#  FT                   /codon_type="Intergenic"
#  FT                   /colour=1
#  FT                   /taxa="12673_2#33"

# Read in the sites
# Shade sites depending on whether they are mostly convergence or reversal


# Reassign the color values in a .tab file depending on the change
tab_file = '/lustre/scratch118/infgen/team81/bb9/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_homoplasies_on_tree.tab'

Iterate through tabfile
    determine if the record is a terminal branch
        lookup SNP location: determine if the loc is homoplasic
        lookup correct colour for the change


# Determine their closest gene
genes_file =  '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies_per_gene.csv'






