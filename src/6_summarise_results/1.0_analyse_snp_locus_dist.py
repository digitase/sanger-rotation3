'''Fisher's test for excess of intergenic hp SNPs
'''

import pandas as pd
import scipy.stats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import operator as op
import sys
plt.ion()

# Read in SNPs homoplasic in multiple sts
summary_sites_by_gene = pd.read_csv(
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_annotate_candidates/hp_sites_summary.known_amr_marked.csv', 
    converters={'sts': eval, 'loc_in_st': eval, 'known_amr_loci': eval, 'intergenic': eval}
)
# Exclude MGEs, repetitive proteins
#
#  phage terminase family protein, transposase
#  spa (immunoglobulin G binding protein A precursor), ebh, clfA, clfB
summary_sites_by_gene = summary_sites_by_gene.loc[ 
    ~summary_sites_by_gene.apply(lambda x: 
        "'ebh'" in x['genes'] or
        "'spa'" in x['genes'] or
        "'clfA'" in x['genes'] or
        "'clfB'" in x['genes'] or
        "'phage terminase family protein'" in x['products'] or
        "transposase" in x['products']
    , axis=1)
]
#  NOTE: annotation is incomplete, not all known amr loci are identified 
summary_sites_by_gene['is_known_amr_locus'] = summary_sites_by_gene['known_amr_loci'].map(lambda sts: not all(len(x) == 0 for x in sts))
#  summary_sites_by_gene = summary_sites_by_gene[~summary_sites_by_gene['is_known_amr_locus']]
#
# Locations of all SNPs
hp_files = {
    'st22': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/st22_homoplasies.csv', 
    'st239': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/st239_homoplasies.csv', 
    #  'st30': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/st30_homoplasies.csv', 
    'st8': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies.csv' 
}
# without known amr sites
#  hp_files = {
    #  'st22': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st22_homoplasies.known_amr_filtered.csv', 
    #  'st239': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st239_homoplasies.known_amr_filtered.csv', 
    #  'st8': '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st8_homoplasies.known_amr_filtered.csv'
#  }

hp_merged = []
for prefix in hp_files.keys():
    #  prefix = 'st22'
    hp = pd.read_csv( 
        hp_files[prefix], 
        converters={'intergenic': eval, 'locus_tag': eval}
    )
    hp['is_intergenic'] = hp['intergenic'].map(lambda x: x != (0, ))
    hp['is_hp'] = hp['n_homoplasic_acctran'].astype(bool) | hp['n_homoplasic_deltran'].astype(bool)
    summary_sites_by_gene_st = summary_sites_by_gene.loc[summary_sites_by_gene['sts'].map(lambda x: prefix in x), ]
    locs_in_st = set()
    for sts, locs in zip(summary_sites_by_gene_st['sts'], summary_sites_by_gene_st['loc_in_st']):
        for st, loc in zip(sts, locs):
            if prefix == st:
                locs_in_st.add(loc)
    hp['is_hp_multi'] = hp['loc'].isin(locs_in_st)
    hp_merged.append(hp)
    # Per strain tests
    #
    #  cont_table = [
        #  [sum(hp['is_hp'] & hp['is_intergenic']), sum(~hp['is_hp'] & hp['is_intergenic'])],
        #  [sum(hp['is_hp'] & ~hp['is_intergenic']), sum(~hp['is_hp'] & ~hp['is_intergenic'])]
    #  ]
    #  print('===== {} ====='.format(prefix))
    #  print('is_hp')
    #  print('[[hi, ni], [hg, ng]]')
    #  print(cont_table)
    #  print(scipy.stats.fisher_exact(cont_table))
    #  cont_table = [
        #  [sum(hp['is_hp_multi'] & hp['is_intergenic']), sum(~hp['is_hp_multi'] & hp['is_intergenic'])],
        #  [sum(hp['is_hp_multi'] & ~hp['is_intergenic']), sum(~hp['is_hp_multi'] & ~hp['is_intergenic'])]
    #  ]
    #  print('is_hp_multi')
    #  print('[[hi, ni], [hg, ng]]')
    #  print(cont_table)
    #  print(scipy.stats.fisher_exact(cont_table))
hp_merged = pd.concat(hp_merged)
hp_merged['dist_to_locus'] = hp_merged['intergenic'].map(op.itemgetter(0))

# Plotting dist to locus
plt.figure(figsize=(11, 6))
sns.distplot( 
    hp_merged.loc[~hp_merged.is_hp, 'dist_to_locus'],
    kde=False, bins=150, label='Non-homoplasic', hist_kws={"alpha": 1, "color": "grey"}
)
sns.distplot( 
    hp_merged.loc[hp_merged.is_hp, 'dist_to_locus'],
    kde=False, bins=150, label='Homoplasic in >=1 ST', hist_kws={"alpha": 1, "color": "lightblue"}
)
sns.distplot( 
    hp_merged.loc[hp_merged.is_hp_multi, 'dist_to_locus'],
    kde=False, bins=150, label='Homoplasic in >1 ST', hist_kws={"alpha": 1, "color": "darkblue"}
)
plt.yscale('log')
plt.ylabel('Frequency')
plt.xlabel('Location relative to nearest locus')
plt.legend()
plt.savefig('.output/analyse_snp_locus_dist.pdf')

#  Are there more hps than expected in intergenic regions?
#
# log output
orig_stdout = sys.stdout
sys.stdout = open('.output/analyse_snp_locus_dist.results.txt','w')
#
cont_table = pd.crosstab(hp_merged['is_hp'], hp_merged['is_intergenic'])
print(cont_table)
print(scipy.stats.fisher_exact(cont_table))
#
cont_table = pd.crosstab(hp_merged['is_hp_multi'], hp_merged['is_intergenic'])
print(cont_table)
print(scipy.stats.fisher_exact(cont_table))

#  Are there more hps than expected in a certain region upstream of genes?
hp_merged['in_region'] = ( hp_merged['dist_to_locus'] < 0 ) & (hp_merged['dist_to_locus'] >= -150)
#
cont_table = pd.crosstab(hp_merged.loc[hp_merged.is_intergenic, 'is_hp_multi'], hp_merged.loc[hp_merged.is_intergenic, 'in_region'])
print(cont_table)
print(scipy.stats.fisher_exact(cont_table))
#
cont_table = pd.crosstab(hp_merged.loc[hp_merged.is_intergenic, 'is_hp'], hp_merged.loc[hp_merged.is_intergenic, 'in_region'])
print(cont_table)
print(scipy.stats.fisher_exact(cont_table))

#
sys.stdout.close()
sys.stdout = orig_stdout

#  hp_merged.is_hp = hp_merged.is_hp.astype('category')
#  hp_merged.is_intergenic = hp_merged.is_intergenic.astype('category')
#
#  g = sns.FacetGrid(hp_merged, row="is_hp", margin_titles=True)
#  bins = np.linspace(-750, 750, 150)
#  g.map(plt.hist, "dist_to_locus", color="steelblue", bins=bins, lw=0)
#  g.set(yscale='log')

