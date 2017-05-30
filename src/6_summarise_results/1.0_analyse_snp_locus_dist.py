'''Fisher's test for excess of intergenic hp SNPs
'''

import pandas as pd
import scipy.stats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import operator as op
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
    #  test both at the one st and multist levels
    #  including and excluding known amrs
    cont_table = [
        [sum(hp['is_hp'] & hp['is_intergenic']), sum(~hp['is_hp'] & hp['is_intergenic'])],
        [sum(hp['is_hp'] & ~hp['is_intergenic']), sum(~hp['is_hp'] & ~hp['is_intergenic'])]
    ]
    print('===== {} ====='.format(prefix))
    print('is_hp')
    print('[[hi, ni], [hg, ng]]')
    print(cont_table)
    print(scipy.stats.fisher_exact(cont_table))
    cont_table = [
        [sum(hp['is_hp_multi'] & hp['is_intergenic']), sum(~hp['is_hp_multi'] & hp['is_intergenic'])],
        [sum(hp['is_hp_multi'] & ~hp['is_intergenic']), sum(~hp['is_hp_multi'] & ~hp['is_intergenic'])]
    ]
    print('is_hp_multi')
    print('[[hi, ni], [hg, ng]]')
    print(cont_table)
    print(scipy.stats.fisher_exact(cont_table))

hp['dist_to_locus'] = hp['intergenic'].map(op.itemgetter(0))

hp['dist_to_locus'].replace(0, 1e-12)

hp[hp.is_hp]['intergenic'].map(op.itemgetter(0)).hist()

hp.is_hp = hp.is_hp.astype('category')
hp.is_intergenic = hp.is_intergenic.astype('category')

tips = sns.load_dataset("tips")
g = sns.FacetGrid(hp, row="is_hp", col="is_intergenic", margin_titles=True)
bins = np.linspace(-750, 750, 150)
g.map(plt.hist, "dist_to_locus", color="steelblue", bins=bins, lw=0)
g.set(ylim=(None, 100))

