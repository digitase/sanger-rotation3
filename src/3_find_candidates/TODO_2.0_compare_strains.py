'''Compare homoplasic features across strains.
'''

import pandas as pd

prefixes = {
    "st22": "S.aureus_ST22_BSAC_Pfizer",
    "st239": "S.aureus_ST239_global_Singapore_Pfizer",
    "st30": "S.aureus_ST30_BSAC_Pfizer",
    "st8": "S.aureus_ST8_BSAC_Pfizer_revised"
}
snp_summary_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_snps/{prefix}/{prefix}.out"
tree_files = {
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/RAxML_bestTree.ml_S.aureus_ST22_BSAC_Pfizer",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/RAxML_bestTree.ml_S.aureus_ST239_global_Singapore_Pfizer",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/RAxML_bestTree.ml_S.aureus_ST30_BSAC_Pfizer",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/revised_lane_list/RAxML_bestTree.ml_S.aureus_ST8_BSAC_Pfizer_revised"
}
genes_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/{prefix}/{short_prefix}_homoplasies_per_gene.csv"
embl_files = {
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}

# Read in gene homoplasies
genes = dict((short_prefix, pd.read_csv(genes_file_template.format(prefix=prefixes[short_prefix], short_prefix=short_prefix))) for short_prefix in prefixes)

genes_ranked = {}
for short_prefix, gene in genes.items():
    # Rank genes by: (agree_prop_acctran_deltran * n_homoplasic_acctran)/n_total
    gene = genes[short_prefix]
    gene.loc[:, 'ranker'] = gene['agree_prop_acctran_deltran'] * gene['n_homoplasic_acctran'] / gene['n_total']
    min_hp_thresh = gene.loc[gene.ranker > 0, 'n_homoplasic_acctran'].quantile(0.90)
    gene = gene.loc[(gene['ranker'] > 0) & (gene['n_homoplasic_acctran'] > min_hp_thresh)].copy()
    gene['rank'] = gene['ranker'].rank(method='average')
    genes_ranked[short_prefix] = gene

#  Sum ranks over all strains
genes_ranked_merged = pd.concat(
    genes_ranked.values(), axis=1
)






