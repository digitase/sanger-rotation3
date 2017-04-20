'''Compare homoplasic features across strains.
'''

import pandas as pd
import collections
import os
import itertools
import networkx as nx

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
    gene['rank_f'] = gene['rank'] / len(gene['rank'])
    # TODO convert to percentile rank as number of ranked genes differ between strains
    genes_ranked[short_prefix] = gene

# Read in reciprocal best hits tables and build adjaceny graph
blast_rbh_dir = '/nfs/users/nfs_b/bb9/workspace/rotation3/src/3_find_candidates/.output/blast/blast_rbh/'
rbh_tables = [f for f in os.listdir(blast_rbh_dir) if f.endswith('.out')]
rbh_edges = []
for table in rbh_tables:
    sp1, sp2 = table.split('.')[1].split('_')
    rbh = pd.read_table(os.path.join(blast_rbh_dir, table))
    # Assert rbhs are unique
    assert(len(rbh['#A_id']) == len(rbh['#A_id'].unique()))
    assert(len(rbh['B_id']) == len(rbh['B_id'].unique()))
    rbh['locus_tag_sp1'] = rbh['#A_id'].map(lambda x: '_'.join(x.split('_')[:-3]))
    rbh['locus_tag_sp2'] = rbh['B_id'].map(lambda x: '_'.join(x.split('_')[:-3]))
    for lt1, lt2 in zip(rbh['locus_tag_sp1'], rbh['locus_tag_sp2']):
        rbh_edges.append((lt1, lt2))
G = nx.from_edgelist(rbh_edges).to_undirected()

#
#  Collate ranks over all strains
#
#  For groups of homologous loci,
#  obtain mean rank (normalised to 0-1, non-homoplasic sites ignored) over strains,
#  ranking by (agree_prop_acctran_deltran * n_homoplasic_acctran)/n_total
#
ranked_merged = pd.concat(genes_ranked.values(), names=genes_ranked.keys())
ranked_merged['locus_tag_simplified'] = ranked_merged['locus_tag'].map(lambda x: x.split("'")[1])
# Process genic and intergenic regions separately
ranked_merged_intergenic = ranked_merged[ranked_merged['intergenic']]
ranked_merged_genic = ranked_merged[~ranked_merged['intergenic']]
#
homolog_genic = []
homolog_ranks_genic = []
for c in nx.connected_components(G):
    homolog_genic.append(c)
    homolog_ranks_genic.append(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'rank_f'].mean())
#
homolog_intergenic = []
homolog_ranks_intergenic = []
for c in nx.connected_components(G):
    homolog_intergenic.append(c)
    homolog_ranks_intergenic.append(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'rank_f'].mean())

pd.DataFrame({
    'homologs': homolog_genic, 
    'rank_f': homolog_ranks_genic, 
}).sort_values('rank_f', ascending=False).head(25)

pd.DataFrame({
    'homologs': homolog_intergenic, 
    'rank_f': homolog_ranks_intergenic, 
}).sort_values('rank_f', ascending=False).head(25)

#  TODO: filter out known loci from 2010 paper

