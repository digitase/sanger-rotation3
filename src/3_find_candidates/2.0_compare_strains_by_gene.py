'''Compare homoplasic features across strains.
'''

from Bio import SeqIO
import pandas as pd
import collections
import os
import itertools
import networkx as nx

prefixes = {
    "st22": "S.aureus_ST22_BSAC_Pfizer",
    "st239": "S.aureus_ST239_global_Singapore_Pfizer",
    #  "st30": "S.aureus_ST30_BSAC_Pfizer",
    "st8": "S.aureus_ST8_BSAC_Pfizer_revised"
}
genes_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/{prefix}/{short_prefix}_homoplasies_per_gene.csv"
embl_files = {
    #  "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/misc/2017-05-23_staph_annot_from_matt/EMRSA15_with_Mels_currated.modified.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}

# Read in gene homoplasies
genes = dict((short_prefix, pd.read_csv(genes_file_template.format(prefix=prefixes[short_prefix], short_prefix=short_prefix))) for short_prefix in prefixes)
# Rank strains separately
print('Ranking loci...')
genes_ranked = {}
for short_prefix, gene in genes.items():
    # Rank genes by: (agree_prop_acctran_deltran * n_homoplasic_sites)/n_sites
    gene = genes[short_prefix]
    gene.loc[:, 'ranker'] = gene['agree_prop_acctran_deltran'] * gene['n_homoplasic_sites'] / gene['n_sites']
    # Filter for loci in the top 20% of number of homoplasies
    min_hp_thresh_quantile = 0.5
    min_hp_thresh = gene.loc[gene.ranker > 0, 'n_homoplasic_sites'].quantile(min_hp_thresh_quantile)
    print('min_hp_thresh: >= {}, equiv. quantile: {}'.format(min_hp_thresh, min_hp_thresh_quantile))
    gene = gene.loc[(gene['ranker'] > 0) & (gene['n_homoplasic_sites'] >= min_hp_thresh)].copy()
    gene['rank'] = gene['ranker'].rank(method='average')
    # Use normalised rank as number of ranked genes differ between strains
    gene['rank_norm'] = gene['rank'] / len(gene['rank'])
    genes_ranked[short_prefix] = gene

with open('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_stats.txt', 'w') as stats_file:
    for short_prefix, gene in genes.items():
        # Get some simple summary stats
        print('~~~ {} ~~~'.format(short_prefix), file=stats_file)
        print(gene.iloc[:, [7, 10, 11, 12, 13, 14]].apply(sum), file=stats_file)
        print('agree_prop_acctran_deltran mean: {}'.format(gene.loc[gene.n_homoplasic_sites.astype(bool), 'agree_prop_acctran_deltran'].mean()), file=stats_file)
        print('n genes with homoplasies: {}'.format(len(gene.loc[(gene.n_homoplasic_sites.astype(bool)) & ~gene.intergenic])), file=stats_file)
        print('n genes: {}'.format(len(~gene.intergenic)), file=stats_file)
        print('n near-gene-regions with homoplasies: {}'.format(len(gene.loc[(gene.n_homoplasic_sites.astype(bool)) & gene.intergenic])), file=stats_file)
        print('n near-gene-regions: {}'.format(len(gene.intergenic)), file=stats_file)

# Read in reciprocal best hits tables and build adjaceny graph
print('Constructing homology graph...')
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
#  ranking by (agree_prop_acctran_deltran * n_homoplasic_sites)/n_total
#
print('Aggregating ranks...')
ranked_merged = pd.concat(genes_ranked.values(), keys=genes_ranked.keys()).reset_index()
ranked_merged.columns = ['st' if x == 'level_0' else x for x in ranked_merged.columns]
ranked_merged['locus_tag_simplified'] = ranked_merged['locus_tag'].map(lambda x: x.split("'")[1])
# Process genic and intergenic regions separately
ranked_merged_intergenic = ranked_merged[ranked_merged['intergenic']]
ranked_merged_genic = ranked_merged[~ranked_merged['intergenic']]

# TODO OOP :<
homolog_genic, homolog_intergenic = [], []
homolog_ranks_genic, homolog_ranks_intergenic = [], []
homolog_ranks_all_genic, homolog_ranks_all_intergenic = [], []
for c in nx.connected_components(G):
    homolog_genic.append(tuple(sorted(c)))
    homolog_intergenic.append(tuple(sorted(c)))
    #
    # Get mean rank of all locus_tags that are in the homologs connected component
    homolog_ranks_genic.append(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'rank_norm'].mean())
    homolog_ranks_intergenic.append(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'rank_norm'].mean())
    #
    homolog_ranks_all_genic.append(
            tuple(zip(
                tuple(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'st']), 
                tuple(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'locus_tag_simplified']), 
                tuple(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'rank_norm']),
                tuple(ranked_merged_genic.loc[ranked_merged_genic['locus_tag_simplified'].isin(c), 'n_homoplasic_sites'])
            ))
    )
    homolog_ranks_all_intergenic.append(
            tuple(zip(
                tuple(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'st']), 
                tuple(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'locus_tag_simplified']), 
                tuple(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'rank_norm']),
                tuple(ranked_merged_intergenic.loc[ranked_merged_intergenic['locus_tag_simplified'].isin(c), 'n_homoplasic_sites'])
            ))
    )

# Get gene and product annotations for locus tags
print('Getting annotations for loci...')
locus_annotations = {}
for short_prefix, embl_file in embl_files.items():
    for feature in next(SeqIO.parse(embl_file, 'embl')).features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers['locus_tag'][0]
            gene = feature.qualifiers['gene'] if 'gene' in feature.qualifiers else None
            product = feature.qualifiers['product'] if 'product' in feature.qualifiers else None 
            locus_annotations[locus_tag] = {'st': short_prefix, 'gene': gene, 'product': product}

#
# Mark known loci from 2010 paper
#
print('Marking known homoplasic loci...')
#
known_hps_annotated = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/misc/homoplasies_st239_harris_et_al_2010_science_table1.annotated.csv', index_col=0)
known_hps_annotated['locus_tag_simplified'] = known_hps_annotated['locus_tag'].map(lambda x: eval(x)[0][0])
#
known_hps_annotated_genic = known_hps_annotated[known_hps_annotated['intergenic'].map(lambda x: all(d == 0 for d in eval(x)))]
known_hps_annotated_intergenic = known_hps_annotated[~known_hps_annotated['intergenic'].map(lambda x: all(d == 0 for d in eval(x)))]
print('Read in {} known homoplasic sites, in {} unique loci, intergenic near {} unique loci.'.format(
    len(known_hps_annotated), 
    len(known_hps_annotated_genic['locus_tag_simplified'].unique()),
    len(known_hps_annotated_intergenic['locus_tag_simplified'].unique())
))

#
print('Writing summary...')
genic_summary = pd.DataFrame({
    'sts': [tuple(locus_annotations[locus_tag]['st'] for locus_tag in homologs) for homologs in homolog_genic], 
    'homologs': homolog_genic, 
    'genes': [tuple(locus_annotations[locus_tag]['gene'] for locus_tag in homologs) for homologs in homolog_genic], 
    'products': [tuple(locus_annotations[locus_tag]['product'] for locus_tag in homologs) for homologs in homolog_genic], 
    'rank_norms': homolog_ranks_all_genic, 
    'n_sts': [len(x) for x in homolog_ranks_all_genic],
    'mean_rank_norm': homolog_ranks_genic, 
    'known_hp_harris_et_al_2010': [any(locus in set(known_hps_annotated_genic['locus_tag_simplified']) for locus in homologs) for homologs in homolog_genic]
}).sort_values('mean_rank_norm', ascending=False)
#
intergenic_summary = pd.DataFrame({
    'sts': [tuple(locus_annotations[locus_tag]['st'] for locus_tag in homologs) for homologs in homolog_intergenic], 
    'homologs': homolog_intergenic, 
    'genes': [tuple(locus_annotations[locus_tag]['gene'] for locus_tag in homologs) for homologs in homolog_intergenic], 
    'products': [tuple(locus_annotations[locus_tag]['product'] for locus_tag in homologs) for homologs in homolog_intergenic], 
    'rank_norms': homolog_ranks_all_intergenic, 
    'n_sts': [len(x) for x in homolog_ranks_all_intergenic],
    'mean_rank_norm': homolog_ranks_intergenic, 
    'known_hp_harris_et_al_2010': [any(locus in set(known_hps_annotated_intergenic['locus_tag_simplified']) for locus in homologs) for homologs in homolog_intergenic]
}).sort_values('mean_rank_norm', ascending=False)
#
print('Marked known {} ranked homolog groups (genic), {} intergenic.'.format(
    genic_summary[~pd.isnull(genic_summary.mean_rank_norm)]['known_hp_harris_et_al_2010'].sum(), 
    intergenic_summary[~pd.isnull(intergenic_summary.mean_rank_norm)]['known_hp_harris_et_al_2010'].sum()
))

# Missing value of mean_rank_norm means none of the locus tags in the homology group were present after filtering
genic_summary.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_genic.csv')
intergenic_summary.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_intergenic.csv')

