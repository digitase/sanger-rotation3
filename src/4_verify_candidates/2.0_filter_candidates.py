'''Filter list of candidates e.g. remove known AMR loci
'''

import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

def flatten_list_or_tuple(xs):
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from flatten_list_or_tuple(x)
    else:
        yield xs

def tag_known_loci(candidates, known_amr_loci):
    for gs in candidates:
        gs = flatten_list_or_tuple(eval(gs))
        gs_result = set()
        for g in gs:
            if g:
                for known in known_amr_loci:
                    if known.startswith(g.lower()[:3]):
                        gs_result.add(known)
                # print(process.extract(g.lower(), known_amr_loci, limit=2))
        yield gs_result

if __name__ == "__main__":

    with open('.output/known_amr_loci.txt', 'r') as infile:
        known_amr_loci = [x.strip().lower() for x in infile.readlines()]

    summary_sites = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_sites_summary.csv')
    summary_sites['known_amr_loci'] = list(tag_known_loci(summary_sites['gene'], known_amr_loci))
    summary_sites_by_gene = pd.DataFrame().assign(
        sts=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['level_0'])),
        loc_in_st=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['loc'])),
        n_sts=summary_sites.groupby('loc_in_alignment').apply(len),
        n_homoplasic_acctran=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_homoplasic_acctran'])),
        n_homoplasic_deltran=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_homoplasic_deltran'])),
        n_total=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_total_acctran'])),
        locus_tag=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['locus_tag'])),
        product=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['product'])),
        gene=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['gene'])),
        intergenic=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['intergenic'])),
        known_amr_loci=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['known_amr_loci']))
    )

    summary_genic = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_genic.csv')
    summary_genic = summary_genic.query('n_sts > 1')
    summary_genic['known_amr_loci'] = list(tag_known_loci(summary_genic['genes'], known_amr_loci))

    summary_intergenic = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_intergenic.csv')
    summary_intergenic = summary_intergenic.query('n_sts > 1')
    summary_intergenic['known_amr_loci'] = list(tag_known_loci(summary_intergenic['genes'], known_amr_loci))

    # Remove records with known amr loci
    summary_sites_by_gene_filtered = summary_sites_by_gene[summary_sites_by_gene['known_amr_loci'].map(lambda sts: all(len(x) == 0 for x in sts))]
    summary_genic_filtered = summary_genic[~summary_genic['known_amr_loci'].astype(bool)]
    summary_intergenic_filtered = summary_intergenic[~summary_intergenic['known_amr_loci'].astype(bool)]

    summary_sites_by_gene_filtered.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_verify_candidates/hp_sites_summary.filtered.csv')
    summary_genic_filtered.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_verify_candidates/hp_homologs_summary_genic.filtered.csv')
    summary_intergenic_filtered.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_verify_candidates/hp_homologs_summary_intergenic.filtered.csv')

