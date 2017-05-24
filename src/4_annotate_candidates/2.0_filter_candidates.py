'''Filter list of candidates e.g. remove known AMR loci
'''

import pandas as pd
#  from fuzzywuzzy import fuzz
#  from fuzzywuzzy import process
from collections import Counter

def flatten_list_or_tuple(xs):
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from flatten_list_or_tuple(x)
    else:
        yield xs

def tag_known_loci(candidates, known_amr_loci):
    for gs in candidates:
        gs = flatten_list_or_tuple(gs)
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

    summary_sites = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_sites_summary.csv', converters={'gene': eval, 'locus_tag': eval, 'product': eval})
    summary_sites['known_amr_loci'] = list(tag_known_loci(summary_sites['gene'], known_amr_loci))

    def get_synonymous_from_changes(counter):
        '''Converts Counter({('G->A', 'Synonymous', 'reverse'): 2, ('A->G', 'Synonymous', 'reverse'): 1}) to Counter({'Synonymous': 3})
        '''
        return Counter(x[1] for x in counter.elements())

    summary_sites_by_gene = pd.DataFrame().assign(
        sts=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['level_0'])),
        loc_in_st=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['loc'])),
        n_sts=summary_sites.groupby('loc_in_alignment').apply(len),
        synonymous_acctran=summary_sites.groupby('loc_in_alignment').apply(lambda x: sum((get_synonymous_from_changes(eval(c)) for c in x['change_acctran']), Counter())),
        synonymous_deltran=summary_sites.groupby('loc_in_alignment').apply(lambda x: sum((get_synonymous_from_changes(eval(c)) for c in x['change_deltran']), Counter())),
        n_homoplasic_acctran=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_homoplasic_acctran'])),
        n_homoplasic_deltran=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_homoplasic_deltran'])),
        n_total=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['n_total_acctran'])),
        locus_tags=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['locus_tag'])),
        products=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['product'])),
        genes=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['gene'])),
        intergenic=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['intergenic'])),
        known_amr_loci=summary_sites.groupby('loc_in_alignment').apply(lambda x: tuple(x['known_amr_loci']))
    )

    summary_genic = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_genic.csv', converters={'genes': eval, 'homologs': eval, 'products': eval})
    summary_genic = summary_genic.query('n_sts > 1')
    summary_genic['known_amr_loci'] = list(tag_known_loci(summary_genic['genes'], known_amr_loci))

    summary_intergenic = pd.read_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_homologs_summary_intergenic.csv', converters={'genes': eval, 'homologs': eval, 'products': eval})
    summary_intergenic = summary_intergenic.query('n_sts > 1')
    summary_intergenic['known_amr_loci'] = list(tag_known_loci(summary_intergenic['genes'], known_amr_loci))

    # Remove records with known amr loci
    summary_sites_by_gene_filtered = summary_sites_by_gene[summary_sites_by_gene['known_amr_loci'].map(lambda sts: all(len(x) == 0 for x in sts))]
    summary_genic_filtered = summary_genic[~summary_genic['known_amr_loci'].astype(bool)]
    summary_intergenic_filtered = summary_intergenic[~summary_intergenic['known_amr_loci'].astype(bool)]

    summary_sites_by_gene_filtered.to_csv('.output/hp_sites_summary.filtered.csv')
    summary_genic_filtered.to_csv('.output/hp_homologs_summary_genic.filtered.csv')
    summary_intergenic_filtered.to_csv('.output/hp_homologs_summary_intergenic.filtered.csv')

    # Output non-annotated locus tags i.e. no gene name OR putative or hypothetical in products
    def add_unannotated_locus_tags(gss, pss, lss, existing_results=None):
        '''Add locus tags without annotation to a set
        '''
        if not existing_results:
            existing_results = set()
        for gs, ps, ls in zip(gss, pss, lss):
            if all(g == None for g in gs):
                existing_results.update(ls)
            for p, l in zip(ps, ls):
                if 'hypothetical' in str(p) or 'putative' in str(p):
                    existing_results.add(l)
        return(existing_results)
    #
    locus_tags_unannotated = add_unannotated_locus_tags(
        summary_sites_by_gene_filtered['genes'].map(lambda x: list(flatten_list_or_tuple(x))),
        summary_sites_by_gene_filtered['products'].map(lambda x: list(flatten_list_or_tuple(x))),
        summary_sites_by_gene_filtered['locus_tags'].map(lambda x: list(flatten_list_or_tuple(x)))
    )
    locus_tags_unannotated = add_unannotated_locus_tags(
        summary_genic_filtered['genes'].map(lambda x: list(flatten_list_or_tuple((x)))),
        summary_genic_filtered['products'].map(lambda x: list(flatten_list_or_tuple((x)))),
        summary_genic_filtered['homologs'].map(lambda x: list(flatten_list_or_tuple((x)))),
        locus_tags_unannotated
    )
    locus_tags_unannotated = add_unannotated_locus_tags(
        summary_intergenic_filtered['genes'].map(lambda x: list(flatten_list_or_tuple((x)))),
        summary_intergenic_filtered['products'].map(lambda x: list(flatten_list_or_tuple((x)))),
        summary_intergenic_filtered['homologs'].map(lambda x: list(flatten_list_or_tuple((x)))),
        locus_tags_unannotated
    )

    with open('.output/locus_tags_unannotated.txt', 'w') as outfile:
        outfile.write('\n'.join(locus_tags_unannotated))

