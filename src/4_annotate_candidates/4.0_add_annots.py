
'''Parse annotation results and add annotations to unannotated candidates'''

import pandas as pd
import collections

summary_sites_by_gene_filtered = pd.read_csv('.output/hp_sites_summary.filtered.csv', converters={'genes': eval, 'locus_tags': eval, 'products': eval})
summary_genic_filtered = pd.read_csv('.output/hp_homologs_summary_genic.filtered.csv', converters={'gene': eval, 'homologs': eval, 'product': eval})
summary_intergenic_filtered = pd.read_csv('.output/hp_homologs_summary_intergenic.filtered.csv', converters={'gene': eval, 'homologs': eval, 'product': eval})

def flatten_list_or_tuple(xs):
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from flatten_list_or_tuple(x)
    else:
        yield xs

# Grouping of sites into loci by adjacent identical locus tags
tags_list = summary_sites_by_gene_filtered['locus_tags'].map(lambda x: set(flatten_list_or_tuple(x)))
group_n = 0
grouping = []
for i, tags in enumerate(tags_list):
    if i == 0:
        pass
    else:
        if tags_list[i] & tags_list[i-1]:
            pass
        else:
            group_n += 1
    grouping.append(group_n)
#
summary_sites_by_gene_filtered['locus_group'] = grouping

# Read in iprscan results
iprscan = pd.read_table('.output/interproscan/interproscan_out.tsv', header=None)
iprscan.columns = ['prot_acc', 'seq_md5', 'seq_len', 'analysis', 'sig_acc', 'sig_desc', 'start_loc', 'stop_loc', 'score', 'status', 'date', 'ipr_annot_acc', 'ipr_annot_desc', 'go_annot', 'path_annot']

# Read in blastp results
blastp = pd.read_table('.output/blastp/st_all.cds.pep.unannotated.blastp.outfmt6cust', header=None)
blastp.columns = "qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()

locus_to_annot_iprscan = collections.defaultdict(set)
for row in iprscan.itertuples():
    locus_to_annot_iprscan[row.prot_acc].add((row.analysis, row.sig_acc, row.sig_desc))
    locus_to_annot_iprscan[row.prot_acc].add(('GO_KEGG', row.go_annot, row.path_annot))

locus_to_annot_blastp = collections.defaultdict(set)
for row in blastp.itertuples():
    locus_to_annot_blastp[row.qseqid].add(('blastp', row.sseqid, row.stitle))

def flatten_list_or_tuple(xs):
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from flatten_list_or_tuple(x)
    else:
        yield xs

def locus_tags_to_annot_str(xs, locus_to_annot):
    '''Get annotations for locus tags
    '''
    def locus_tag_to_annot_str(locus_tag, locus_to_annot):
        annot = locus_to_annot[locus_tag]
        if annot:
            return '; '.join('{}={},{}'.format(*a) for a in annot if a)
    if isinstance(xs, (list, tuple)):
        for x in xs:
            yield from locus_tags_to_annot_str(x, locus_to_annot)
    else:
        yield locus_tag_to_annot_str(xs, locus_to_annot)

summary_sites_by_gene_filtered['annot_iprscan'] = summary_sites_by_gene_filtered['locus_tags'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_iprscan)))
summary_genic_filtered['annot_iprscan'] = summary_genic_filtered['homologs'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_iprscan)))
summary_intergenic_filtered['annot_iprscan'] = summary_intergenic_filtered['homologs'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_iprscan)))

summary_sites_by_gene_filtered['annot_blastp'] = summary_sites_by_gene_filtered['locus_tags'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_blastp)))
summary_genic_filtered['annot_blastp'] = summary_genic_filtered['homologs'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_blastp)))
summary_intergenic_filtered['annot_blast'] = summary_intergenic_filtered['homologs'].map(lambda ls: tuple(locus_tags_to_annot_str(ls, locus_to_annot_blastp)))

summary_sites_by_gene_filtered.to_csv('.output/hp_sites_summary.filtered.annotated.csv')
summary_genic_filtered.to_csv('.output/hp_homologs_summary_genic.filtered.annotated.csv')
summary_intergenic_filtered.to_csv('.output/hp_homologs_summary_intergenic.filtered.annotated.csv')

