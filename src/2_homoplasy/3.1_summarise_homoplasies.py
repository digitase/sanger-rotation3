#!/usr/bin/env python
'''Summarise homoplasic SNPs
'''

import os
import pandas as pd
import numpy as np
import collections
import itertools
import intervaltree
from Bio import SeqIO

def get_homoplasies(df, keep_all_sites=False):
    '''Find homoplasic sites
    '''
    df = df.copy()
    if not keep_all_sites:
        # Filter for SNPs that occur at the same position
        df = df[df['loc'].duplicated(keep=False)]
    df = df.sort_values(['loc', 'strand', 'SNP'])
    # Refactor some values
    df['strand'].replace(np.nan, '"forward"', inplace=True)
    df['homoplasy'].replace(np.nan, '', inplace=True)
    df['codon_type'].replace({'SNOP (non-stop codon to stop condon)': 'STOPgain', 'STIP (stop codon to non-stop codon)': 'STOPloss'}, inplace=True)
    # Strip quotes from strings
    for col in ['codon_type', 'homoplasy', 'taxa', 'strand', 'node', 'SNP']:
        df[col] = df[col].apply(lambda x: x.replace('"', ''))
    # Combine SNP, change type, and strand into one
    df = df.assign(change=df[['SNP', 'codon_type', 'strand']].apply(tuple, axis=1))
    # Summarise info at unique positions
    # 
    # Count number of changes with the following properties:
    # convergence = SNP final base is the same as another SNP's final base
    # reversal = SNP final base reverts to ancestral base of another change e.g. loc=56478, st22
    # reversed = Opposite change to reversal
    # other = non convergence, non reversal/reversed changes
    hp = pd.DataFrame()
    hp = hp.assign(
        change=df.groupby('loc').apply(lambda x: collections.Counter(x['change'])),
        n_convergence=df.groupby('loc').apply(lambda x: sum('convergence' in snp for snp in x['homoplasy'])),
        n_reversal=df.groupby('loc').apply(lambda x: sum('reversal' in snp for snp in x['homoplasy'])),
        n_reversed=df.groupby('loc').apply(lambda x: sum('reversed' in snp for snp in x['homoplasy'])),
        n_other=df.groupby('loc').apply(lambda x: sum(not snp for snp in x['homoplasy'])),
        n_total=df.groupby('loc').apply(lambda x: len(x['homoplasy'])),
        #  node=df.groupby('loc').apply(lambda x: list(x['node']))
    )
    hp = hp.assign(n_nonother=hp['n_total']-hp['n_other'])
    #
    return(hp)

def build_intervaltree(features):
    '''
    '''
    tree = intervaltree.IntervalTree()
    for f in features:
        # Add each part of a CompoundLocation
        for part in f.location.parts:
            tree.addi(int(part.start), int(part.end), f)
    return(tree)


def annotate_homoplasies(hp, embl_file):
    '''Annotate homoplasies with gene info
    '''
    print('Parsing EMBL file... {}'.format(embl_file))
    annot_record = next(SeqIO.parse(embl_file, 'embl'))
    #
    # Build dictionary of SNP location -> gene
    #
    # Note: loc is 1-indexed, biopython converts embl coords to 0 indexed on parse
    # TODO check if we need to check SNP strand
    print('Building SNP -> gene dict...')
    loc_features = collections.defaultdict(list)
    #
    # Brute force O(n^2) solution
    #  for feature in annot_record.features:
        #  if feature.type == 'CDS':
            #  for loc in hp.index:
                #  if int(loc-1) in feature:
                    #  loc_features[loc].append(feature)
    # intervaltree solution
    tree = build_intervaltree(annot_record.features)
    for loc in hp.index:
        # Find tree nodes with features containing loc
        intervals = tree.search(int(loc-1))
        intervals = filter(lambda x: x.data.type == 'CDS', intervals)
        # Add unique features to dict
        loc_features[loc] = list(set(
            i.data for i in intervals
        ))
    #
    # Append annotations to df
    #
    # None indicates lack of qualifier value
    # [] Indicates lack of feature
    print('Adding annotations...')
    loc_genes = [loc_features[loc] for loc in hp.index]
    # These are nested lists, as qualifier values are lists e.g. 'gene': ['dnaA', 'dnaH']
    hp['gene'] = [list(map(lambda x: x.qualifiers['gene'] if 'gene' in x.qualifiers else None, loc_gene)) for loc_gene in loc_genes]
    hp['product'] = [list(map(lambda x: x.qualifiers['product'] if 'product' in x.qualifiers else None, loc_gene)) for loc_gene in loc_genes]
    hp['locus_tag'] = [list(map(lambda x: x.qualifiers['locus_tag'] if 'locus_tag' in x.qualifiers else None, loc_gene)) for loc_gene in loc_genes]
    hp['pseudo'] = [list(map(lambda x: '1' if 'pseudo' in x.qualifiers else 0, loc_gene)) for loc_gene in loc_genes]
    # These are non-nested lists, as a feature has a single strand/type
    hp['strand'] = [list(map(lambda x: x.strand, loc_gene)) for loc_gene in loc_genes]
    hp['type'] = [list(map(lambda x: x.type, loc_gene)) for loc_gene in loc_genes]
    #
    return(hp)

if __name__ == "__main__":
    csv_files = {
        "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_homoplasies_on_tree.csv",
        "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_homoplasies_on_tree.csv",
        "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_homoplasies_on_tree.csv",
        "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_homoplasies_on_tree.csv"
    }
    embl_files = {
        "st22": '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/CC22_EMRSA15.embl',
        "st239": '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st239/CC8_TW20.embl',
        "st30": '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/CC30_MRSA252.embl',
        "st8": '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/CC8_USA300_FPR3757.embl'
    }
    out_file_prefixes = {
        "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/",
        "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/",
        "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/",
        "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/"
    }
    dfs = pd.concat([pd.read_csv(f) for f in csv_files.values()], keys=csv_files.keys(), names=["st", "i"])

    for prefix in csv_files.keys():

        #  prefix = 'st30'

        df = dfs.loc[prefix]
        print('{}: {} features in .tab file.'.format(prefix, len(df)))
        hp = get_homoplasies(df, keep_all_sites=True)
        print('{}: {} sites in .tab file.'.format(prefix, len(hp)))
        print('{}: {} homoplasic sites in .tab file.'.format(prefix, len(hp.query('n_total > 1'))))
        hp = annotate_homoplasies(hp, embl_files[prefix])

        #  Write out homoplasies
        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasies.csv")
        hp[hp['gene'].astype(bool)].sort_values('n_convergence', ascending=False).to_csv(outfile)

        #  Get genes with the highest ratio of convergences to total changes
        genes = pd.DataFrame()
        genes['product'] = hp.groupby(hp['gene'].apply(str))['product'].apply(lambda x: list(x)[0])
        genes['n_total'] = hp.groupby(hp['gene'].apply(str))['n_total'].apply(sum)
        genes['n_nonother'] = hp.groupby(hp['gene'].apply(str))['n_nonother'].apply(sum)
        genes['n_convergence'] = hp.groupby(hp['gene'].apply(str))['n_convergence'].apply(sum)
        genes['n_c/n_t'] = genes['n_convergence']/genes['n_total']

        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasic_genes.csv")
        genes.sort_values('n_c/n_t', ascending=False).to_csv(outfile)

        #  It is not the case that all sites with a gene annotation are not tagged Intergenic,
        #  as the 'Intergenic' tag is given if the SNP strand does not match the gene strand
        #  e.g. st22:
        #  hp.ix[98558]
        #  hp.ix[1430575]
        #  hp.ix[1864667]
        #  hp.ix[2776899]

        #  Nor is it the case that all sites without gene annotation are tagged Intergenic,
        #  as STIP/SNOP mutations that represent a break in the EMBL annotation of a gene
        #  will not be counted as in a gene.
        #  e.g. st30:
        #  hp.ix[2474464]
        #  Nonsense mutation in .embl: complement(join(2474081..2474461,2474465..2474998))

        #  if prefix == 'st22':
        #    df[df['loc'] == 2103759] #  An st22 convergence and other
        #    df[df['loc'] == 2592254] #  An st22 reversal
        #    df.query('loc == 56478') #  Convergence and reversal
        #    df.query('loc == 2809075') #  Multiple changes with no homoplasy

