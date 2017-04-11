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
    #
    df = df.copy()
    if not keep_all_sites:
        # Filter for SNPs that occur at the same position
        df = df[df['loc'].duplicated(keep=False)]
    #
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
    #
    # Summarise info at unique positions
    # 
    # Count number of changes with the following properties:
    # convergence = SNP final base is the same as another SNP's final base
    # reversal = SNP final base reverts to ancestral base of another change e.g. loc=56478, st22
    # reversed = An ancestral change that is reversed in a descendent. Opposite to reversal.
    # other = non convergence, non reversal/reversed changes
    hp = pd.DataFrame()
    hp = hp.assign(
        change=df.groupby('loc').apply(lambda x: collections.Counter(x['change'])),
        n_convergence=df.groupby('loc').apply(lambda x: sum('convergence' in snp for snp in x['homoplasy'])),
        n_reversal=df.groupby('loc').apply(lambda x: sum('reversal' in snp for snp in x['homoplasy'])),
        n_both=df.groupby('loc').apply(lambda x: sum('convergence' in snp and 'reversal' in snp for snp in x['homoplasy'])),
        n_homoplasic=df.groupby('loc').apply(lambda x: sum('convergence' in snp or 'reversal' in snp for snp in x['homoplasy'])),
        n_non_homoplasic=df.groupby('loc').apply(lambda x: sum(not 'convergence' in snp and not 'reversal' in snp for snp in x['homoplasy'])),
        n_total=df.groupby('loc').apply(lambda x: len(x['homoplasy'])),
        branches=df.groupby('loc').apply(lambda x: list(x['homoplasy']))
    )
    #
    return(hp)

def build_intervaltree(features):
    '''Interval tree for fast search for containing intervals
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
    # No need check SNP strand, as convergence/reversal is tagged regardless of strand of SNP.
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
    # [] Indicates lack of feature containing a site
    # None indicates lack of a particular qualifier value describing a contained site
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

    transformations = ['acctran', 'deltran']
    csv_files = {
        "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_homoplasies_on_tree.{}.csv",
        "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_homoplasies_on_tree.{}.csv",
        "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_homoplasies_on_tree.{}.csv",
        "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_homoplasies_on_tree.{}.csv"
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

    dfs = {}
    for t in transformations:
        print("Reading {} .tab files...".format(t))
        dfs[t] = pd.concat([pd.read_csv(f.format(t)) for f in csv_files.values()], keys=csv_files.keys(), names=["st", "i"])

    for prefix in csv_files.keys():

        #  prefix = 'st239'

        hps = {}
        for t in transformations:
            print(t)
            df = dfs[t].loc[prefix]
            print('{}: {} features in .tab file.'.format(prefix, len(df)))
            hp = get_homoplasies(df, keep_all_sites=True)
            print('{}: {} sites in .tab file.'.format(prefix, len(hp)))
            print('{}: {} homoplasic sites in .tab file.'.format(prefix, len(hp.query('n_total > 1'))))
            # Check the counts are partitioned right
            assert(all(hp[hp['n_convergence'] + hp['n_reversal'] - hp['n_both'] + hp['n_non_homoplasic'] == hp['n_total']]))
            assert(all(hp[hp['n_homoplasic'] + hp['n_non_homoplasic'] == hp['n_total']]))
            #
            hps[t] = hp

        # Compare acctran and deltran agreement
        print('Comparing acctran and deltran...')
        hp_merged = hps['acctran'].merge(hps['deltran'], left_index=True, right_index=True, suffixes=['_acctran', '_deltran'])
        # Assert that total number of changes at each site is the same
        assert(all(hp_merged['n_total_acctran'] == hp_merged['n_total_deltran']))
        #  Assert that homoplasy assignments are identical between transformations
        #
        #  NOTE: this is not true if reversed is not counted e.g. st239, loc=1524413
        #  acctran
        #  ['reversed in branch leading to 6949_5#5, reversed in branch leading to 6949_5#14',
        #  'reversal from branch leading to 365, convergence with branch leading to 6949_5#14',
        #  'reversal from branch leading to 365, convergence with branch leading to 6949_5#5']
        #  deltran
        #  ['reversed in branch leading to 6949_5#5, convergence with branch leading to 366',
        #  'convergence with branch leading to 362',
        #  'reversal from branch leading to 362'],
        #
        #  assert(all(hp_merged['n_non_homoplasic_acctran'] == hp_merged['n_non_homoplasic_deltran']))
        
        hp_merged['agree_acctran_deltran'] = ((hp_merged['n_convergence_acctran'] == hp_merged['n_convergence_deltran']) &
                                        (hp_merged['n_reversal_acctran'] == hp_merged['n_reversal_deltran']) &
                                        (hp_merged['n_both_acctran'] == hp_merged['n_both_deltran']))

        # Annotate sites
        hp_merged = annotate_homoplasies(hp_merged, embl_files[prefix])
        print('{}: {} annotations added.'.format(prefix, len(hp_merged[~hp_merged['gene'].apply(len).astype(bool)])))

        #  Write out homoplasies
        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasies.csv")
        hp_merged[hp_merged['gene'].astype(bool)].sort_values(['n_convergence_acctran', 'n_convergence_deltran'], ascending=False).to_csv(outfile)

        def nested_lists_to_tuples(ll):
            '''Convert nested lists to nested tuples to allow hashing
            '''
            if type(ll) == list:
                return(tuple(nested_lists_to_tuples(l) for l in ll))
            else:
                return(ll)

        def get_ordered_unique_values(x):
            '''Return unique values in an iterable of hashable elements, retaining order
            '''
            return(tuple(collections.OrderedDict(itertools.zip_longest(x, (None, ))).keys()))
        
        #  Get genes with the highest ratio of convergences to total changes
        print('Summarising by gene...')
        summary_series = hp_merged['locus_tag'].apply(nested_lists_to_tuples)
        genes = pd.DataFrame()
        genes['gene'] = hp_merged.groupby(summary_series)['gene'].apply(lambda x: get_ordered_unique_values(nested_lists_to_tuples(list(x))))
        genes['product'] = hp_merged.groupby(summary_series)['product'].apply(lambda x: get_ordered_unique_values(nested_lists_to_tuples(list(x))))
        genes['pseudo'] = hp_merged.groupby(summary_series)['pseudo'].apply(lambda x: get_ordered_unique_values(nested_lists_to_tuples(list(x))))
        genes['strand'] = hp_merged.groupby(summary_series)['strand'].apply(lambda x: get_ordered_unique_values(nested_lists_to_tuples(list(x))))
        #
        genes['n_sites'] = hp_merged.groupby(summary_series)['strand'].apply(len)
        genes['agree_acctran_deltran'] = hp_merged.groupby(summary_series).apply(lambda x: collections.Counter(x['agree_acctran_deltran']))
        genes['agree_prop_acctran_deltran'] = genes['agree_acctran_deltran'].apply(lambda x: x[True]/(x[True]+x[False]))
        #
        genes['n_convergence_acctran'] = hp_merged.groupby(summary_series)['n_convergence_acctran'].apply(sum)
        genes['n_convergence_deltran'] = hp_merged.groupby(summary_series)['n_convergence_deltran'].apply(sum)
        genes['n_reversal_acctran'] = hp_merged.groupby(summary_series)['n_reversal_acctran'].apply(sum)
        genes['n_reversal_deltran'] = hp_merged.groupby(summary_series)['n_reversal_deltran'].apply(sum)
        genes['n_both_acctran'] = hp_merged.groupby(summary_series)['n_both_acctran'].apply(sum)
        genes['n_both_deltran'] = hp_merged.groupby(summary_series)['n_both_deltran'].apply(sum)
        genes['n_homoplasic_acctran'] = hp_merged.groupby(summary_series)['n_homoplasic_acctran'].apply(sum)
        genes['n_homoplasic_deltran'] = hp_merged.groupby(summary_series)['n_homoplasic_deltran'].apply(sum)
        genes['n_total'] = hp_merged.groupby(summary_series)['n_total_acctran'].apply(sum)
        #
        genes['n_h_acctran/n_t'] = genes['n_homoplasic_acctran']/genes['n_total']
        genes['n_h_deltran/n_t'] = genes['n_homoplasic_deltran']/genes['n_total']

        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasies_per_gene.csv")
        genes.sort_values(['n_h_acctran/n_t', 'n_h_deltran/n_t'], ascending=False).to_csv(outfile)

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

        if prefix == 'st22' and False:

            #  An st22 convergence and other
            df[df['loc'] == 2103759] 

            #  An st22 reversal
            df[df['loc'] == 2592254] 

            #  Convergence and reversal
            df.query('loc == 56478')

            #  Multiple changes with no homoplasy
            df.query('loc == 2809075') 

            #  Sites where one snp is on the reverse, one is on the forward,
            #  but the change is a convergence/reversal, are still marked as convergence/reversal.
            #  It's just that the tag Intergenic is applied to one of them
            df.groupby('loc').apply(lambda x: len(x['strand'].unique())>1).sort_values().tail(15)
            #  1864667     True
            #  1430573     True
            #  1430575     True
            #  1431929     True
            #  2737583     True
            #  1432868     True
            #  1430576     True
            #  1432301     True
            #  1432697     True
            #  2737568     True
            #  2737547     True
            #  2737559     True

