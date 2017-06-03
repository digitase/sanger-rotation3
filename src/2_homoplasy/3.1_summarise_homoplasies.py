#!/usr/bin/env python
'''Summarise homoplasic SNPs
'''

import os
import pandas as pd
import numpy as np
import collections
import itertools
#  import intervaltree
from Bio import SeqIO
import pybedtools

def get_homoplasies(df, keep_all_sites=False):
    '''Find homoplasic sites

    df: Dataframe of every change on the tree
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
    # convergence = SNP final base is the same as another SNP's final base. Use n-1, as in a pair of changes, 1 change alone is not a convergence
    # reversal = SNP final base reverts to ancestral base of another change e.g. loc=56478, st22
    # reversed = An ancestral change that is reversed in a descendent. Opposite to reversal.
    # other = non convergence, non reversal/reversed changes

    #  def interpret_hps(loc_df):
        #  '''Given a df of changes at a location, parse and tally homoplasies
        #  '''
        #  n_convergence, n_reversal, n_both, n_non_homoplasic, n_total = (0, 0, 0, 0, 0)
        # Group by final base (convergence groups)
        #  loc_df['SNP_after'] = loc_df['SNP'].map(lambda x: x.split('->')[1])
        #  for snp, snp_df in loc_df.groupby('SNP_after'):
            #  n_total_snp = len(snp_df)
            #  n_convergence_snp = len([x for x in snp_df['homoplasy'] if 'convergence' in x])
            # Not possible to have a convergence with only a single change
            #  assert(n_convergence_snp == 0 or n_convergence_snp >= 2)
            # n-1 to correct for first change in a convergence group
            #  if n_convergence_snp: n_convergence_snp -= 1
            # no need for n-1 as 'reversed is not counted'
            #  n_reversal_snp = len([x for x in snp_df['homoplasy'] if 'reversal' in x])
            # n_both_snp is ambiguous, as it depends on which convergence is not counted in the n-1
            # hence n_non_homoplasic is also ambiguous
            #
            #  n_convergence += n_convergence_snp
            #  n_reversal += n_reversal_snp
            #  n_total += n_total_snp
        #  return({'n_convergence': n_convergence, 'n_reversal': n_reversal, 'n_total': n_total})
    
    # Group by final base (convergence groups)
    df['SNP_after'] = df['SNP'].map(lambda x: x.split('->')[1])

    print('Summarising homoplasies...')
    hp = pd.DataFrame()
    hp = hp.assign(
        change=df.groupby('loc').apply(lambda x: collections.Counter(x['change'])),
        # Subtract 1 for each convergence group, as the first change in a group is not a convergence
        n_convergence=df.groupby('loc').apply(lambda x: 
            max(
                0,
                sum('convergence' in snp for snp in x['homoplasy']) - len(x[x['SNP_after'].duplicated()]['SNP_after'].unique())
            )
        ),
        n_reversal=df.groupby('loc').apply(lambda x: sum('reversal' in snp for snp in x['homoplasy'])),
        #  n_both=df.groupby('loc').apply(lambda x: sum('convergence' in snp and 'reversal' in snp for snp in x['homoplasy'])),
        #  n_homoplasic=df.groupby('loc').apply(lambda x: sum('convergence' in snp or 'reversal' in snp for snp in x['homoplasy'])),
        #  n_non_homoplasic=df.groupby('loc').apply(lambda x: sum(not 'convergence' in snp and not 'reversal' in snp for snp in x['homoplasy'])),
        n_total=df.groupby('loc').apply(lambda x: len(x)),
    )
    #
    hp = hp.assign(
        n_homoplasic=hp['n_convergence'] + hp['n_reversal']
    )
    return(hp)

#  def build_intervaltree(features):
    #  '''Interval tree for fast search for containing intervals
    #  '''
    #  tree = intervaltree.IntervalTree()
    #  for f in features:
        # Add each part of a CompoundLocation
        #  for part in f.location.parts:
            #  tree.addi(int(part.start), int(part.end), f)
    #  return(tree)

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
    accepted_types = {'CDS', 'ncRNA', 'rRNA', 'tRNA'}
    #
    # pybedtools solution
    #
    # Create Interval for each feature part
    uid = 0
    # maintain mapping of feature uids back to feature parts
    uid_to_feature = dict()
    intervals = []
    for feature in [f for f in annot_record.features if f.type in accepted_types]:
        for part in feature.location.parts:
            # coords are already 0 based once read in by Bio
            interval = pybedtools.Interval(chrom=1, start=int(part.start), end=int(part.end), name=uid, strand=int(part.strand))
            intervals.append(interval)
            uid_to_feature[uid] = feature
            uid += 1
    #  Build BedTool from intervals and locs
    features = pybedtools.BedTool(intervals).sort()
    # locs are 1 indexed positions
    locs = pybedtools.BedTool(pybedtools.Interval(chrom=1, start=l-1, end=l, name="loc", strand=1) for l in hp.index)
    #
    # Find closest feature to each loc, ignore strand when searching for overlap
    # NOTE: does not account for circular chrom
    loc_to_features = collections.defaultdict(set)
    for result in locs.closest(features, d=True, t='all'):
        # This loc is 1 indexed
        _, _, loc, _, _, _, _, f_start, f_end, f_uid, _, f_strand, dist = result.fields
        loc, f_start, f_end, f_uid, f_strand, dist = map(int, (loc, f_start, f_end, f_uid, f_strand, dist))
        # Fix the sign of dist, because bedtools is STUPID
        # Upstream of a gene is negative
        if dist == 0:
            pass
        elif dist > 0:
            if f_strand == 1:
                # -1 to compare with 0 indexed coords
                if loc-1 < f_start:
                    sign = -1
                else:
                    assert(loc-1 >= f_end)
                    sign = 1
            elif f_strand == -1:
                if loc-1 < f_start:
                    sign = 1
                else:
                    assert(loc-1 >= f_end)
                    sign = -1
            else:
                raise(ValueError(f_strand))
            dist = dist * sign 
        else:
            raise(ValueError(dist))
        feature = uid_to_feature[f_uid]
        loc_to_features[loc].add((feature, dist))

        # test case with combos of plus and minus strand genes and left and right exact positions (st239)
        #  ([(x,y) for (x,y) in loc_to_features.items() if x in [655367, 2850827, 1179671, 2064502]])

    def get_qualifier_value_with_dist(feature_parts, qualifier):
        '''Summarise unique qualifier values within the given feature_parts, along with distance to closest feature
        '''
        dists = tuple(x[1] for x in feature_parts)
        genes = tuple(tuple(x[0].qualifiers[qualifier]) if qualifier in x[0].qualifiers else None for x in feature_parts)
        intergenic = all(dist != 0 for dist in dists)
        # Include distance to closest feature for intergenic locs
        if intergenic:
            result = tuple(set(zip(genes, dists)))
        else:
            result = tuple(set(genes))
        return(result)

    def get_qualifier_value(feature_parts, qualifier):
        '''Summarise unique qualifier values within the given feature_parts
        '''
        result = tuple(tuple(x[0].qualifiers[qualifier]) if qualifier in x[0].qualifiers else None for x in feature_parts)
        result = tuple(set(result))
        return(result)

    def get_dists_for_intergenic(feature_parts):
        '''Return distances to closest feature parts for intergenic locs
        '''
        dists = tuple(x[1] for x in feature_parts)
        intergenic = all(dist != 0 for dist in dists)
        # Report distance to closest feature(s) for intergenic locs
        if intergenic:
            return(dists)
        else:
            # Instead of returning None, do this to allow for case with empty .embl file
            return(tuple(set(dists)))

    # Get just the features for each loc, sorted by loc
    print('Adding annotations...')
    loc_features_sorted = [tuple(loc_to_features[loc]) for loc in hp.index]
    #
    # These are nested tuples, as qualifier values can take multiple values e.g. 'gene': ['dnaA', 'dnaH']
    hp['locus_tag'] = [get_qualifier_value(feature_parts, 'locus_tag') for feature_parts in loc_features_sorted]
    hp['gene'] = [get_qualifier_value(feature_parts, 'gene') for feature_parts in loc_features_sorted]
    hp['product'] = [get_qualifier_value(feature_parts, 'product') for feature_parts in loc_features_sorted]
    #
    hp['intergenic'] = [get_dists_for_intergenic(feature_parts) for feature_parts in loc_features_sorted]
    #
    hp['pseudo'] = [tuple(map(lambda x: '1' if 'pseudo' in x[0].qualifiers else 0, feature_parts)) for feature_parts in loc_features_sorted]
    # These are non-nested lists, as a feature has a single strand/type
    hp['strand'] = [tuple(map(lambda x: x[0].strand, feature_parts)) for feature_parts in loc_features_sorted]
    hp['type'] = [tuple(map(lambda x: x[0].type, feature_parts)) for feature_parts in loc_features_sorted]

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
        #  "st22": '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/CC22_EMRSA15.embl',
        "st22": '/nfs/users/nfs_b/bb9/workspace/rotation3/misc/2017-05-23_staph_annot_from_matt/EMRSA15_with_Mels_currated.modified.embl',
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
    #  for prefix in ['st22']:

        #  prefix = 'st239'

        hps = {}
        for t in transformations:
            print('Processing: {}, {}'.format(t, prefix))
            df = dfs[t].loc[prefix]
            print('{}: {} changes in .tab file.'.format(prefix, len(df)))
            hp = get_homoplasies(df, keep_all_sites=True)
            print('{}: {} sites in .tab file.'.format(prefix, len(hp)))
            print('{}: {} homoplasic sites in .tab file.'.format(prefix, len(hp.query('n_convergence > 0 | n_reversal > 0'))))
            #
            # Check the counts are partitioned right
            #  assert(all(hp[hp['n_convergence'] + hp['n_reversal'] - hp['n_both'] + hp['n_non_homoplasic'] == hp['n_total']]))
            #  assert(all(hp[hp['n_homoplasic'] + hp['n_non_homoplasic'] == hp['n_total']]))
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
        hp_merged['agree_acctran_deltran'] = (
            (hp_merged['n_convergence_acctran'] == hp_merged['n_convergence_deltran']) 
            & (hp_merged['n_reversal_acctran'] == hp_merged['n_reversal_deltran']) 
            #  & (hp_merged['n_both_acctran'] == hp_merged['n_both_deltran'])
        )
        print('n_convergence_acctran: {}, n_convergence_deltran: {}'.format(
            len(hp_merged.query('n_convergence_acctran > 0')),
            len(hp_merged.query('n_convergence_deltran > 0'))
        ))
        print('n_reversal_acctran: {}, n_reversal_deltran: {}'.format(
            len(hp_merged.query('n_reversal_acctran > 0')),
            len(hp_merged.query('n_reversal_deltran > 0'))
        ))
        print('Homoplasic sites with acctran/deltran reconstruction agreement: {}/{}'.format(
            sum(hp_merged.query('n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0')['agree_acctran_deltran']),
            len(hp_merged.query('n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0'))
        ))

        # Annotate sites
        hp_merged = annotate_homoplasies(hp_merged, embl_files[prefix])
        print('{}: {} within-feature annotations added.'.format(prefix, sum(hp_merged['intergenic'].apply(lambda x: x != (0, )))))

        #  Write out homoplasies
        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasies.csv")
        hp_merged[hp_merged['gene'].astype(bool)].sort_values(['n_homoplasic_acctran', 'n_homoplasic_deltran'], ascending=False).to_csv(outfile)

        #  def nested_lists_to_tuples(ll):
            #  '''Convert nested lists to nested tuples to allow hashing
            #  '''
            #  if type(ll) == list:
                #  return(tuple(nested_lists_to_tuples(l) for l in ll))
            #  else:
                #  return(ll)

        def get_ordered_unique_values(x):
            '''Return unique values in an iterable of hashable elements, retaining order
            '''
            return(tuple(collections.OrderedDict(itertools.zip_longest(x, (None, ))).keys()))

        #  Get genes with the highest ratio of convergences to total changes
        print('Summarising by feature...')
        #  summary_series = hp_merged['locus_tag'].apply(nested_lists_to_tuples)
        # Group intergenic and non intergenic separately
        summary_series = hp_merged.apply(lambda x: str((x['locus_tag'], 'intergenic')) if x['intergenic'] != (0, ) else str(x['locus_tag']), axis=1)
        #
        genes = pd.DataFrame()
        genes['intergenic'] = hp_merged.groupby(summary_series).apply(lambda x: any(dist != (0, ) for dist in x['intergenic']))
        genes['locus_tag'] = hp_merged.groupby(summary_series)['locus_tag'].apply(lambda x: get_ordered_unique_values(x))
        genes['gene'] = hp_merged.groupby(summary_series)['gene'].apply(lambda x: get_ordered_unique_values(x))
        genes['product'] = hp_merged.groupby(summary_series)['product'].apply(lambda x: get_ordered_unique_values(x))
        genes['pseudo'] = hp_merged.groupby(summary_series)['pseudo'].apply(lambda x: get_ordered_unique_values(x))
        genes['strand'] = hp_merged.groupby(summary_series)['strand'].apply(lambda x: get_ordered_unique_values(x))
        #
        # n_sites is the number of sites within the gene
        genes['n_sites'] = hp_merged.groupby(summary_series)['strand'].apply(len)
        genes['agree_acctran_deltran'] = hp_merged.groupby(summary_series).apply(lambda x: collections.Counter(x['agree_acctran_deltran']))
        genes['agree_prop_acctran_deltran'] = genes['agree_acctran_deltran'].apply(lambda x: x[True]/(x[True]+x[False]))
        #
        #  genes['n_convergence_acctran'] = hp_merged.groupby(summary_series)['n_convergence_acctran'].apply(sum)
        #  genes['n_convergence_deltran'] = hp_merged.groupby(summary_series)['n_convergence_deltran'].apply(sum)
        #  genes['n_reversal_acctran'] = hp_merged.groupby(summary_series)['n_reversal_acctran'].apply(sum)
        #  genes['n_reversal_deltran'] = hp_merged.groupby(summary_series)['n_reversal_deltran'].apply(sum)
        #  genes['n_both_acctran'] = hp_merged.groupby(summary_series)['n_both_acctran'].apply(sum)
        #  genes['n_both_deltran'] = hp_merged.groupby(summary_series)['n_both_deltran'].apply(sum)
        #  genes['n_homoplasic_acctran'] = hp_merged.groupby(summary_series)['n_homoplasic_acctran'].apply(sum)
        #  genes['n_homoplasic_deltran'] = hp_merged.groupby(summary_series)['n_homoplasic_deltran'].apply(sum)
        # n_total is number of changes at a site (summed for all sites in the gene), giving number of changes in a gene
        #  i.e. assert(all(hp[hp['n_convergence'] + hp['n_reversal'] - hp['n_both'] + hp['n_non_homoplasic'] == hp['n_total']]))
        #  genes['n_total'] = hp_merged.groupby(summary_series)['n_total_acctran'].apply(sum)
        #
        # Get number of sites with homoplasies in a gene
        #
        genes['n_convergence_sites_acctran'] = hp_merged.groupby(summary_series)['n_convergence_acctran'].apply(lambda x: x.astype(bool).sum())
        genes['n_convergence_sites_deltran'] = hp_merged.groupby(summary_series)['n_convergence_deltran'].apply(lambda x: x.astype(bool).sum())
        genes['n_reversal_sites_acctran'] = hp_merged.groupby(summary_series)['n_reversal_acctran'].apply(lambda x: x.astype(bool).sum())
        genes['n_reversal_sites_deltran'] = hp_merged.groupby(summary_series)['n_reversal_deltran'].apply(lambda x: x.astype(bool).sum())
        #  genes['n_homoplasic_sites_acctran'] = hp_merged.groupby(summary_series)['n_homoplasic_acctran'].apply(lambda x: x.astype(bool).sum())
        #  genes['n_homoplasic_sites_deltran'] = hp_merged.groupby(summary_series)['n_homoplasic_deltran'].apply(lambda x: x.astype(bool).sum())
        ##
        #
        #  genes['n_h_acctran/n_t'] = genes['n_homoplasic_acctran']/genes['n_total']
        #  genes['n_h_deltran/n_t'] = genes['n_homoplasic_deltran']/genes['n_total']
        #
        genes['n_homoplasic_sites'] = hp_merged.groupby(summary_series)['n_homoplasic_acctran'].apply(lambda x: x.astype(bool).sum())
        assert(all(
            genes['n_homoplasic_sites'] ==
            hp_merged.groupby(summary_series)['n_homoplasic_deltran'].apply(lambda x: x.astype(bool).sum())
        ))
        #
        genes['n_hp_sites/n_sites'] = genes['n_homoplasic_sites']/genes['n_sites']

        # Write out summary
        #
        #  Why are the tuples so deeply nested?
        #  single feature: tuple, feature can have multiple qualifiers
        #  single location (hp): tuple of tuples, a location can have overlapping features
        #  grouping of locations in features (genes): tuple of tuple of tuples, a single feature contains multiple locations
        genes.index.name = 'feature'
        outfile = os.path.join(out_file_prefixes[prefix], prefix + "_homoplasies_per_gene.csv")
        #  genes.sort_values(['n_h_acctran/n_t', 'n_h_deltran/n_t'], ascending=False).to_csv(outfile)
        genes.sort_values(['n_hp_sites/n_sites', 'n_sites'], ascending=False).to_csv(outfile)

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

