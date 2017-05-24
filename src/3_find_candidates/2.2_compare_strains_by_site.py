'''Compare homoplasic sites across strains

Identify sites where homoplasies have occured in multiple strains
'''

import pandas as pd
import pickle

prefixes = {
    "st22": "S.aureus_ST22_BSAC_Pfizer",
    "st239": "S.aureus_ST239_global_Singapore_Pfizer",
    #  "st30": "S.aureus_ST30_BSAC_Pfizer",
    "st8": "S.aureus_ST8_BSAC_Pfizer_revised"
}
id_to_short_prefix = {
    "CC22_EMRSA15": "st22",
    "TW20": "st239",
    "CC30_MRSA252": "st30",
    "CC8_USA300_FPR3757": "st8"
}
hp_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/{prefix}/{short_prefix}_homoplasies.csv"
embl_files = {
    #  "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/misc/2017-05-23_staph_annot_from_matt/EMRSA15_with_Mels_currated.modified.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    #  "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}
#  snps_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_snps/{prefix}/{prefix}.out"

print('Reading hps...')
hps = dict((short_prefix, pd.read_csv(hp_file_template.format(prefix=prefixes[short_prefix], short_prefix=short_prefix))) for short_prefix in prefixes)
#  print('Reading snp summaries...')
#  snps = dict((short_prefix, pd.read_table(snps_file_template.format(prefix=prefixes[short_prefix], short_prefix=short_prefix))) for short_prefix in prefixes)

print('Mapping st loc to alignment loc...')
#  pos_to_alignment_pos = {}
#  for short_prefix in prefixes:
    #  short_prefix = 'st22'
    #  print(short_prefix)
    #  hp = hps[short_prefix]
    #  snp = snps[short_prefix]
    #  pos_to_alignment_pos[short_prefix] = dict(zip(snp.iloc[:, 1], snp.iloc[:, 0]))
    #  hp['loc_in_alignment'] = hp['loc'].apply(lambda loc: pos_to_alignment_pos[short_prefix][loc])
# Read in mapping
with open('.output/pos_to_pos_aln.mauveAligner.pk', 'rb') as pk_file:
    pos_to_alignment_pos = pickle.load(pk_file)
# Rename seq ids to short_prefixes
ids = list(pos_to_alignment_pos.keys())
for k in ids:
    pos_to_alignment_pos[id_to_short_prefix[k]] = pos_to_alignment_pos.pop(k)
# Apply mapping
for short_prefix in prefixes:
    hp = hps[short_prefix]
    hp['loc_in_alignment'] = hp['loc'].apply(lambda loc: pos_to_alignment_pos[short_prefix][loc])

print('Merging st datasets...')
hps_merged = pd.concat(hps.values(), keys=hps.keys()).reset_index()

# Filter out non-homplasic sites
print('Filtering out non homoplasic sites...')
hps_merged = hps_merged.loc[hps_merged['n_homoplasic_acctran'].astype(bool) | hps_merged['n_homoplasic_deltran'].astype(bool)]

print('Writing output...')
dupe_sites = hps_merged[hps_merged.duplicated('loc_in_alignment', keep=False)].sort_values(['loc_in_alignment', 'level_1'])
dupe_sites.to_csv('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/3_find_candidates/hp_sites_summary.csv')

