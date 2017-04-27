#!/nfs/users/nfs_b/bb9/miniconda3/envs/py35/bin/python
'''Plot the tip states for homoplasic sites in and around genes against the tree.
'''

import collections
import itertools
import functools
import pandas as pd
import sys
import os
import tempfile
from Bio import SeqIO
sys.path.append(os.path.abspath("/nfs/users/nfs_b/bb9/workspace/rotation3/src/scripts/"))
import parse_ft_tab_file

reportlabtest_script = '/nfs/users/nfs_b/bb9/workspace/rotation3/src/scripts/reportlabtest_modified.py'
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
# Alternate acctran/deltran tree
#  transformations = ['acctran', 'deltran']
#  tree_file_template = "/lustre/scratch118/infgen/team81/bb9/2_homoplasy/reconstruct_snps_on_tree/{transformation}/{prefix}/{prefix}_{transformation}_steps.tre"
hp_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/{prefix}/{short_prefix}_homoplasies.csv"
embl_files = {
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}
#  out_tab_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_verify_candidates/plots/{prefix}/{prefix}.locus_tag_{locus_tag}.tab"
out_pdf_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_verify_candidates/plots_reportlabtest/homologs/{prefix}/{prefix}.locus_tag_{locus_tag}.pdf"
# Dir to store intermediate files in
tmp_dir = '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/'

# NOTE: make sure there are no spaces in the locus tags
#
# Top 3 ranked loci for each st
#  to_plot = {
    #  'st22': ['SAEMRSA1525350', 'SAEMRSA1504890', 'SAEMRSA1517670'],
    #  'st239': ['SATW20_27680', 'SATW20_26650', 'SATW20_06320'],
    #  'st30': ['SAR1827', 'SAR2709', 'SAR0375'],
    #  'st8': ['SAUSA300_0148', 'SAUSA300_1981', 'SAUSA300_2565']
#  }
#

# Loci in top 10 homologous sets
to_plot_str = '''
('SAR0370', 'SATW20_08950')
('SAR1827', 'SATW20_27080')
('SAR0373', 'SAUSA300_0808')
('SAR0372', 'SAUSA300_0807')
('SAEMRSA1519360', 'SAR0381', 'SATW20_09070', 'SATW20_20110')
('SAR0378', 'SATW20_09030')
('SAR1826', 'SATW20_27070')
('SAR0377', 'SATW20_09020')
('SAEMRSA1519350', 'SAR0382', 'SATW20_09080', 'SATW20_20100', 'SAUSA300_1981')
('SAEMRSA1512890', 'SAR1439', 'SATW20_14270', 'SAUSA300_1319')
'''
# Lazy way of plotting copypastas from .csv files open in Excel
to_plot_tuple = functools.reduce(lambda x, y: x+y, [eval(x) for x in to_plot_str.strip().split('\n')])
#
to_plot = {
    'st22' : [x for x in to_plot_tuple if x.startswith('SAEMRSA')],
    'st239': [x for x in to_plot_tuple if x.startswith('SATW20_')],
    'st30' : [x for x in to_plot_tuple if x.startswith('SAR')],
    'st8'  : [x for x in to_plot_tuple if x.startswith('SAUSA300_')],
}

#
print(to_plot)
for short_prefix, locus_tags in to_plot.items():

    #  short_prefix = 'st8'
    #  locus_tags = ['SAUSA300_2565']
    print((short_prefix, locus_tags))

    #  Read in snp state for each taxon
    print('{}: Reading in snps...'.format(short_prefix))
    snp_summary = pd.read_table(snp_summary_file_template.format(prefix=prefixes[short_prefix]))
    snp_summary.columns = map(str.strip, snp_summary.columns)

    # Create taxa by snps matrix
    snps = snp_summary.iloc[:, 10:].transpose()
    snps.index.name = 'taxon'
    snps.columns = snp_summary.iloc[:, 1]
    snps.columns.name = 'loc'

    # Read in closest locus for each loc
    print('{}: Getting locus bounds...'.format(short_prefix))
    hp = pd.read_csv(hp_file_template.format(prefix=prefixes[short_prefix], short_prefix=short_prefix))
    # Check that all sites in hp appear in the summary
    assert(all(hp['loc'].isin(snps.columns)))

    # Get start and end coords of each locus
    embl_record = next(SeqIO.parse(embl_files[short_prefix], 'embl'))
    locus_tag_to_bounds = {}
    for feature in embl_record.features:
        if feature.type == 'CDS':
            b, e = feature.location.start.position, feature.location.end.position
            locus_tag_to_bounds[feature.qualifiers['locus_tag'][0]] = (b, e)

    for locus_tag in locus_tags:

        #  locus_tag = 'SAUSA300_2565'

        # Keep only locs closest to each locus 
        print('{}: {}: Gathering closest snps...'.format(short_prefix, locus_tag))
        locs_to_keep = hp[hp['locus_tag'].apply(lambda x: locus_tag in x)]['loc']
        snps_filtered = snps.copy().loc[:, snps.columns.isin(locs_to_keep)]
        # Convert to long form
        snps_filtered['taxon'] = snps.index
        meta = pd.melt(snps_filtered, id_vars=['taxon'], value_name='snp')

        # Fill in reference base
        print('{}: {}: Filling in ref bases...'.format(short_prefix, locus_tag))
        loc_to_ref_base = dict(zip(snp_summary.iloc[:, 1], snp_summary['Ref_base']))
        meta['snp_full'] = meta['snp']
        # Check if there are any snps associated with the locus at all
        if len(meta):
            meta.loc[meta['snp'] == '.', 'snp_full'] = meta.loc[meta['snp'] == '.'].apply(lambda row: loc_to_ref_base[row['loc']], axis=1)
        else:
            print('~~~ WARNING: no sites associated with locus_tag: {} ~~~'.format(locus_tag))

        # Mark homoplasic sites
        print('{}: {}: Marking homoplasic sites...'.format(short_prefix, locus_tag))
        loc_to_homoplasic = dict(zip(hp['loc'], hp['n_homoplasic_acctran'] + hp['n_homoplasic_deltran'] > 0))
        #  meta.loc[:, 'homoplasic'] = meta.apply(lambda row: loc_to_homoplasic[row['loc']], axis=1)
        meta.loc[:, 'homoplasic'] = meta['loc'].map(lambda x: loc_to_homoplasic[x])

        # Colour loc based on whether a loc is homoplasic
        print('{}: {}: Colouring locs...'.format(short_prefix, locus_tag))
        #  colours = {
            #  ('.', False): 'white', 
            #  ('N', False): 'white', 
            #  ('A', False): 'lightgreen', 
            #  ('C', False): 'lightblue', 
            #  ('G', False): 'gray', 
            #  ('T', False): 'lightred',
            #  ('.', True): 'white', 
            #  ('N', True): 'white', 
            #  ('A', True): 'green', 
            #  ('C', True): 'blue', 
            #  ('G', True): 'black', 
            #  ('T', True): 'red'
        #  }
        colours = {
            ('.', False): 13, 
            ('N', False): 13, 
            ('A', False): 13, 
            ('C', False): 13, 
            ('G', False): 13, 
            ('T', False): 13,
            ('.', True): 1, 
            ('N', True): 13, 
            ('A', True): 3, 
            ('C', True): 4, 
            ('G', True): 14, 
            ('T', True): 2 
        }
        #
        meta.loc[:, 'colour'] = meta.apply(lambda x: colours[(x['snp_full'], x['homoplasic'])], axis=1, reduce=True)

        # Generate .tab file features
        def meta_row_to_tab_ft(x):
            '''Convert metadata table row to .tab file feature.
            '''
            ft = parse_ft_tab_file.Feature(
                    key='SNP',
                    loc=x['loc'], 
                    qualifs={'taxa': x['taxon'], 'colour': x['colour']}
            )
            return(ft)
        #
        print('{}: {}: Generating .tab file...'.format(short_prefix, locus_tag))
        ft_tab_file = parse_ft_tab_file.FtTabFile()
        ft_tab_file.features = list(meta.apply(meta_row_to_tab_ft, axis=1, reduce=True))
        
        #  ft_tab_file.write_tab(out_tab_template.format(prefix=prefixes[short_prefix], locus_tag=locus_tag))
        # Temporary .tab file for reportlabtest
        # Must end in .tab to be recognised.
        _, tmp_tab_file = tempfile.mkstemp(suffix='.{}.{}.reportlabtest.tab'.format(short_prefix, locus_tag), dir=tmp_dir)
        ft_tab_file.write_tab(tmp_tab_file)

        try:
            print('{}: {}: Running reportlabtest...'.format(short_prefix, locus_tag))

            #  bsub = ('', )
            bsub = (r'bsub -G team81 -q normal -R "select[mem>2000] rusage[mem=2000]" -M 2000 -J "reportlabtest.{short_prefix}.{locus_tag}" -o "/dev/null" -e "/dev/null"'.format(short_prefix=short_prefix, locus_tag=locus_tag), )

            out_pdf_dir, out_pdf = os.path.split(out_pdf_template.format(prefix=prefixes[short_prefix], locus_tag=locus_tag))
            os.makedirs(out_pdf_dir, exist_ok=True)

            start = min(locus_tag_to_bounds[locus_tag][0], min(meta['loc'])) if len(meta['loc']) else locus_tag_to_bounds[locus_tag][0]
            end = max(locus_tag_to_bounds[locus_tag][1], max(meta['loc'])) if len(meta['loc']) else locus_tag_to_bounds[locus_tag][1]

            cmd = ' '.join(
                bsub + (
                    r'python2 {reportlabtest_script}',
                    r'-t "{tree}"',
                    r'-q taxa',
                    r'-b {start} -e {end}',
                    r'-l 2 -A 2 -E 15',
                    r'-a 2 -L left --proportion 0.4',
                    r'-o "{pdf}"',
                    r'"{tab}" "{embl}"'
                )
            ).format(
                reportlabtest_script=reportlabtest_script,
                tree=tree_files[short_prefix],
                start=start,
                end=end, 
                #  pdf='test.pdf',
                pdf=os.path.join(out_pdf_dir, out_pdf),
                #  tab='test.tab', 
                #  tab=out_tab_template.format(prefix=prefixes[short_prefix], locus_tag=locus_tag),
                tab=tmp_tab_file,
                embl=embl_files[short_prefix]
            )
            #
            #  with open('test.sh', 'w') as test_fhandle:
                #  test_fhandle.write(cmd)
            # 

            print(cmd)
            _ = os.system(cmd)

        finally:
            # TODO Can't remove .tab files before bsub job completes...
            if not len(bsub[0]):
                os.remove(tmp_tab_file)

