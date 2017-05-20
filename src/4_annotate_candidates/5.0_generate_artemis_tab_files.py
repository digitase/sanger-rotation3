'''Write out locations of sites with homoplasic SNPs in multiple lineages, using coords for each lineage
'''

import pandas as pd
import sys, os
import collections
#  sys.path.append(os.path.abspath("/nfs/users/nfs_b/bb9/workspace/rotation3/src/scripts/"))
#  import parse_ft_tab_file

summary_sites_by_gene_filtered = pd.read_csv('.output/hp_sites_summary.filtered.annotated.csv', converters={'sts': eval, 'loc_in_st': eval})

# For each strain, get locations of hp SNPs within that strain

tabfiles = {
    'st22': "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_snps_on_tree.tab",
    'st239': "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_snps_on_tree.tab",
    #  'st30': "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_snps_on_tree.tab",
    'st8': "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_snps_on_tree.tab"
}

#  Get locations of hp snps in each st
st_snp_locs = collections.defaultdict(set)
st_snp_to_hp_sts = {}
for st in tabfiles:
    for sts, locs in zip(summary_sites_by_gene_filtered['sts'], summary_sites_by_gene_filtered['loc_in_st']):
        for snp_st, loc in zip(sts, locs):
           if st == snp_st: 
                st_snp_locs[st].add(loc)
                st_snp_to_hp_sts[(st, loc)] = sorted(sts)

#  Filter tabfiles for hp snps
for st, tabfile in tabfiles.items():
    output_lines = []
    with open(tabfile, 'r') as tabfile_handle:
        for line in tabfile_handle:
            if line.startswith('FT   SNP'):
                ft_snp_loc = int(line.split()[2])
                if ft_snp_loc in st_snp_locs[st]:
                    output_lines.append(line)
                    #  Remove the loc to avoid adding duplicate locations
                    st_snp_locs[st].remove(ft_snp_loc)
                    #  Process the rest of the entrie's lines
                    #
                    #  Example entry format
                    #  FT   SNP             2737412
                    #  FT                   /strand="reverse"
                    #  FT                   /node="1655->1657"
                    #  FT                   /SNP="T->C"
                    #  FT                   /homoplasy="reversal from branch leading to 1654"
                    #  FT                   /codon_type="Synonymous"
                    #  FT                   /colour=3
                    #  FT                   /taxa="7414_8#71, 7414_8#51, 7521_6#60"
                    while True:
                        line = next(tabfile_handle)
                        #  Skip these lines
                        if '/node=' in line  or '/taxa=' in line: 
                            continue
                        #  Simplify line
                        elif '/homoplasy' in line:
                            output_lines.append('FT                   /homoplasy=1\n')
                        #  Add a qualifier with the sts that the snp is homoplasic in 
                        elif line.startswith('FT   SNP'): 
                            output_lines.append('FT                   /homoplasy_sts="{}"\n'.format','.join(st_snp_to_hp_sts[(st, ft_snp_loc)])))
                            break
                        else:
                            output_lines.append(line)
    # 
    with open('.output/hp_sites_summary.filtered.annotated.{}.tab'.format(st), 'w') as outfile:
        for line in output_lines:
            outfile.write(line)

