#!/usr/bin/env python
'''Convert tab files to csv format
'''

import sys 
import os
import errno
sys.path.append(os.path.abspath("/nfs/users/nfs_b/bb9/workspace/rotation3/src/scripts/"))
import parse_ft_tab_file

transformations = ["acctran", "deltran"]
tabfiles = [
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/{}/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_homoplasies_on_tree.tab",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/{}/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_homoplasies_on_tree.tab",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/{}/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_homoplasies_on_tree.tab",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/{}/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_homoplasies_on_tree.tab"
]

out_dir = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/"

for t in transformations:
    for f in tabfiles:
        #  Parse tab file
        tabfile = parse_ft_tab_file.FtTabFile()
        f = f.format(t)
        print("Parsing... " + f)
        tabfile.parse_file(f)
        # Make outdir
        st_name = os.path.basename(os.path.dirname(f))
        try:
            os.makedirs(os.path.join(out_dir, st_name))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # Write out csv
        f_base_prefix = os.path.splitext(os.path.basename(f))[0]
        tabfile.write_csv(os.path.join(out_dir, st_name, f_base_prefix + ".{}.csv".format(t)))

