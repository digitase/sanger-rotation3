

#  Remove non homoplasic snps, recolor homoplasic snps by substitution

features = []

class Feature:
    qualifier = "FT"

    def __init__(self):
        self.type
        self.position
        self.props

    

tab_filename = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_homoplasies_on_tree.tab"
with open(tab_filename, 'r') as tab_file:
    for line in tab_file:
        # parse feature table
        if line.startwith("FT"):
            if 





#  FT   SNP             2592254
#  FT                   /strand="reverse"
#  FT                   /node="root->608"
#  FT                   /SNP="G->A"
#  FT                   /homoplasy="reversed in branch leading to 12971_1#50"
#  FT                   /codon_type="Synonymous"
#  FT                   /colour=2
#  FT                   /taxa="12673_2#62, 12625_5#43, 12673_3#89, 12971_2#49, 14672_6#18, 14324_2#58, 14355_3#4, 7521_5#67, 14355_3#30, 14355_3#38, 12673_2#43, 12641_3#52, 12673_4#75, 7748_6#16, 7564_8#14, 12971_2#34, 12673_2#29, 14672_7#33, 14672_6#41, 14672_8#86, 14672_8#92, 14672_5#85, 14672_5#68, 14672_6#89, 14672_5#72, 14672_7#20, 14672_5#75, 14672_8#44, 14324_3#48, 7480_8#69, 14672_7#51, 14672_8#60, 14672_8#85, 14672_7#66, 7748_6#43, 14672_6#63, 14672_7#31, 14672_6#51"
