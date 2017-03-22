#!/bin/bash

export PATH="/usr/bin/:$PATH"

#-a: whole genome alignment file
#-t: tree file
#-p: output file prefix
#-e: EMBL annotation file
#-r: name of reference sequence related to EMBL file

outdir="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/st8"
mkdir -p "$outdir"

aln="/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/S.aureus_ST8_BSAC_Pfizer_no_mge.aln"
tree="/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/RAxML_bestTree.ml_S.aureus_ST8_BSAC_Pfizer"
prefix="st8_bsac_pfizer"
embl="/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl"
reference="CC8_USA300_FPR3757"

#
cd "$outdir"
~sh16/scripts/reconstruct_snps_on_tree.py \
    -a "$aln" \
    -t "$tree" \
    -p "$prefix" \
    -e "$embl" \
    -r "$reference" \
    --dNdS
cd -
    
#2g memory is plenty
