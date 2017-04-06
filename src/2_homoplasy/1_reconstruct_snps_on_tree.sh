#!/bin/bash
#
# Map SNPs to tree per st

#-a: whole genome alignment file
#Use st239 aln without duplicate of ref seq
declare -a alns=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/S.aureus_ST22_BSAC_Pfizer_no_mge.aln' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/data/S.aureus_ST239_global_Singapore_Pfizer_no_mge.deduped.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/S.aureus_ST30_BSAC_Pfizer_no_mge.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/revised_lane_list/S.aureus_ST8_BSAC_Pfizer_revised_no_mge.aln'
)
#-t: tree file
declare -a trees=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/RAxML_bestTree.ml_S.aureus_ST22_BSAC_Pfizer' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st239/RAxML_bestTree.ml_S.aureus_ST239_global_Singapore_Pfizer' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/RAxML_bestTree.ml_S.aureus_ST30_BSAC_Pfizer' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/revised_lane_list/RAxML_bestTree.ml_S.aureus_ST8_BSAC_Pfizer_revised'
)
#-p: output file prefix
declare -a prefixes=(
    'S.aureus_ST22_BSAC_Pfizer' \
    'S.aureus_ST239_global_Singapore_Pfizer' \
    'S.aureus_ST30_BSAC_Pfizer' \
    'S.aureus_ST8_BSAC_Pfizer_revised'
)
#-e: EMBL annotation file
declare -a embls=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/CC22_EMRSA15.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st239/CC8_TW20.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/CC30_MRSA252.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/CC8_USA300_FPR3757.embl'
)
#-r: name of reference sequence related to EMBL file
#just the .id field of alignment object, the first whitespace delimited part
declare -a refs=(
    'CC22_EMRSA15' \
    'TW20' \
    'CC30_MRSA252' \
    'CC8_USA300_FPR3757'
)

# transformation="acctran"
transformation="deltran"
# transformation="ML"
#
for (( i = 0; i < ${#alns[@]}; i++  )); do
    aln="${alns[$i]}"
    tree="${trees[$i]}"
    prefix="${prefixes[$i]}"
    embl="${embls[$i]}"
    ref="${refs[$i]}"
    # Make st specific output dir
    outdir="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/$transformation/$prefix"
    mkdir -p "$outdir"
    #
    cd "$outdir"
    echo Starting $prefix at $(date)
    bsub \
        -G team81 -q long \
        -R "select[mem>16000] rusage[mem=16000]" -M 16000 \
        -J "reconstruct_snps_on_tree.$prefix" \
        -o "$outdir/jobid_%J.bsub_o.log" \
        -e "$outdir/jobid_%J.bsub_e.log" \
            "python2 ~sh16/scripts/reconstruct_snps_on_tree.py -a $aln -t $tree -p $prefix -T $transformation -e $embl -r $ref"
    cd -
    echo Finished bsub of $prefix at $(date)
done

