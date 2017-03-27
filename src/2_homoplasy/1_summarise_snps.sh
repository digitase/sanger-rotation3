#!/bin/bash
#
# Get snp summary from alignment
#

declare -a alns=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/S.aureus_ST22_BSAC_Pfizer_no_mge.aln' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/data/S.aureus_ST239_global_Singapore_Pfizer_no_mge.deduped.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/S.aureus_ST30_BSAC_Pfizer_no_mge.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/revised_lane_list/S.aureus_ST8_BSAC_Pfizer_revised_no_mge.aln'
)
declare -a embls=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/CC22_EMRSA15.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st239/CC8_TW20.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/CC30_MRSA252.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/CC8_USA300_FPR3757.embl'
)
declare -a prefixes=(
    'S.aureus_ST22_BSAC_Pfizer' \
    'S.aureus_ST239_global_Singapore_Pfizer' \
    'S.aureus_ST30_BSAC_Pfizer' \
    'S.aureus_ST8_BSAC_Pfizer_revised'
)
declare -a refs=(
    'CC22_EMRSA15' \
    'TW20' \
    'CC30_MRSA252' \
    'CC8_USA300_FPR3757'
)

for (( i = 0; i < ${#alns[@]}; i++ )); do
    aln="${alns[$i]}"
    prefix="${prefixes[$i]}"
    embl="${embls[$i]}"
    ref="${refs[$i]}"
    # Make st specific output dir
    outdir="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_snps/$prefix"
    mkdir -p "$outdir"
    #
    cd "$outdir"
    echo Starting $prefix at $(date)
    bsub \
        -G team81 -q long \
        -R "select[mem>12000] rusage[mem=12000]" -M 12000 \
        -J "summarise_snps.$prefix" \
        -o "$outdir/jobid_%J.bsub_o.log" \
        -e "$outdir/jobid_%J.bsub_e.log" \
            python2 /nfs/users/nfs_s/sh16/scripts/summarise_snps.py \
                -i "$aln" -r "$ref" \
                -e "$embl" -o "$prefix" -t
    echo Finished $prefix at $(date)
done

#summarise_snps_from_aln.py Usage:
#summarise_snps_from_aln.py [options] <input alignment file>
#Options:
#-r              name of reference strain
#-e <filename>   embl file of reference strain [optional]
#-x              produce moving window snp plot and chi-squared test plot
#-c              produce recombination hotspot tab files for each sequence
#-o              prefix for output file names
#-t              produce tab files of snp locations
#-a              produce alignment of snp locations
#-g              produce mapping graph files
#-p              run phylogeny of snp locations using RAxML
#-m <model>      model of evolution for phylogeny [GTRGAMMA/GTRGAMMAI/GTRCAT]
#-b <int>        number of bootstrap replicates [0 = do not run bootstrap]
#-h              show this help
#Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009

