#!/bin/bash
#
# Plot annotation in tab file against tree 
#

# tab files from reconstruct_snps_on_tree.py
declare -a tabs=(
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_homoplasies_on_tree.tab" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_homoplasies_on_tree.tab" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_homoplasies_on_tree.tab" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_homoplasies_on_tree.tab" 
)
declare -a embls=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/CC22_EMRSA15.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st239/CC8_TW20.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/CC30_MRSA252.embl' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/CC8_USA300_FPR3757.embl'
)
# acctran trees
declare -a trees=(
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST22_BSAC_Pfizer/S.aureus_ST22_BSAC_Pfizer_acctran_steps.tre" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST239_global_Singapore_Pfizer/S.aureus_ST239_global_Singapore_Pfizer_intergenic_acctran_steps.tre" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST30_BSAC_Pfizer/S.aureus_ST30_BSAC_Pfizer_acctran_steps.tre" 
    "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reconstruct_snps_on_tree/acctran/S.aureus_ST8_BSAC_Pfizer_revised/S.aureus_ST8_BSAC_Pfizer_revised_acctran_steps.tre" 
)
declare -a prefixes=(
    'S.aureus_ST22_BSAC_Pfizer' \
    'S.aureus_ST239_global_Singapore_Pfizer' \
    'S.aureus_ST30_BSAC_Pfizer' \
    'S.aureus_ST8_BSAC_Pfizer_revised'
)

#Reportlabtest also allows metadata to be added to trees from a csv file. The
#csv file must be comma separated, with the first column containing taxon names
#equivalent to those in the tree. The csv file should also contain a header row
#containing headings for each column of data. The metadata file is provided with
#the -m option. 

                # -b 2941669 -e 2944488 \
for (( i = 0; i < ${#trees[@]}; i++ )); do
    tree="${trees[$i]}"
    embl="${embls[$i]}"
    prefix="${prefixes[$i]}"
    tab="${tabs[$i]}"
    # Make st specific output dir
    outdir="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/reportlabtest/$prefix"
    mkdir -p "$outdir"
    #
    cd "$outdir"
    echo Starting $prefix at $(date)
    bsub \
        -G team81 -q normal \
        -R "select[mem>2000] rusage[mem=2000]" -M 2000 \
        -J "reportlabtest.$prefix" \
        -o "$outdir/jobid_%J.bsub_o.log" \
        -e "$outdir/jobid_%J.bsub_e.log" \
            python2 /nfs/users/nfs_b/bb9/workspace/rotation3/src/2_homoplasy/reportlabtest_modified.py \
                -t "$tree" \
                -q taxa \
                -l 2 -A 2 -E 5 \
                -a 2 -L left \
                -o "$prefix.homoplasies_on_acctran_steps_tree.pdf" \
                "$tab" "$embl"
    echo Finished $prefix at $(date)
done
    
