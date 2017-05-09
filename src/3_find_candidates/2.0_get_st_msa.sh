#
#
# Using MAFFT, construct MSA of reference ST seqs from each of the within ST MSAs
#

declare -a alns=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/S.aureus_ST22_BSAC_Pfizer_no_mge.aln' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/data/S.aureus_ST239_global_Singapore_Pfizer_no_mge.deduped.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/S.aureus_ST30_BSAC_Pfizer_no_mge.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/revised_lane_list/S.aureus_ST8_BSAC_Pfizer_revised_no_mge.aln'
)
declare -a refs=(
    'CC22_EMRSA15' \
    'TW20' \
    'CC30_MRSA252' \
    'CC8_USA300_FPR3757'
)

outdir='.output/mafft/'
mkdir -p "$outdir"

refs_fa="$outdir/refs.fa"
echo -n > "$refs_fa"

for (( i = 0; i < ${#alns[@]}; i++ )); do
    aln="${alns[$i]}"
    ref="${refs[$i]}"
    # Get ref sequence for st from within st .msa
    fastagrep.sh "$ref" "$aln" >> "$refs_fa"
done

# ------------------------------------------------------------------------------
  # MAFFT v7.205 (2014/10/20)
  # http://mafft.cbrc.jp/alignment/software/
  # MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
# ------------------------------------------------------------------------------
# High speed:
  # % mafft in > out
  # % mafft --retree 1 in > out (fast)

# High accuracy (for <~200 sequences x <~2,000 aa/nt):
  # % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)
  # % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)
  # % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)

# If unsure which option to use:
  # % mafft --auto in > out

# --op # :         Gap opening penalty, default: 1.53
# --ep # :         Offset (works like gap extension penalty), default: 0.0
# --maxiterate # : Maximum number of iterative refinement, default: 0
# --clustalout :   Output: clustal format, default: fasta
# --reorder :      Outorder: aligned, default: input order
# --quiet :        Do not report progress
# --thread # :     Number of threads (if unsure, --thread -1)

bsub \
    -G team81 -q long \
    -n4 \
    -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" -M 4000 \
    -J "mafft.$prefix" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        "mafft --auto --thread 4 $refs_fa > $outdir/refs_msa.fa"

