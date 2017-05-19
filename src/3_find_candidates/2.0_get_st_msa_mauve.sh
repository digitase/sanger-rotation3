#
#
# Using mauve, construct XMFA alignment of reference ST seqs from each of the within ST MSAs
#

declare -a alns=(
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st22/S.aureus_ST22_BSAC_Pfizer_no_mge.aln' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/data/S.aureus_ST239_global_Singapore_Pfizer_no_mge.deduped.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st30/S.aureus_ST30_BSAC_Pfizer_no_mge.aln' \
    '/lustre/scratch118/infgen/team81/dj9/staph.aureus/ben/st8/revised_lane_list/S.aureus_ST8_BSAC_Pfizer_revised_no_mge.aln'
)
declare -a ref_prefixes=(
    'CC22_EMRSA15' \
    'TW20' \
    'CC30_MRSA252' \
    'CC8_USA300_FPR3757'
)
declare -a refs=(
    '/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.fasta' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.dna' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.fasta' \
    '/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.dna'
)

outdir='.output/mauveAligner/'
mkdir -p "$outdir"

refs_fa="$outdir/refs.fa"
cat ${refs[@]} > "$refs_fa"

# mauveAligner [options] <seq1 filename> <sml1 filename> ...  <seqN filename> <smlN filename>
# Options:
            # --output=<file> Output file name.  Prints to screen by default
            # --mums Find MUMs only, do not attempt to determine locally collinear blocks (LCBs)
            # --no-recursion Don't perform recursive anchor identification (implies --no-gapped-alignment)
            # --no-lcb-extension If determining LCBs, don't attempt to extend the LCBs
            # --seed-size=<number> Initial seed match size, default is log_2( average seq. length )
            # --max-extension-iterations=<number> Limit LCB extensions to this number of attempts, default is 4
            # --eliminate-inclusions Eliminate linked inclusions in subset matches.
            # --weight=<number> Minimum LCB weight in base pairs per sequence
            # --match-input=<file> Use specified match file instead of searching for matches
            # --lcb-match-input  Indicates that the match input file contains matches that have been clustered into LCBs
            # --lcb-input=<file> Use specified lcb file instead of constructing LCBs (skips LCB generation)
            # --scratch-path=<path>  For large genomes, use a directory for storage of temporary data.  Should be given two or more times to with different paths.
            # --id-matrix=<file> Generate LCB stats and write them to the specified file
            # --island-size=<number> Find islands larger than the given number
            # --island-output=<file> Output islands the given file (requires --island-size)
            # --backbone-size=<number> Find stretches of backbone longer than the given number of b.p.
            # --max-backbone-gap=<number> Allow backbone to be interrupted by gaps up to this length in b.p.
            # --backbone-output=<file> Output islands the given file (requires --island-size)
            # --coverage-output=<file> Output a coverage list to the specified file (- for stdout)
            # --repeats Generates a repeat map.  Only one sequence can be specified
            # --output-guide-tree=<file> Write out a guide tree to the designated file
            # --collinear Assume that input sequences are collinear--they have no rearrangements

# Gapped alignment controls:
            # --no-gapped-alignment Don't perform a gapped alignment
            # --max-gapped-aligner-length=<number> Maximum number of base pairs to attempt aligning with the gapped aligner
            # --min-recursive-gap-length=<number> Minimum size of gaps that Mauve will perform recursive MUM anchoring on (Default is 200)

# Signed permutation matrix options:
            # --permutation-matrix-output=<file> Write out the LCBs as a signed permutation matrix to the given file
            # --permutation-matrix-min-weight=<number> A permutation matrix will be written for every set of LCBs with weight between this value and the value of --weight

# Alignment output options:
            # --alignment-output-dir=<directory> Outputs a set of alignment files (one per LCB) to a given directory
            # --alignment-output-format=<directory> Selects the output format for --alignment-output-dir
            # --output-alignment=<file> Write out an XMFA format alignment to the designated file

# Supported alignment output formats are: phylip, clustal, msf, nexus, mega, codon

bsub \
    -G team81 -q long \
    -R "select[mem>2000] rusage[mem=2000]" -M 2000 \
    -J "mauveAligner.$prefix" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        mauveAligner \
            --output="$refs_fa.mauve" \
            --output-alignment="$refs_fa.alignment.xmfa" \
            "$refs_fa"

