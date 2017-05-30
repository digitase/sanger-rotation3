#
# Run iprscan on all CDSs for all STs
#

IPR_DIR="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/packages/interproscan-5.23-62.0/"
# Meet installation requirements
PY27_DIR="/software/pathogen/external/apps/usr/local/Python-2.7.13/bin/"
export JAVA_HOME="/nfs/users/nfs_b/bb9/packages/jre1.8.0_121/"
export PATH="$JAVA_HOME/bin/:$PY27_DIR:$PATH"

outdir=".output/interproscan_all_array/"
mkdir -p "$outdir"
outdir="$(readlink -f $outdir)"

# Split locus tag from rest of header line
# stop_symbol - Single character string, what to use for terminators. This defaults to the asterisk, "*".
cat "/lustre/scratch118/infgen/team81/bb9/3_find_candidates/blast/fasta/"st*.cds.pep.fasta | \
    sed '/^>SAEMRSA15\|^>SAR/ s/_/ /' | \
    sed '/>SATW20_\|^>SAUSA300_/ s/_/ /2' | \
    sed 's/*//g' \
    > "$outdir/st_all.cds.pep.fasta"

mkdir -p "$outdir/fasta_chunks"
cat "$outdir/st_all.cds.pep.fasta" | parallel --pipe -N200 --recstart '>' "cat > $outdir/fasta_chunks/st_all.cds.pep.fasta.{#}"
n_chunks=$(ls "$outdir/fasta_chunks"/st_all.cds.pep.fasta.* | wc -l)
echo Running $n_chunks chunks.

# Excluded analyses:
# Gene3D
# bin/gene3d/4.1.0/cath-resolve-hits: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.18' not found (required by bin/gene3d/4.1.0/cath-resolve-hits)`

mkdir -p "$outdir/interproscan_out/"
mkdir -p "$outdir/logs/"
bsub \
    -G team81 -q normal \
    -n 2 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" -M 4000 \
    -J "interproscan_all_array[1-$n_chunks]" \
    -o "$outdir/logs/jobid_%J.%I.bsub_o.log" \
    -e "$outdir/logs/jobid_%J.%I.bsub_e.log" \
        "bash _3.0_run_iprscan_all_helper.sh $IPR_DIR $outdir"

# NOTE: analyses below run 1 st at a time (around 2.5k loci per ST)
exit

    # -R "select[mem>8000] rusage[mem=8000]" -M 8000 \
bsub \
    -G team81 -q long \
    -n 4 -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" -M 24000 \
    -J "interproscan_all_st22" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        "$IPR_DIR/interproscan.sh" \
            -i "$outdir/st22.cds.pep.fasta" \
            --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
            --goterms --pathways \
            --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
            --output-file-base "$outdir/interproscan_out_st22"

bsub \
    -G team81 -q long \
    -n 4 -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" -M 24000 \
    -J "interproscan_all_st239" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        "$IPR_DIR/interproscan.sh" \
            -i "$outdir/st239.cds.pep.fasta" \
            --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
            --goterms --pathways \
            --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
            --output-file-base "$outdir/interproscan_out_st239"

bsub \
    -G team81 -q long \
    -n 4 -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" -M 24000 \
    -J "interproscan_all_st8" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        "$IPR_DIR/interproscan.sh" \
            -i "$outdir/st8.cds.pep.fasta" \
            --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
            --goterms --pathways \
            --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
            --output-file-base "$outdir/interproscan_out_st8"

# TODO what a mess

