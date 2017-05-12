



IPR_DIR="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/packages/interproscan-5.23-62.0/"
# Meet installation requirements
PY27_DIR="/software/pathogen/external/apps/usr/local/Python-2.7.13/bin/"
export JAVA_HOME="/nfs/users/nfs_b/bb9/packages/jre1.8.0_121/"
export PATH="$JAVA_HOME/bin/:$PY27_DIR:$PATH"

outdir="$(readlink -f .output/interproscan/)"
mkdir -p "$outdir"

# Extract peptide seqs of unannotated loci
#
unannotated_loci=".output/locus_tags_unannotated.txt"
# Split locus tag from rest of header line
# stop_symbol - Single character string, what to use for terminators. This defaults to the asterisk, "*".
cat "/lustre/scratch118/infgen/team81/bb9/3_find_candidates/blast/fasta/"st*.cds.pep.fasta | \
    sed '/^>SAEMRSA15\|^>SAR/ s/_/ /' | \
    sed '/>SATW20_\|^>SAUSA300_/ s/_/ /2' | \
    sed 's/*//g' \
    > "$outdir/st_all.cds.pep.fasta"
#
fasta_extractor.pl -i "$outdir/st_all.cds.pep.fasta" -s "$outdir/st_all.cds.pep.unannotated.fasta" -f "$unannotated_loci"
wc -l "$unannotated_loci"

# Excluded analyses:
# Gene3D
# bin/gene3d/4.1.0/cath-resolve-hits: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.18' not found (required by bin/gene3d/4.1.0/cath-resolve-hits)`

bsub \
    -G team81 -q long \
    -R "select[mem>8000] rusage[mem=8000]" -M 8000 \
    -J "interproscan_unannotated" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        "$IPR_DIR/interproscan.sh" \
            -i "$outdir/st_all.cds.pep.unannotated.fasta" \
            --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
            --goterms --pathways \
            --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
            --output-file-base "$outdir/interproscan_out"

