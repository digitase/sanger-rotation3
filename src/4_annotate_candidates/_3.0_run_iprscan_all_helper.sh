# Helper for ./3.0_run_iprscan_all.sh

IPR_DIR="$1"
outdir="$2"

echo "$IPR_DIR/interproscan.sh" \
    -i "$outdir/fasta_chunks/st_all.cds.pep.fasta.$LSB_JOBINDEX" \
    --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
    --goterms --pathways \
    --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
    --output-file-base "$outdir/interproscan_out/st_all.cds.pep.fasta.interproscan_out.$LSB_JOBINDEX"

"$IPR_DIR/interproscan.sh" \
    -i "$outdir/fasta_chunks/st_all.cds.pep.fasta.$LSB_JOBINDEX" \
    --applications TIGRFAM,SFLD,SUPERFAMILY,Hamap,Coils,ProSiteProfiles,SMART,CDD,PRINTS,PIRSF,ProSitePatterns,Pfam,ProDom,MobiDBLite \
    --goterms --pathways \
    --tempdir "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/tmp/" \
    --output-file-base "$outdir/interproscan_out/st_all.cds.pep.fasta.interproscan_out.$LSB_JOBINDEX"

