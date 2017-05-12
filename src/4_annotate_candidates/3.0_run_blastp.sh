
# blastp against nr protein db

input_pep_fa=".output/interproscan/st_all.cds.pep.unannotated.fasta"
blastdb_nr='/data/blastdb/Supported/NR/nr'
outdir=".output/blastp"

mkdir -p "$outdir"

bsub \
    -G team81 -q long \
    -n 4 -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
    -J "blastp_unannotated" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        blastp \
            -query "$input_pep_fa" \
            -db "$blastdb_nr" \
            -out "$outdir/st_all.cds.pep.unannotated.blastp.outfmt6.tsv" \
            -outfmt 6 \
            -max_target_seqs 1 \
            -num_threads 4

