
# blastp against nr protein db

# input_pep_fa=".output/interproscan/st_all.cds.pep.unannotated.test.fasta"
input_pep_fa=".output/interproscan/st_all.cds.pep.unannotated.fasta"
blastdb_nr='/data/blastdb/Supported/NR/nr'
outdir=".output/blastp"

mkdir -p "$outdir"

outfmt=11
nthreads=4
outblastpfile="$outdir/st_all.cds.pep.unannotated.blastp.outfmt$outfmt"

bsub \
    -G team81 -q long \
    -n "$nthreads" -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" -M 24000 \
    -J "blastp_unannotated" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        blastp \
            -query "$input_pep_fa" \
            -db "$blastdb_nr" \
            -out "$outblastpfile" \
            -outfmt "$outfmt" \
            -num_threads "$nthreads" \
            -max_target_seqs 3

bsub \
    -G team81 -q normal \
    -R "select[mem>16000] rusage[mem=16000]" -M 16000 \
    -J "blast_formatter_unannotated" \
    -o "$outdir/jobid_%J.bsub_o.log" \
    -e "$outdir/jobid_%J.bsub_e.log" \
        blast_formatter \
            -archive "$outblastpfile" \
            -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 1 \
            -out "$outdir/st_all.cds.pep.unannotated.blastp.outfmt6cust"

