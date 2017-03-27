#!/bin/bash
#

outdir="/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/1_reads_to_trees/1_mapping/"
datadir="/nfs/users/nfs_b/bb9/workspace/rotation3/data/"

#
# Link in fastqs
#

while read seq_tag; do
    pathfind -t lane -id "$seq_tag" -f fastq -symlink "$outdir/fastq/"
done < "./seq_tags_st8.txt"

#
# Mapping
#

# Copy in dna and embl file
# Write permissions are required on these files
cp "$datadir/st8/CC8_USA300_FPR3757.dna" "$outdir"
cp "$datadir/st8/CC8_USA300_FPR3757.embl" "$outdir"


# Map with SMALT
# .dna: reference seq
# .embl: annotations
# -v: SMALT aligner version
#    
# TODO try using the v1.2 script version, usng BWA over SMALT. Use -y dirty flag to track errors
cd "$outdir/bam"
export PATH="/usr/bin/:$PATH"
python2 ~sh16/scripts/multiple_mappings_to_bam.py \
     -r "$outdir/CC8_USA300_FPR3757.dna" \
     -v 0.7.4 \
     -E -G \
     -o st8_bsac \
     -F 2 \
     -g -t -a \
     -M 2 \
     -k -I \
     "$outdir/fastq/7414_7#18_1.fastq.gz" "$outdir/fastq/7414_7#18_2.fastq.gz"
cd -

#
# Contamination check TODO
#

 #All bacteria samples are automatically assembled, using Velvet, Sspace and
 #GapFiller as soon as they come off the sequencer. We have seen that e.g. an
 #unusal total length or an unusually high number of contigs can be a hint for
 #contamination. To get this information the best way is to create a spreadsheet
 #of statistics on each assembly.  
 #assemblyfind -t lane -i 8113_4 -stats

 #...look at heterozygous SNPs
 #This approach needs the mapping to be done first. The easiest way then is to
 #link or copy the *.bcf files into a new directory, ignoring the *_variant.bcf
 #files. Then run the reportlabtest script for every lane using the following
 #command: 

 #bsub -o out -e err -R "select[mem > 8000] rusage[mem=8000]" -M 8000000 ~sh16/scripts/reportlabtest.py -n -o <output.pdf> -v H <*.bcf> <reference.embl>

