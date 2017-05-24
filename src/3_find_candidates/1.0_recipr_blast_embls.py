'''Get homology between features across strains
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import itertools
import sys
import os

out_base_dir = '/nfs/users/nfs_b/bb9/workspace/rotation3/src/3_find_candidates/.output/'
#
embl_files = {
    #  "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/misc/2017-05-23_staph_annot_from_matt/EMRSA15_with_Mels_currated.modified.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}

# For each embl file, write all CDS sequences, headed by locus_tags
cds_fastas = {}
for short_prefix, embl_file in embl_files.items():
    print('Writing CDSs...: ' + short_prefix)
    embl_record = next(SeqIO.parse(embl_files[short_prefix], 'embl'))
    #
    def feature_to_SeqRecord(feature):
        try:
            seq = feature.extract(embl_record.seq).translate(table=11, cds=True, to_stop=False)
        except CodonTable.TranslationError:
            seq = feature.extract(embl_record.seq).translate(table=11, cds=False, to_stop=False)
        rec = SeqRecord(
            seq,
            id='{locus_tag}_{begin}_{end}_{strand}'.format(
                locus_tag=feature.qualifiers['locus_tag'][0],
                begin=feature.location.start.position,
                end=feature.location.end.position,
                strand=feature.strand
            ),
            name=feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else '',
            description=feature.qualifiers['product'][0] if 'product' in feature.qualifiers else ''
        )
        return rec
    #
    os.makedirs(os.path.join(out_base_dir, 'blast', 'fasta'), exist_ok=True)
    out_fasta_file = os.path.join(out_base_dir, 'blast', 'fasta', '{short_prefix}.cds.pep.fasta'.format(short_prefix=short_prefix))
    SeqIO.write((feature_to_SeqRecord(f) for f in embl_record.features if f.type == 'CDS'), handle=out_fasta_file, format='fasta')
    #
    cds_fastas[short_prefix] = out_fasta_file

# Run reciprocal blasts

os.makedirs(os.path.join(out_base_dir, 'blast', 'blast_rbh'), exist_ok=True)
#
# Run for each pair of strains
for sp1, sp2 in itertools.combinations(cds_fastas, 2):
    fa1, fa2 = cds_fastas[sp1], cds_fastas[sp2]
    #
    bsub_cmd = (r'bsub -G team81 -q normal -n4 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" -M 1000 -J "blast_rbh.{sp1}.{sp2}" -o "{olog}" -e "{elog}"'.format(
        sp1=sp1, 
        sp2=sp2,
        olog=os.path.join(out_base_dir, 'blast', 'blast_rbh', 'bsub_o.blast_rbh.{sp1}_{sp2}.log'.format(sp1=sp1, sp2=sp2)),
        elog=os.path.join(out_base_dir, 'blast', 'blast_rbh', 'bsub_e.blast_rbh.{sp1}_{sp2}.log'.format(sp1=sp1, sp2=sp2)))
    , )
    #
    #  Options:
      #  -h, --help            show this help message and exit
      #  -a DBTYPE, --alphabet=DBTYPE
                            #  Alphabet type (nucl or prot; required)
      #  -t TASK, --task=TASK  BLAST task (e.g. blastp, blastn, megablast; required)
      #  -i MIN_IDENTITY, --identity=MIN_IDENTITY
                            #  Minimum percentage identity (optional, default 70)
      #  -c MIN_COVERAGE, --coverage=MIN_COVERAGE
                            #  Minimum HSP coverage (optional, default 50)
      #  --nr                  Preprocess FASTA files to collapse identitical entries
                            #  (make sequences non-redundant)
      #  -o FILE, --output=FILE
                            #  Output filename (required)
      #  --threads=THREADS     Number of threads when running BLAST. Defaults to the
                            #  $GALAXY_SLOTS environment variable if set, or 1.
    cmd = ' '.join(
        bsub_cmd + (
            r'python {blast_rbh_script}',
            r'-a prot -t blastp',
            r'--identity=70 --coverage=50',
            r'-o {out} --threads=4',
            r'{fa1} {fa2}'
        )
    ).format(
        blast_rbh_script='/nfs/users/nfs_b/bb9/packages/galaxy_blast/tools/blast_rbh/blast_rbh.py',
        out=os.path.join(out_base_dir, 'blast', 'blast_rbh', 'blast_rbh.{sp1}_{sp2}.out'.format(sp1=sp1, sp2=sp2)),
        fa1=fa1, fa2=fa2
    )
    print('Executing...:')
    print(cmd)
    #
    os.system(cmd)

