'''Get homology between features across strains
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import sys


embl_files = {
    "st22": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.embl",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.embl",
    "st30": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.embl",
    "st8": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.embl",
}
#
ref_files = {
    "st22" : "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/CC22_EMRSA15.fasta",
    "st239": "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/CC8_TW20.dna",
    "st30" : "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st30/CC30_MRSA252.fasta",
    "st8"  : "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/CC8_USA300_FPR3757.dna",
}


# For each embl file, write all CDS sequences, headed by locus_tags
for short_prefix, embl_file in embl_files.items():

    short_prefix = 'st22'

    embl_record = next(SeqIO.parse(embl_files[short_prefix], 'embl'))

    out_fasta_records = []

    def feature_to_SeqRecord(feature):
        print(feature)
        try:
            seq = feature.extract(embl_record.seq).translate(table=11, cds=True, to_stop=False)
        except CodonTable.TranslationError:
            seq = feature.extract(embl_record.seq).translate(table=11, cds=False, to_stop=False)
        rec = SeqRecord(
            seq,
            id=feature.qualifiers['locus_tag'][0],
            name=feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else '',
            description=feature.qualifiers['product'][0] if 'product' in feature.qualifiers else ''
        )
        return rec

    SeqIO.write((feature_to_SeqRecord(f) for f in embl_record.features if f.type == 'CDS'), handle=sys.stdout, format='fasta')

