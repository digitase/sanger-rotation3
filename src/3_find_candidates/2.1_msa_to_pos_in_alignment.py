'''Convert MSA to position in alignment of all STs
'''

import pandas as pd
from Bio import SeqIO, AlignIO
import pickle

aln = AlignIO.read(open('.output/mafft/refs_msa.fa'), 'fasta')

# Create 1-indexed mappings of pos of non gaps to pos in alignment
pos_to_pos_aln = {}
for seq in list(aln):
    print(seq.id)
    pos_alns = (pos_aln+1 for pos_aln, base in enumerate(str(seq.seq)) if not base == '-')
    pos_to_pos_aln[seq.id] = dict(zip(range(1, len(seq)-seq.seq.count('-')+1), pos_alns))

with open('.output/pos_to_pos_aln.pk', 'wb') as pk_file:
    pickle.dump(pos_to_pos_aln, pk_file)

