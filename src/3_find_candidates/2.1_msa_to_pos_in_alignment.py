'''Convert MSA to position in alignment of all STs
'''

import pandas as pd
from Bio import SeqIO, AlignIO
import pickle
import collections
import itertools

mafft = False
mauveAligner = True

if mafft:
    # MAFFT alignment
    aln = AlignIO.read(open('.output/mafft/refs_msa.fa'), 'fasta')

    # Create 1-indexed mappings of pos of non gaps to pos in alignment
    pos_to_pos_aln = {}
    for seq in list(aln):
        print(seq.id)
        pos_alns = (pos_aln+1 for pos_aln, base in enumerate(str(seq.seq)) if not base == '-')
        pos_to_pos_aln[seq.id] = dict(zip(range(1, len(seq)-seq.seq.count('-')+1), pos_alns))

    with open('.output/pos_to_pos_aln.mafft.pk', 'wb') as pk_file:
        pickle.dump(pos_to_pos_aln, pk_file)

if mauveAligner:

    def add_mauve_entry_to_mapping(current_aln_pos, mauve_pos_str, seq_ids, aligned_seqs=None, pos_to_pos_aln=None):
        '''Add mauve entry to mapping
        Example format: 153 292618 -2687311 294793
        #
        In a GappedAlignment, we need to consider the position of the gaps
        e.g.
            gagaggatgatgattataagggagtgtt----
            G----------GATTAT-CGAAAGTGTTAGAG
            G----------GATTAT-CGAAAGTGTTAGAG
            GAGGG---GATGATTATAAGGGAGTGGT----
        '''
        # Initialise mapping if one doesn't exist
        if not pos_to_pos_aln:
            pos_to_pos_aln = collections.defaultdict(list)
        # Get length of MUM
        anchor_locs  = [int(l) for l in mauve_pos_str.strip().split()]
        interval_len = anchor_locs.pop(0)
        assert(len(seq_ids) == len(anchor_locs))
        # In a GappedAlignment, we pass in aligned_seqs with gaps
        # If there are no aligned_seqs, generate dummy ungapped seqs for each seq_id
        if not aligned_seqs:
            aligned_seqs = ['X' * interval_len] * len(seq_ids)
        assert(len(seq_ids) == len(aligned_seqs))
        # Enumerate positions in the overall MSA
        msa_is = list(range(current_aln_pos, current_aln_pos+interval_len))
        for seq_id, anchor_loc, aligned_seq in zip(seq_ids, anchor_locs, aligned_seqs):
            if anchor_loc > 0:
                aligned_seq_i = itertools.count(anchor_loc, 1)
                aligned_seq_is = list(next(aligned_seq_i) if not base == '-' else None for base in aligned_seq)
            # Negative location indicates inversion:
            # Abs location counts down from previous anchor loc
            elif anchor_loc < 0:
                aligned_seq_i = itertools.count(-anchor_loc, 1)
                aligned_seq_is = list(reversed(list(next(aligned_seq_i) if not base == '-' else None for base in aligned_seq)))
            else:
                continue
            assert(len(aligned_seq_is) == len(msa_is))
            #  mapping = dict((aligned_seq_i, msa_i) for (aligned_seq_i, msa_i) in zip(aligned_seq_is, msa_is) if aligned_seq_i)
            mapping = ((aligned_seq_i, msa_i) for (aligned_seq_i, msa_i) in zip(aligned_seq_is, msa_is) if aligned_seq_i)
            #  assert(not (pos_to_pos_aln[seq_id].keys() & mapping.keys())), pos_to_pos_aln[seq_id].keys() & mapping.keys()
            #  pos_to_pos_aln[seq_id].update(mapping)
            pos_to_pos_aln[seq_id].extend(mapping)
        # Update the alignment pos
        current_aln_pos += interval_len
        return current_aln_pos, pos_to_pos_aln
    
    # Create 1-indexed mappings of pos of non gap bases in each sequence to pos in overall multiple seq alignment
    pos_to_pos_aln = None
    current_aln_pos = 1 
    seq_names = ['CC22_EMRSA15', 'TW20', 'CC30_MRSA252', 'CC8_USA300_FPR3757']
    #  seq_names = [0,1,2,3]
    seq_lens = []
    n_intervals, n_gapped_alignments, n_anchors = 0, 0, 0
    with open('.output/mauveAligner/refs.fa.mauve', 'r') as mauve_file:
        for line in mauve_file:
            # For gapped alignments, feed in the local CLUSTALW alignment
            if line.startswith('GappedAlignment'):
                current_aln_pos, pos_to_pos_aln = add_mauve_entry_to_mapping(current_aln_pos, next(mauve_file).strip(), seq_names, [next(mauve_file).strip() for _ in seq_names], pos_to_pos_aln)
                n_gapped_alignments += 1
                print('Processed: {} Intervals, {} GappedAlignment, {} MUM anchors ...'.format(n_intervals, n_gapped_alignments, n_anchors))
            # For anchors, just add a contiguous sequence of indices
            elif line[0] in '0123456789':
                current_aln_pos, pos_to_pos_aln = add_mauve_entry_to_mapping(current_aln_pos, line.strip(), seq_names, None, pos_to_pos_aln)
                n_anchors += 1
                print('Processed: {} Intervals, {} GappedAlignment, {} MUM anchors ...'.format(n_intervals, n_gapped_alignments, n_anchors))
            elif line.startswith('Interval'):
                n_intervals += 1
                print('Processed: {} Intervals, {} GappedAlignment, {} MUM anchors ...'.format(n_intervals, n_gapped_alignments, n_anchors))
            # Read total length for each seq from headers
            elif 'Sequence' in line and 'Length' in line:
                seq_lens.append(int(line.strip().split()[1]))

    # Test case
    #
    #  current_aln_pos=34352
    #  line='121 34558 -2232571 34517 34569'
    #  seqs=[
    #  '------------tgtagttaaaaaatctgattgggataaaggtgatctatataaaactttagtccatgataagttacccaagcagt----taaa--agtgcata----taaaagaa-----',
    #  'GCATAAAATCTTAATTGTTAATTCAACTGGTATTGTTGAAAATAATGGTTTAAATGGTTTAATTGCTATGAAATTAGCTAGAGAATATGGTAGACCAGTGCTAATGGTTAAAAGAATAGGG',
    #  '------------TGTAGTTAAAAAATCTGATTGGGATAAAGGTGATCTATATAAAACTTTAGTCCATGATAAGTTACCCAAGCAGT----TAAA--AGTGCATA----TAAAAGAA-----',
    #  '------------TGTAGTTAAAAAATCTGATTGGGATAAAGGTGATCTATATAAAACTTTAGTCCATGATAAGTTACCCAAGCAGT----TAAA--AGTGCATA----TAAAAGAA-----']
    #  current_aln_pos, pos_to_pos_aln = add_mauve_entry_to_mapping(current_aln_pos, line.strip(), seq_names, None, pos_to_pos_aln)

    # Convert lists to mappings
    seq_lens = dict(zip(seq_names, seq_lens))
    for seq_name in pos_to_pos_aln.keys():
        mapping_list = pos_to_pos_aln[seq_name]
        mapping = dict(mapping_list)
        # Check there are no doubly-mapped positions
        assert(len(mapping_list) == len(mapping))
        # Check every base is mapped
        assert(len(mapping_list) == seq_lens[seq_name])
        #
        pos_to_pos_aln[seq_name] = mapping

    with open('.output/pos_to_pos_aln.mauveAligner.pk', 'wb') as pk_file:
        pickle.dump(pos_to_pos_aln, pk_file)

