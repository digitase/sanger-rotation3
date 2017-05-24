'''Add sRNAs to .embl annotation for st22
'''

def add_ncRNAs(merged):
    #  Changed misc_feature Feature Key to ncRNA
    # Adds locus tag qualifiers
    with open('./Mels_currated_EMRSA15.tab', 'r') as infile:
        for line in infile:
            # position line for entry
            if line.startswith('FT   misc_feature'):
                line = line.replace('misc_feature', 'ncRNA       ')
                merged.append(line)
                # get next line with locustag
                #  e.g.
                #  FT                   /note="SAUSA300s001	T-box riboswitch"
                line = next(infile)
                locus_tag = line.split('"')[1].split()[0]
                product = ' '.join(line.split('"')[1].split()[1:])
                merged.append('FT                   /locus_tag="{}"\n'.format(locus_tag))
                merged.append('FT                   /product="{}"\n'.format(product))
                merged.append(line)
            else:
                merged.append(line)

#  Fixed the header line for the main annotation
merged = []
with open('./EMRSA15.art', 'r') as infile:
    for line in infile:
        if line.startswith('ID'):
            merged.append('ID   EM0000001; SV 1; circular; genomic DNA; STD; PRO; 2832299 BP.\n')
        #  Appended sRNAs at the bottom of the main annotation
        elif line.startswith('FT                   /systematic_id="'):
            merged.append(line)
            merged.append(line.replace('systematic_id', 'locus_tag'))
        #  FT                   /note="tRNA Leu anticodon GAG, Cove score 40.45"
        elif line.startswith('FT                   /note="tRNA '):
            aa, _, codon = line.split('"')[1].split()[1:4]
            merged.append('FT                   /locus_tag="SAEMRSA15tRNA{}{}"\n'.format(aa, codon))
            merged.append(line)
        #  FT                   /gene="tRNA_47"
        #  FT   tRNA            complement(1945407..1945478)
        #  FT                   /domain="Rfam:RF00005;tRNA;Score=38.40;positions 1 to 59"
        #  FT                   /gene="tRNA_50"
        elif line.startswith('FT                   /gene="tRNA_'):
            merged.append(line)
            gene = line.split('"')[1]
            merged.append('FT                   /locus_tag="SAEMRSA15tRNA{}"\n'.format(gene.replace('_', '')))
            merged.append('FT                   /product="{}"\n'.format(gene))
        elif line.startswith('FT   rRNA'):
            merged.append(line)
            #  FT   rRNA            complement(1951221..1952775)
            #  FT                   /label=23S
            loc_str = line.strip().split()[-1]
            start, end = loc_str.replace('complement', '').replace('(', '').replace(')', '').split('..')
            line = next(infile)
            label = line.strip().replace('FT                   /label=', '')
            if 'complement' in loc_str:
                merged.append('FT                   /locus_tag="SAEMRSA15rRNA{}{}{}"\n'.format(label, start, end))
            else:
                merged.append('FT                   /locus_tag="SAEMRSA15rRNA{}{}{}"\n'.format(label, end, start))
            merged.append('FT                   /product="{} rRNA st22 loc: {}"\n'.format(label, loc_str))
            merged.append('FT                   /gene="{}"\n'.format(label))
            merged.append(line)
        elif line.startswith('SQ'):
            add_ncRNAs(merged)
            merged.append(line)
        else:
            merged.append(line)
with open('./EMRSA15_with_Mels_currated.modified.embl', 'w') as outfile:
    for line in merged:
        outfile.write(line)

