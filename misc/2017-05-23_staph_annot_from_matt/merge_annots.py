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
        elif line.startswith('SQ'):
            add_ncRNAs(merged)
            merged.append(line)
        else:
            merged.append(line)
with open('./EMRSA15_with_Mels_currated.modified.embl', 'w') as outfile:
    for line in merged:
        outfile.write(line)

