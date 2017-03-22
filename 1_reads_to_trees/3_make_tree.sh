
#3) Generating a phylogenetic tree TODO

#3a) Using RAxML to generate a maximum likelihood tree

#~sh16/scripts/run_RAxML.py -a <alignment file>  -s -o <new aligment file name>

#or extract SNPs into an alignment using 

#snp_sites -o <output filename> inputfile name

#and then use this snps.aln in RAxML without the -s flag

#There are several files produced with the prefix ‘RAxML_’, open the bipartition file with FigTree to look at the tree.

#3b) Using FastTree to generate a neighbour joining tree

#Extraction of SNPs 

#~sh16/scripts/summarise_snps.py -i <alignment file> -r <name of your reference strain> -o <outfile prefix> -t -a

#Generating tree 

#~tc7/scripts/pub/run_FastTree.pl -o <output.tree> -i <alignmentfile.aln> -m 30 -f -n -j
#running FastTree on farm3: FastTree -nt <input.aln> > <output.tee>

