# Using Homoplasy to Study Recent Selection in *Staphylococcus aureus*

### Informatics rotation project for 3rd rotation (2016-17)

Repo for project-associated code.

### Supervisors: Julian Parkhill and Simon Harris

## Project description

Within closely related populations of bacterial species with low recombination rates, homoplasy in phylogenetic reconstructions is rare and indicative of strong selective pressures acting at the homoplasic sites. In a previous analysis of a small number of genomes from one sequence type of *Staphylococcus aureus*, it was noted that most of the small number of homoplasic mutations identified were associated with antimicrobial resistance, but others were unexplained. In this project we will examine a large collection of *Staphylococcus aureus* genomes sequenced previously in collaboration with the British Society for Antimicrobial Chemotherapy to search for sites under selection using homoplasy analysis. The BSAC collection included over 1000 isolates in many sequence types. By looking for similarities in the sites and genes exhibiting homoplasy in these different lineages we aim to identify novel loci which may be important for virulence, disease severity or transmissibility. Once candidate loci have been identified, pathway analysis will be used to give insight into their function, and where possible, analysis of the possible effects of any identified mutations on protein structure will be examined.

## Notes

`./data/`

- "There is a separate directory for each ST/lineage that we discussed. Within
each directory you will find the reference strain used for mapping to generate
the alignment (e.g. For ST22: CC22_EMRSA15.fasta), the whole genome alignment
(e.g. S.aureus_ST22_BSAC_Pfizer.aln), core genome alignment (MGE blocks masked,
e.g. S.aureus_ST22_BSAC_Pfizer_no_mge.aln) and the core genome snp alignment
(e.g. S.aureus_ST22_BSAC_Pfizer_no_mge_snp.aln). I’ve generated a tree for each
ST (except ST8, still running) for you to use (RAxML_bestTree*). For each
lineage the reference strain should be included in the alignment and tree
files, this is needed if you want to know the context of a SNP position (e.g.
Intergenic, synonymous, non-synonymous)."
