# GO Enrichment of loci near homoplasic sites

library(AnnotationDbi, lib.loc='/nfs/users/nfs_b/bb9/R/x86_64-pc-linux-gnu-library/3.3')
library(topGO)
library(Rgraphviz)
library(data.table)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)

# Read in GO mappings
iprscan_files <- list.files( 
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_annotate_candidates/interproscan_all_array/interproscan_out/', 
    pattern='*.tsv', full.names=T
)
mappings <- list()
for (f in iprscan_files) {
    print(f)
    annots <- fread(f, sep='\t', fill=T)
    new_mappings <- annots[V14 != "", .(categories=paste(unique(strsplit(paste(V14, collapse='|'), '|', fixed=T)[[1]]), collapse='|')), by=V1]
    new_mappings = setNames(lapply(new_mappings$categories, function(x) strsplit(x, '\\|')[[1]]), new_mappings$V1) 
    mappings <- c(mappings, new_mappings)
}

# Gene universe is the set from which the genes of interest are taken:
# Coding sequences from the genome
geneUniverse <- sapply(
    names(readDNAStringSet('/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/4_annotate_candidates/interproscan_all_array/st_all.cds.pep.fasta')),
    function(x) strsplit(x, ' ')[[1]][1]
)

# Read in locations of hp snps
# hp_files = c(
    # '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/st22_homoplasies.csv',
    # '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/st239_homoplasies.csv',
    # '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies.csv'
# )
# without known amr sites
hp_files = c(
             '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st22_homoplasies.known_amr_filtered.csv',
             '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st239_homoplasies.known_amr_filtered.csv',
             '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/6_summarise_results/st8_homoplasies.known_amr_filtered.csv')

names(hp_files) <- c(
    'st22', 
    'st239', 
    'st8'
)

genesOfInterest <- c()
for (hp_file in hp_files) {
    print(hp_file)
    hp <- fread(hp_file)
    hp$dist_to_locus <- unlist(lapply(hp$intergenic, function(x) strtoi(gsub('\\(|,|\\)', '', as.character(x)))))
    hp[, geneNames := gsub("\\(|,|\\)|'", '', locus_tag)]
    # Define set of interest to be those loci that are just downstream of hp sites
    genesOfInterest.new <- unique(hp[dist_to_locus < 0 & dist_to_locus >= -200, gsub("\\(|,|\\)|'", '', locus_tag)])
    # genesOfInterest.new <- unique(gsub("\\(|,|\\)|'", '', hp$locus_tag))
    genesOfInterest <- c(
        genesOfInterest,
        unlist(strsplit(genesOfInterest.new, split=' '), recursive=F)
    )
}

outdir <- '.output/topGO/'
dir.create(outdir)
setwd(outdir)

# Enrich
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
#
for (onto in c('BP', 'MF', 'CC')) {
    GOdata <- new("topGOdata", ontology=onto, allGenes=geneList, annot=annFUN.gene2GO, gene2GO=mappings, nodeSize=5)
    classic.fisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
    weight01.fisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
    # hist(score(weight01.fisher), 50, xlab = "p-values")
    allRes <- GenTable(GOdata, classic=classic.fisher, weight=weight01.fisher,
        orderBy="weight", ranksOf="classic", topNodes=50, numChar=999)
    write.csv(allRes, sprintf('%s.genTable.csv', onto))
    printGraph(GOdata, classic.fisher, firstSigNodes=5, fn.prefix=onto, useInfo='all', pdfSW=T)
    printGraph(GOdata, weight01.fisher, firstSigNodes=5, fn.prefix=onto, useInfo='all', pdfSW=T)
}

