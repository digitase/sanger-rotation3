# get homoplasy indices for trees
library(ape)
library(phangorn)

alns <- c(
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/S.aureus_ST22_BSAC_Pfizer_no_mge_snp.aln",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/S.aureus_ST239_global_Singapore_Pfizer_no_mge_snp.aln",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/revised_lane_list/S.aureus_ST8_BSAC_Pfizer_revised_no_mge_snp.aln"
)

trees <- c(
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st22/RAxML_bestTree.ml_S.aureus_ST22_BSAC_Pfizer",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st239/RAxML_bestTree.ml_S.aureus_ST239_global_Singapore_Pfizer",
    "/nfs/users/nfs_b/bb9/workspace/rotation3/data/st8/revised_lane_list/RAxML_bestTree.ml_S.aureus_ST8_BSAC_Pfizer_revised"
)

sink('.output/hp_indices.results.txt', split=T)
for (i in 1:length(alns)) {
    print(alns[i])
    aln <- read.phyDat(alns[i], format="fasta")
    print(trees[i])
    tree <- read.tree(trees[i], comment.char='~')
    print('Consistency index')
    ci <- CI(tree, aln, cost=NULL)
    print(ci)
    print('Retention index')
    ri <- RI(tree, aln, cost=NULL)
    print(ri)
    print('Rescaled consistency index = ci * ri')
    print(ci*ri)
    print('Homoplasy index = 1-ci')
    print(1-ci)
}
sink()

