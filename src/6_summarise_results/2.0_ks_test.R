# Determine whether locations of snps and hp snps relative to loci are drawn from the same distribution

library(Matching)
library(data.table)
library(afex)

hp_files = c(
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST22_BSAC_Pfizer/st22_homoplasies.csv',
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST239_global_Singapore_Pfizer/st239_homoplasies.csv',
    # '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST30_BSAC_Pfizer/st30_homoplasies.csv',
    '/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/S.aureus_ST8_BSAC_Pfizer_revised/st8_homoplasies.csv'
)
names(hp_files) <- c(
    'st22', 
    'st239', 
    # 'st30', 
    'st8'
)

# read in hp locations
hps <- list()
i <- 1
for (hp_file in hp_files) {
    hp <- fread(hp_file)
    hp$dist_to_locus <- unlist(lapply(hp$intergenic, function(x) strtoi(gsub('\\(|,|\\)', '', as.character(x)))))
    hps[[i]] <- hp
    i <- i + 1
}
hps <- rbindlist(hps)

sink('.output/ks_test.results.txt', split=T)
print('ks.boot package:Matching Bootstrap Kolmogorov-Smirnov')
ks.boot(
    hps[(n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0), 'dist_to_locus'], 
    hps[!(n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0), 'dist_to_locus'], 
    nboots=1000, alternative="two.sided", print.level=1
)
# Two-sample Kolmogorov-Smirnov test

# data:  Tr and Co
# D = 0.10313, p-value = 5.607e-14
# alternative hypothesis: two-sided

hps[is.na(dist_to_locus)]
print('compare.2.vectors package:afex Compare two vectors using various tests.')
compare.2.vectors(
    hps[(n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0), 'dist_to_locus'], 
    hps[!(n_homoplasic_acctran > 0 | n_homoplasic_deltran > 0), 'dist_to_locus'],
    na.rm=T
)
# $parametric
   # test test.statistic test.value  test.df            p
# 1     t              t   4.464178 28285.00 8.068698e-06
# 2 Welch              t   2.604557  1603.28 9.284332e-03

# $nonparametric
             # test test.statistic    test.value test.df            p
# 1 stats::Wilcoxon              W  2.000169e+07      NA 0.0006644721
# 2     permutation              Z  4.462685e+00      NA 0.0000100000
# 3  coin::Wilcoxon              Z -3.403836e+00      NA 0.0007300000
# 4          median              Z  8.676931e+00      NA 0.0000000000

sink()

