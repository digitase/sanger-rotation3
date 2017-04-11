
import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo
import pandas as pd
import itertools
import collections
import os
import sys
import copy

prefixes = {
    "st22": "S.aureus_ST22_BSAC_Pfizer",
    "st239": "S.aureus_ST239_global_Singapore_Pfizer",
    "st30": "S.aureus_ST30_BSAC_Pfizer",
    "st8": "S.aureus_ST8_BSAC_Pfizer_revised"
}
tree_file_template = "/lustre/scratch118/infgen/team81/bb9/2_homoplasy/reconstruct_snps_on_tree/{transformation}/{prefix}/{prefix}_{transformation}_steps.tre"
changes_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_homoplasies/{prefix}/{prefix}_homoplasies_on_tree.{transformation}.csv"
snps_file_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/summarise_snps/{prefix}/{prefix}.out"
out_pdf_template = "/nfs/users/nfs_b/bb9/workspace/rotation3/lustre/2_homoplasy/plots/{prefix}/{prefix}.{transformation}.loc_{loc}.pdf"

short_prefix = 'st239'
colors = {'A': 'green', 'C': 'blue', 'G': 'black', 'T': 'red'}
transformations = ['acctran', 'deltran']
locs = [7255]

print('Parsing snps file...')
prefix = prefixes[short_prefix]
snps = pd.read_table(snps_file_template.format(prefix=prefix))

for loc in locs:
    for transformation in transformations:
        print('Analysing {}: prefix {}, loc {}, transformation {} ...'.format(short_prefix, prefix, loc, transformation))
        changes = pd.read_csv(changes_file_template.format(prefix=prefix, transformation=transformation))
        tree = next(Phylo.parse(tree_file_template.format(prefix=prefix, transformation=transformation), 'newick'))
        
        # Color code bases
        changes_loc = changes[changes['loc'] == loc]

        # Color code homoplasic changes
        print('Coloring branches...')
        node_colors = dict()
        for i, change in changes_loc.iterrows():
            snp0, snp1 = change.ix['SNP'].replace('"', '').split('->')
            node0, node1 = change.ix['node'].replace('"', '').split('->')
            node_colors[node0] = colors[snp0]
            node_colors[node1] = colors[snp1]

        # Color tree branches with homoplasic changes
        tree_copy = copy.deepcopy(tree)
        for clade in tree_copy.find_clades():
            if not clade.name:
                clade.name = str(clade.confidence)
            if clade.name in node_colors:
                clade.color = node_colors[clade.name]
            else:
                clade.color = 'gray'
        tree_copy.root_at_midpoint()
        tree_copy.ladderize()

        # Build dictionary for tip label colors
        print('Getting tip colors...')
        snps.columns = map(str.strip, snps.columns)
        # Get reference prefix and base
        snps_refname = snps.columns[1].replace('Position_in_', '')
        snps_refbase = snps[snps.iloc[:, 1] == loc]['Ref_base'].values[0]
        # Get other prefixes and bases, then color
        snps_taxa = [snps_refname] + list(snps.columns[10:])
        snps_base = [snps_refbase] + list(snps.loc[snps.iloc[:, 1] == loc, snps.columns[10:]].values[0]) 
        snps_base = [snps_refbase if base == '.' else base for base in snps_base]
        snps_color = [colors[base] if base in colors else 'grey' for base in snps_base]

        # Plot tree
        plt.rcParams['font.size'] = 6
        plt.rcParams['lines.linewidth'] = 2
        plt.rcParams['figure.figsize'] = [24, 48]
        #
        print('Plotting...')
        #  fig = plt.figure()
        Phylo.draw(tree_copy, do_show=False, show_confidence=False, label_colors=dict(zip(snps_taxa, snps_color)))
        #
        plt_file = out_pdf_template.format(prefix=prefix, loc=loc, transformation=transformation)
        os.makedirs(os.path.dirname(plt_file), exist_ok=True)
        plt.savefig(plt_file)

