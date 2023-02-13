import os, glob, re
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def find_outer(tree_file, method):
    if method == 'concatenation':
        with open(tree_file) as tree_handle:
            tree_str_line = tree_handle.readline()
            stripped_tree_str = tree_str_line.strip().strip(';').strip('(').strip(')')
            in_parentheses_str = re.search(r'\(.+\)\d+\.?\d*', stripped_tree_str).group()
            print(in_parentheses_str)
            in_parentheses_list = re.findall(r'n\d', in_parentheses_str)
            bs = float((in_parentheses_str.split(')')[1]))
            if 'n0' in in_parentheses_list:
                in_parentheses_list.remove('n0')
                outer = in_parentheses_list[0]
            else:
                outer = set(['n1', 'n2', 'n3']).difference(set(in_parentheses_list)).pop()
            # print(tree_str_line, outer)
    if method == 'alstral':
        with open(tree_file) as tree_handle:
            tree_str_line = tree_handle.readline()
            if not tree_str_line.strip():
                return None, None
            in_parentheses_str = re.search(r'\(n\d\,n\d\)', tree_str_line).group()
            print(in_parentheses_str)
            in_parentheses_list = re.findall(r'n\d', in_parentheses_str)
            outer = set(['n1', 'n2', 'n3']).difference(set(in_parentheses_list)).pop()
            bs = float(re.search(r'\)\d+\.?\d*', tree_str_line).group().strip(')'))

    return outer, bs


def plot_species_tree(result_root, ax, method):
    color_list = ["#75b84f", "cornflowerblue", "#ffad0f",]

    outer_list = []
    if method == 'concatenation':
        treefile_list = glob.glob(os.path.join(result_root, '*.fasta.treefile'))
    elif method == 'alstral':
        treefile_list = glob.glob(os.path.join(result_root, '*.species.treefile'))
    else:
        raise Exception('Unknown phylogeny reconstruction method')

    for tree_file in treefile_list:
        outer, bs = find_outer(tree_file, method)
        if outer:
            outer_list.append(outer)
    counter = Counter(outer_list)

    y = []
    for i in range(3):
        species = 'n'+str(i+1)
        if species in counter:
            y.append(counter[species]/len(outer_list))
        else:
            y.append(0)
    bar_list = ax.bar(range(3), y, edgecolor="", width=0.8, ecolor="#202020")
    for bar, color in zip(bar_list, color_list):
        bar.set(color=color)
    ax.set_xticks(range(3))
    ax.set_xticklabels(['(1)','(2)','(3)'])
    

def plot_all():
    result_root_list = [\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80',\
    ]

    phylogeny_fig_file = 'phylogeny.pdf'

    with PdfPages(phylogeny_fig_file) as phylogeny_fig_handle:   
        fig, axes = plt.subplots(2, len(result_root_list), figsize=(7,5), sharey=True)
        for i, result_root in enumerate(result_root_list):
            concatenation_ax, alstral_ax = axes[:,i]
            plot_species_tree(result_root, concatenation_ax, 'concatenation')
            plot_species_tree(result_root, alstral_ax, 'alstral')
        phylogeny_fig_handle.savefig()


if __name__ == '__main__':
    plot_all()