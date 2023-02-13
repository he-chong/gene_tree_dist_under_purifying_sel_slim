import os, glob, time, re
import subprocess
from collections import defaultdict


def concatenation(repeat_dir_list, out_dir, total):
    cat_fasta_file = os.path.join(out_dir, 'repeat_{}.fasta'.format(total))
    if not os.path.exists(cat_fasta_file): 
        out_aln_dict = defaultdict(str)
        for repeat_dir in repeat_dir_list:
            for tree_file in glob.glob(os.path.join(repeat_dir, '*.fasta')):        
                with open(tree_file) as tree_handle:
                    for line in tree_handle:
                        if '>' in line:
                            species = line.strip('>').strip('\n')
                        else:
                            seq = line.strip('\n')
                            out_aln_dict[species] += seq

        with open(cat_fasta_file, 'w') as cat_fasta_handle:
            for species, seq in sorted(out_aln_dict.items(), key=lambda i:int(i[0].strip('n'))):
                cat_fasta_handle.write('>')
                cat_fasta_handle.write(species)
                cat_fasta_handle.write('\n')
                cat_fasta_handle.write(seq)
                cat_fasta_handle.write('\n\n')

    cat_tree_file = cat_fasta_file + '.treefile'
    if not os.path.exists(cat_tree_file):
        p = subprocess.Popen(['iqtree', '-s', cat_fasta_file, '-m', 'JC', '-o', 'n0', '-bb', '1000', '-nt', '15'], stdout=subprocess.PIPE)
        p.communicate()


def read_tree_file(tree_file):
    with open(tree_file) as tree_handle:
        tree_str_line = tree_handle.readline()
        stripped_tree_str = tree_str_line.strip().strip(';').strip('(').strip(')')
        in_parentheses_str = re.search(r'\(.+\)\d+', stripped_tree_str).group()
        bs = int(in_parentheses_str.split(')')[1])
    return tree_str_line, bs


def alstra(iqtree_repeat_dir_list, out_dir, total, bs_cutoff=50):
    gene_tree_file = os.path.join(out_dir, 'repeat_{}.genes.treefile'.format(total))
    if not os.path.exists(gene_tree_file):
        with open(gene_tree_file, 'w') as gene_tree_handle:
            for iqtree_repeat_dir in iqtree_repeat_dir_list:
                for single_tree_file in sorted(glob.glob(os.path.join(iqtree_repeat_dir, '*.fasta.treefile')), key=lambda i:int(os.path.basename(i).split('.')[0].split('_')[-1])):
                    tree_str_line, bs = read_tree_file(single_tree_file)
                    if bs >= bs_cutoff:
                        gene_tree_handle.write(tree_str_line)

    species_tree_file = os.path.join(out_dir, 'repeat_{}.species.treefile'.format(total))
    if not os.path.exists(species_tree_file):
        p = subprocess.Popen(['astral', '-i', gene_tree_file, '-o', species_tree_file], stdout=subprocess.PIPE)
        p.communicate()


def infer_a_scenario(result_dir, iqtree_dir, out_dir, repeat_num):
    result_dir = os.path.abspath(result_dir)
    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    total_repeat_dir_list = []
    total_iqtree_repeat_dir_list = []
    for subset_base in os.listdir(result_dir):
        for repeat_base in os.listdir(os.path.join(result_dir, subset_base, 'sequence_fasta')):
            repeat_dir = os.path.join(result_dir, subset_base, 'sequence_fasta', repeat_base)
            iqtree_repeat_dir = os.path.join(iqtree_dir, subset_base, repeat_base)
            total_repeat_dir_list.append(repeat_dir)
            total_iqtree_repeat_dir_list.append(iqtree_repeat_dir)

    repeat_size = int(len(total_repeat_dir_list)/repeat_num)
    repeat_dir_list_reshape = []
    iqtree_repeat_dir_list_reshape = []
    for i in range(repeat_num):
        repeat_dir_list_each_repeat = []
        iqtree_repeat_dir_list_each_repeat = []
        for j in range(repeat_size):
            repeat_dir = total_repeat_dir_list.pop()
            iqtree_repeat_dir = total_iqtree_repeat_dir_list.pop()
            repeat_dir_list_each_repeat.append(repeat_dir)
            iqtree_repeat_dir_list_each_repeat.append(iqtree_repeat_dir)
        repeat_dir_list_reshape.append(repeat_dir_list_each_repeat)
        iqtree_repeat_dir_list_reshape.append(iqtree_repeat_dir_list_each_repeat)

    total = 0
    for repeat_dir_list_each_repeat, iqtree_repeat_dir_list_each_repeat in zip(repeat_dir_list_reshape, iqtree_repeat_dir_list_reshape):
        total += 1
        concatenation(repeat_dir_list_each_repeat, out_dir, total)
        alstra(iqtree_repeat_dir_list_each_repeat, out_dir, total)


def infer_all():
    result_dir_list = [\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400/',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80/',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400/',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80/',
    ]
    iqtree_dir_list = [\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400/',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80/',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400/',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80/',
    ]
    out_dir_list = [\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400/',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80/',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400/',\
    'species_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80/',\
    ]
    repeat_num = 20
    for result_dir, iqtree_dir, out_dir in zip(result_dir_list, iqtree_dir_list, out_dir_list):
        infer_a_scenario(result_dir, iqtree_dir, out_dir, repeat_num)


if __name__ == '__main__':
    infer_all()


