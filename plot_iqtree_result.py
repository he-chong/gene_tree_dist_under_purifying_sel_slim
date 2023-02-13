import os, glob, re
from collections import Counter, defaultdict
# from scipy.stats import binom_test
from scipy.stats import chisquare
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def find_outer(tree_file):
    with open(tree_file) as tree_handle:
        tree_str_line = tree_handle.readline()
        stripped_tree_str = tree_str_line.strip().strip(';').strip('(').strip(')')
        in_parentheses_str = re.search(r'\(.+\)\d+', stripped_tree_str).group()
        in_parentheses_list = re.findall(r'n\d+', in_parentheses_str)
        bs = int(in_parentheses_str.split(')')[1])
        if 'n0' in in_parentheses_list:
            in_parentheses_list.remove('n0')
            outer = in_parentheses_list[0]
        else:
            outer = set(['n1', 'n2', 'n3']).difference(set(in_parentheses_list)).pop()
        # print(tree_str_line, outer)
    return outer, bs


def find_outer_no_bs(tree_file):
    with open(tree_file) as tree_handle:
        tree_str_line = tree_handle.readline()
        stripped_tree_str = tree_str_line.strip().strip(';').strip('(').strip(')')
        in_parentheses_str = re.search(r'\(.+\)', stripped_tree_str).group()
        in_parentheses_list = re.findall(r'n\d+', in_parentheses_str)
        if 'n0' in in_parentheses_list:
            in_parentheses_list.remove('n0')
            outer = in_parentheses_list[0]
        else:
            outer = set(['n1', 'n2', 'n3']).difference(set(in_parentheses_list)).pop()
    return outer


def analyze_a_repeat(repeat_dir_list, true_topology_file_list, bs_cutoff=50):
    # print(repeat_dir)
    # print(true_topology_file)
    true_outer_list = []
    for true_topology_file in true_topology_file_list:
        with open(true_topology_file) as true_topology_handle:
            for line in true_topology_handle:
                if line.strip():
                    true_outer = line.strip().split()[1].replace('#', 'n')
                    true_outer_list.append(true_outer)

    inferred_outer_list = []
    selected_list = []
    i = 0
    for repeat_dir in repeat_dir_list:
        for tree_file in sorted(glob.glob(os.path.join(repeat_dir, '*.treefile')), key=lambda i:int(os.path.basename(i).split('.')[0].split('_')[-1])):
            outer, bs = find_outer(tree_file)
            if bs >= bs_cutoff:
                inferred_outer_list.append(outer)
                selected_list.append(i)
            i += 1

    inferred_freq_dict = defaultdict(float)
    inferred_outer_counter = Counter(inferred_outer_list)
    for inferred_outer, count in inferred_outer_counter.items():
        inferred_freq_dict[inferred_outer] = count/len(inferred_outer_list)

    selected_true_outer_list = []
    for selected in selected_list:
        selected_true_outer_list.append(true_outer_list[selected])

    true_freq_dict = defaultdict(float)
    true_outer_counter = Counter(selected_true_outer_list)
    for true_outer, count in true_outer_counter.items():
        true_freq_dict[true_outer] = count/len(selected_true_outer_list)

    correct = 0
    incorrect = 0
    nonconcordant_2 = 0
    nonconcordant_3 = 0
    for inferred_outer, true_outer in zip(inferred_outer_list, selected_true_outer_list):
        if inferred_outer == 'n2':
            nonconcordant_2 += 1
        if inferred_outer == 'n3':
            nonconcordant_3 += 1
        if inferred_outer == true_outer:
            correct += 1
        else:
            incorrect += 1
    delta = (nonconcordant_3-nonconcordant_2) / (nonconcordant_3+nonconcordant_2)
    # p_value = binom_test(x=nonconcordant_3, n=nonconcordant_3+nonconcordant_2, p=0.5, alternative='two-sided')
    chi, p_value = chisquare([nonconcordant_2, nonconcordant_3])
    correct_rate = correct/(correct+incorrect)
    # print(len(selected_true_outer_list)/len(true_outer_list))
    return inferred_freq_dict, true_freq_dict, correct_rate, delta, p_value


def plot_results(inferred_top_freq_dict, true_top_freq_dict, correct_rate_list, true_ax, inferred_ax):
    color_list = ["#75b84f", "cornflowerblue", "#ffad0f"]
    # inferred_ax.set_title('Inferred topology frequencies')
    # true_ax.set_title('True topology frequencies')

    inferred_top_data = []
    true_top_data = []
    for top in ['n1','n2','n3']:
        inferred_top_data.append(inferred_top_freq_dict[top])
        true_top_data.append(true_top_freq_dict[top])

    inferred_boxplot = inferred_ax.boxplot(inferred_top_data, patch_artist = True)
    true_boxplot = true_ax.boxplot(true_top_data, patch_artist = True)
    for boxplot in (inferred_boxplot, true_boxplot):
        for each, color in zip(boxplot['boxes'], color_list):
            each.set(color=color, alpha=0.45)

        for each, color in zip(boxplot['medians'], color_list):
            each.set(color=color)

        for each, color in zip(boxplot['fliers'], color_list):
            each.set(marker='x', markeredgecolor=color, markersize=4)

        for key in ('whiskers', 'caps'):
            for i, each in enumerate(boxplot[key]):
                if i % 2:
                    each.set(color=color_list[int((i-1)/2)])
                else:
                    each.set(color=color_list[int((i)/2)])
    mean_correct_rate = round(sum(correct_rate_list)/len(correct_rate_list), 2)
    x = inferred_ax.get_xlim()[0] + (inferred_ax.get_xlim()[1] - inferred_ax.get_xlim()[0])*0.1
    inferred_ax.text(x, 0.42, 'Mean correct rate:\n{}'.format(mean_correct_rate), fontsize=8)

    for ax in [inferred_ax, true_ax]:
        ax.set_xticklabels(['(1)','(2)','(3)'])

    inferred_ax.set_ylim(0.24, 0.44)

    # plt.bar(range(len(delta_list)), delta_list)
    # plt.show()

def plot_delta(delta_list, p_value_list, ax):
    plus_x_list = []
    minus_x_list = []
    plus_delta_list = []
    minus_delta_list = []
    asterisk_x = []
    asterisk_y = []
    for i, delta in enumerate(delta_list):
        if delta >= 0:
            plus_x_list.append(i)
            plus_delta_list.append(delta)
            if p_value_list[i] < 0.05:
                asterisk_x.append(i)
                asterisk_y.append(delta+0.02)
        else:
            minus_x_list.append(i)
            minus_delta_list.append(delta)
            if p_value_list[i] < 0.05:
                asterisk_x.append(i)
                asterisk_y.append(delta-0.02)
    ax.bar(plus_x_list, plus_delta_list, color='#ffad0f', width=0.6)
    ax.bar(minus_x_list, minus_delta_list, color='cornflowerblue', width=0.6)
    y_down, y_up = ax.get_ylim()
    x_down, x_up = ax.get_xlim()
    ax.plot([x_down, x_up], [0, 0], linestyle='--')
    ax.plot(asterisk_x, asterisk_y, linestyle='', marker='*', markersize=5, markeredgewidth=0.3, color='darkslategrey')
    ax.set_ylim([-0.24, 0.24])


def analyze_a_scenario(iqtree_result_root, slim_result_root, true_ax, inferred_ax, repeat_num):
    inferred_top_stat_dict = defaultdict(list)
    true_top_stat_dict = defaultdict(list)
    incorrect_count_dict = defaultdict(list)
    correct_rate_list = []
    delta_list = []
    p_value_list = []

    total_repeat_dir_list = []
    total_true_topology_file_list = []
    for subset_base in os.listdir(iqtree_result_root):
        subset_fasta_dir = os.path.join(iqtree_result_root, subset_base)
        for repeat_base in os.listdir(subset_fasta_dir):
            repeat_dir = os.path.join(subset_fasta_dir, repeat_base)
            true_topology_file = os.path.join(slim_result_root, subset_base, 'true_topology', repeat_base+'.txt')
            total_repeat_dir_list.append(repeat_dir)
            total_true_topology_file_list.append(true_topology_file)

    repeat_size = int(len(total_repeat_dir_list)/repeat_num)
    repeat_dir_list_reshape = []
    true_topology_file_list_reshape = []
    for i in range(repeat_num):
        repeat_dir_list_each_repeat = []
        true_topology_file_list_each_repeat = []
        for j in range(repeat_size):
            repeat_dir = total_repeat_dir_list.pop()
            true_topology_file = total_true_topology_file_list.pop()
            repeat_dir_list_each_repeat.append(repeat_dir)
            true_topology_file_list_each_repeat.append(true_topology_file)
        repeat_dir_list_reshape.append(repeat_dir_list_each_repeat)
        true_topology_file_list_reshape.append(true_topology_file_list_each_repeat)

    for repeat_dir_list_each_repeat, true_topology_file_list_each_repeat in zip(repeat_dir_list_reshape, true_topology_file_list_reshape):
        inferred_freq_dict, true_freq_dict, correct_rate, delta, p_value = analyze_a_repeat(repeat_dir_list_each_repeat, true_topology_file_list_each_repeat)
        for inferred_outer, inferred_count in inferred_freq_dict.items():
            inferred_top_stat_dict[inferred_outer].append(inferred_count)

        for true_outer, true_count in true_freq_dict.items():
            true_top_stat_dict[true_outer].append(true_count)

        correct_rate_list.append(correct_rate)
        delta_list.append(delta)
        p_value_list.append(p_value)

    plot_results(inferred_top_stat_dict, true_top_stat_dict, correct_rate_list, true_ax, inferred_ax)
    return delta_list, p_value_list


def analyze_all():
    iqtree_result_root_list = [\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400',\
    'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80',\
    ]

    slim_root_list = [\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400',\
    'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80',\
    ]
    repeat_num = 20
    iqtree_fig_file = 'iqtree.pdf'
    delta_fig_file = 'delta.pdf'

    delta_list_list = []
    p_value_list_list = []
    with PdfPages(iqtree_fig_file) as iqtree_fig_handle:        
        fig, axes = plt.subplots(2,4, sharey=True, figsize=(9,8))
        for i, (iqtree_result_root, slim_result_root) in enumerate(zip(iqtree_result_root_list, slim_root_list)):
            true_ax, inferred_ax = axes[:,i]
            print(i)
            print(iqtree_result_root, slim_result_root)
            delta_list, p_value_list = analyze_a_scenario(iqtree_result_root, slim_result_root, true_ax, inferred_ax, repeat_num)
            delta_list_list.append(delta_list)
            p_value_list_list.append(p_value_list)
        iqtree_fig_handle.savefig()

    with PdfPages(delta_fig_file) as delta_fig_handle:    
        fig, axes = plt.subplots(1, 4, figsize=(12,4), sharey=True)
        for delta_list, p_value_list, ax in zip(delta_list_list, p_value_list_list, axes):
            plot_delta(delta_list, p_value_list, ax)
        delta_fig_handle.savefig()


if __name__ == '__main__':
    analyze_all()
