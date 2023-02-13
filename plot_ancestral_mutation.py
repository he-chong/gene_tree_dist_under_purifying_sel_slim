import os
import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib.backends.backend_pdf import PdfPages


def plot_a_result(result_dir, repeat_size, repeat_num, ax):
    color_list = ['royalblue', 'mediumseagreen']
    total_basal_count_list = []
    total_inner_count_list = []
    for each_result_base in os.listdir(result_dir):
        each_result_dir = os.path.join(result_dir, each_result_base)
        print(each_result_base)
        for repeat_base in os.listdir(each_result_dir):
            repeat_file = os.path.join(each_result_dir, repeat_base)
            with open(repeat_file) as repeat_handle:
                for line in repeat_handle:
                    basal_count = int(line.strip().split()[2])
                    inner_count = int(line.strip().split()[3])
                    total_basal_count_list.append(basal_count)
                    total_inner_count_list.append(inner_count)
    random.shuffle(total_basal_count_list)
    random.shuffle(total_inner_count_list)
    print(len(total_basal_count_list))
    print(len(total_inner_count_list))
    basal_freq_list = []
    inner_freq_list = []
    for i in range(repeat_num):
        basal_count_repeat = []
        inner_count_repeat = []
        for j in range(repeat_size):
            basal_count = total_basal_count_list.pop()
            inner_count = total_inner_count_list.pop()
            basal_count_repeat.append(basal_count)
            inner_count_repeat.append(inner_count)
        basal_sum = sum(basal_count_repeat)
        inner_sum = sum(inner_count_repeat)
        basal_freq = basal_sum / (basal_sum+inner_sum)
        inner_freq = inner_sum / (basal_sum+inner_sum)
        basal_freq_list.append(basal_freq)
        inner_freq_list.append(inner_freq)

    basal_inner_data = [basal_freq_list, inner_freq_list]
    boxplot = ax.boxplot(basal_inner_data, widths=[0.6]*len(basal_inner_data), patch_artist = True)
    for each, color in zip(boxplot['boxes'], color_list):
        each.set(color=color,  alpha=0.45)

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
    x_down, x_up = ax.get_xlim()
    ax.plot([x_down, x_up], [2/3]*2, color=color_list[0], linestyle='--', linewidth=1)
    ax.plot([x_down, x_up], [1/3]*2, color=color_list[1], linestyle='--', linewidth=1)
    ax.set_xticklabels(('A+B', 'C'))


def plot_ancestral_mutation():
    para_list = [
        ['sim_results_three_species_1bp_biallelic/s=0.0_200_200_200', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=0.0_1000_1000_40', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=-0.0075_200_200_200', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=-0.0075_1000_1000_40', 6000, 20],
    ]

    out_fig = 'ancestral_mutation.pdf'
    with PdfPages(out_fig) as out_handle:
        fig, axes = plt.subplots(1, len(para_list), figsize=(9, 3), sharey=True)
        for para, ax in zip(para_list, axes):
            plot_a_result(*(para+[ax]))
        out_handle.savefig()


if __name__ == '__main__':
    plot_ancestral_mutation()
