import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import numpy as np
import random
from matplotlib.backends.backend_pdf import PdfPages


def plot_a_result(result_dir, repeat_size, repeat_num, ax):
    total_topology_list = []
    for each_result_base in os.listdir(result_dir):
        each_result_dir = os.path.join(result_dir, each_result_base)
        print(each_result_base)
        for repeat_base in os.listdir(each_result_dir):
            repeat_file = os.path.join(each_result_dir, repeat_base)
            with open(repeat_file) as repeat_handle:
                for line in repeat_handle:
                    topology_label = line.strip().split()[1].strip('n')
                    total_topology_list.append(topology_label)
    random.shuffle(total_topology_list)
    print(len(total_topology_list))
    topology_freq_dict_reshape = defaultdict(list)
    for i in range(repeat_num):
        each_repeat = []
        for j in range(repeat_size):
            t = total_topology_list.pop()
            each_repeat.append(t)
        counter = Counter(each_repeat)
        for topology, count in counter.items():
            topology_freq_dict_reshape[topology].append(count/repeat_size)



    # df = pd.DataFrame(topology_freq_dict_reshape)
    # print(df)
    # # df = df.reindex(columns=['0','1','2'])
    # df = df.reindex(columns=['1','2','3'])
    # print(df)
    
    # df.plot(kind='box', showmeans=True)
    # df.plot(kind='bar')
    
    # plt.figure(figsize=[5,8])
    x_list = []
    y_list = []
    yerr_list = []
    color_list = ["#75b84f", "cornflowerblue", "#ffad0f",]
    for topology, freq_list in sorted(topology_freq_dict_reshape.items(), key=lambda i:int(i[0])):
        x_list.append(topology)
        y_list.append(np.mean(freq_list))
        yerr_list.append(np.std(freq_list, ddof=1))
    print(x_list)
    print(y_list)
    print(yerr_list)
    print(result_dir)
    for i, x in enumerate(x_list):
        ax.bar(x, y_list[i], yerr=yerr_list[i], width=0.8, color=color_list[i], edgecolor="", ecolor="#202020", capsize=6, error_kw={"elinewidth":1.5})
    ax.set_ylim([0.28, 0.41])
    f_concordant = 1 - 2*np.exp(-0.01)/3
    f_disconcordant = np.exp(-0.01)/3
    x_down, x_up = ax.get_xlim()
    ax.plot([x_down, x_up], [f_concordant]*2, color="#75b84f", linestyle='--', linewidth=1.5)
    ax.plot([x_down, x_up], [f_disconcordant]*2, color="dimgray", linestyle='--', linewidth=1.5)
    ax.set_xticklabels(['(1)','(2)','(3)'])
    ax.set_xlim([x_down, x_up])


def plot_all_results():
    fig, axes = plt.subplots(1, 4, sharey=True)
    para_list = [
        ['sim_results_three_species_1bp_biallelic/s=0.0_200_200_200', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=0.0_1000_1000_40', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=-0.0075_200_200_200', 6000, 20],
        ['sim_results_three_species_1bp_biallelic/s=-0.0075_1000_1000_40', 6000, 20],
    ]

    biallelic_fig_file = 'biallelic.pdf'

    with PdfPages(biallelic_fig_file) as biallelic_fig_handle:   
        for para, ax in zip(para_list, axes):
            print(para)
            plot_a_result(*(para+[ax]))
        biallelic_fig_handle.savefig()



if __name__ == '__main__':
    plot_all_results()