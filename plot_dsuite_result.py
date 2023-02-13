import os, glob
from collections import defaultdict
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_a_dsuite_result(result_dir, ax):
    p_value_list = []
    result_dict = defaultdict(list)
    plus_x_list = []
    plus_d_statistic_list = []
    minus_x_list = []    
    minus_d_statistic_list = []
    asterisk_x = []
    asterisk_y = []
    for x, result_base in enumerate(sorted(os.listdir(result_dir), key=lambda i: int(i.split('_')[-1]))):
        result_file = os.path.join(result_dir, result_base, 'DTparallel_sim_four_species_chromosome_{}_combined_tree.txt'.format(result_base))
        if os.path.exists(result_file):
            with open(result_file) as result_handle:
                result_line_split = result_handle.readlines()[1].strip().split()
                if result_line_split[0] =='s3':
                    abba, baba = [float(i) for i in result_line_split[-2:]]
                    d_statistic = float(result_line_split[3])
                    plus_x_list.append(x)
                    plus_d_statistic_list.append(d_statistic)
                else:
                    baba, abba = [float(i) for i in result_line_split[-2:]]
                    d_statistic = -float(result_line_split[3])
                    minus_x_list.append(x)
                    minus_d_statistic_list.append(d_statistic)
                result_dict['ABBA'].append(abba)
                result_dict['BABA'].append(baba)
                p_value = float(result_line_split[5])
                print(p_value)
                if p_value < 0.05:
                    print(result_base)
                    print(p_value)
                    print(d_statistic)
                    asterisk_x.append(x)
                    if d_statistic > 0:
                        asterisk_y.append(d_statistic+0.004)
                    else:
                        asterisk_y.append(d_statistic-0.004)
                p_value_list.append(p_value)
            
                
    # result_df = pd.DataFrame(result_dict)
    a = np.array(p_value_list)
    print(sum(np.where(a<0.05, 1, 0)))
    # d_statistic_series = pd.Series(d_statistic_list)
    print(len(plus_d_statistic_list)/(len(plus_d_statistic_list)+len(minus_d_statistic_list)))
    # print(result_df)
    # result_df.plot(kind='box')
    # d_statistic_series.plot(kind='bar', ax=ax)
    ax.bar(plus_x_list, plus_d_statistic_list, color='#ffad0f', width=0.6)
    ax.bar(minus_x_list, minus_d_statistic_list, color='cornflowerblue', width=0.6)
    y_down, y_up = ax.get_ylim()
    x_down, x_up = ax.get_xlim()
    ax.plot(asterisk_x, asterisk_y, linestyle='', marker='*', markersize=5, markeredgewidth=0.3, color='darkslategrey')
    ax.plot([x_down, x_up], [0, 0], linestyle='--')
    ax.set_ylim([-0.03, 0.03])


def plot_all_results():
    result_dir_list = [
    'Dsuite_results_four_species_chromosome_u=2.4e-6/neutral_400_400_400',
    'Dsuite_results_four_species_chromosome_u=2.4e-6/neutral_2000_2000_80',
    'Dsuite_results_four_species_chromosome_u=2.4e-6/deleterious_400_400_400',
    'Dsuite_results_four_species_chromosome_u=2.4e-6/deleterious_2000_2000_80',
    ]

    dsuite_fig_file = 'dsuite.pdf'

    with PdfPages(dsuite_fig_file) as dsuite_fig_handle:
        fig, axes = plt.subplots(1, len(result_dir_list), figsize=(12,4), sharey=True)
        for result_dir, ax in zip(result_dir_list, axes):
            plot_a_dsuite_result(result_dir, ax)
        dsuite_fig_handle.savefig()


if __name__ == '__main__':
    plot_all_results()
