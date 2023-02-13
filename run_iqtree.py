import os, glob, multiprocessing, subprocess, time, shutil
from collections import defaultdict


def run_iqtree(sequence_fasta, out_dir):
    print(sequence_fasta)
    tree_file = sequence_fasta + '.treefile'
    if not os.path.exists(tree_file):
        p = subprocess.Popen(['iqtree', '-s', sequence_fasta, '-m', 'JC', '-o', 'n0', '-bb', '1000', '-nt', '2'], stdout=subprocess.PIPE)
        p.communicate()
        shutil.copy(tree_file, out_dir)


def run_a_scenario(result_root, out_root):
    for a_subset in sorted(os.listdir(result_root), key=lambda i: int(i.split('_')[0])):
        for repeat_base in os.listdir(os.path.join(result_root, a_subset, 'sequence_fasta')):
            repeat_dir = os.path.join(result_root, a_subset, 'sequence_fasta', repeat_base)
            out_dir = os.path.join(out_root, a_subset, repeat_base)
            sequence_fasta_list = sorted(glob.glob(os.path.join(repeat_dir, '*.fasta')), key=lambda i: int(i.split('_')[-1].split('.')[0]))
            if len(sequence_fasta_list) == 1000:
                if not os.path.isdir(out_dir):
                    os.makedirs(out_dir)
                pool = multiprocessing.Pool(processes=15)
                for sequence_fasta in sequence_fasta_list:
                    pool.apply_async(run_iqtree, (sequence_fasta, out_dir))
                pool.close()
                pool.join()


def run_all():
    result_root_list = [
        'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400',
        'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80',
        'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400',
        'sim_results_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80'
    ]
    out_root_list = [
        'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_400_400_400',
        'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/neutral_2000_2000_80',
        'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_400_400_400',
        'iqtree_trees_four_species_1000bp_nucleotide_u=2.4e-6/deleterious_2000_2000_80'
    ]
    for result_root, out_root in zip(result_root_list, out_root_list):
        run_a_scenario(result_root, out_root)


if __name__ == '__main__':
    run_all()