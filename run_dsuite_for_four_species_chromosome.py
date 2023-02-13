import os, glob, subprocess, multiprocessing, shutil
import random


def merge_for_a_repeat(repeat_dir, sample_name_dir):
    original_working_dir = os.getcwd()
    os.chdir(repeat_dir)
    reheadered_vcf_gz_base_list = []
    # print(repeat_dir)
    repeat_base = os.path.basename(repeat_dir)
    merged_vcf_gz = repeat_base+'_merged.vcf.gz'
    if not os.path.isfile(merged_vcf_gz):
        for vcf_base in glob.glob('*.vcf'):
            sample_name_file = os.path.join(sample_name_dir, vcf_base.replace('.vcf', '_sample.txt'))
            print(sample_name_file)
            reheadered_vcf_base = vcf_base.replace('.vcf', '_reheadered.vcf')
            reheadered_vcf_gz_base = vcf_base.replace('.vcf', '_reheadered.vcf.gz')
            p_reheader = subprocess.Popen(['bcftools', 'reheader', '-s', sample_name_file, '-o', reheadered_vcf_base, vcf_base], stdout=subprocess.PIPE)
            p_reheader.communicate()
            p_bgzip = subprocess.Popen(['bgzip', '-f', reheadered_vcf_base], stdout=subprocess.PIPE)
            p_bgzip.communicate()
            p_tabix = subprocess.Popen(['tabix', '-f', reheadered_vcf_gz_base], stdout=subprocess.PIPE)
            p_tabix.communicate()
            reheadered_vcf_gz_base_list.append(reheadered_vcf_gz_base)
        p_merge = subprocess.Popen(['bcftools', 'merge', '-Oz', '-0', '-o', merged_vcf_gz]+reheadered_vcf_gz_base_list, stdout=subprocess.PIPE)
        p_merge.communicate()
    # else:
    #     print(merged_vcf_gz)
    merged_vcf_gz_abspath = os.path.abspath(merged_vcf_gz)
    os.chdir(original_working_dir)
    return merged_vcf_gz_abspath

def run_dsuite(result_root, sample_name_dir, population_file, tree_file, out_root):
    group_num = 20
    result_root = os.path.abspath(result_root)
    sample_name_dir = os.path.abspath(sample_name_dir)
    working_dir = os.path.abspath('.')
    merged_vcf_gz_list = []
    for repeat_index, repeat_base in enumerate(sorted(os.listdir(result_root),key=lambda i:int(i.split('_')[-1]))):
        if int(repeat_base.split('_')[-1]) != repeat_index+1:
            print(repeat_base)
            raise Exception()
        repeat_dir = os.path.join(result_root, repeat_base)
        if any(os.listdir(repeat_dir)):
            merged_vcf_gz = merge_for_a_repeat(repeat_dir, sample_name_dir)
            merged_vcf_gz_list.append(merged_vcf_gz)
        else:
            print('Empty:', repeat_dir)
    random.shuffle(merged_vcf_gz_list)
    base_name = os.path.basename(result_root)
    group_size = int(len(merged_vcf_gz_list)/group_num)
    for i in range(group_num):
        name = base_name + '_' + str(i+1)
        out_dir = os.path.join(out_root, name)
        # print(out_dir)
        begin_data = i*group_size
        end_data = (i+1)*group_size
        data_in_a_group = merged_vcf_gz_list[begin_data:end_data]
        print(data_in_a_group)
        print(group_size, len(merged_vcf_gz_list))
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
            shutil.copy(population_file, out_dir)    
            used_population_file = os.path.join(out_dir, os.path.basename(population_file))
            p_dsuite = subprocess.Popen(['/home/hc/Dsuite/utils/DtriosParallel', '-n', name, '-t', tree_file, used_population_file]+data_in_a_group, stdout=subprocess.PIPE)
            p_dsuite.communicate()

def run_all():
    result_root = 'sim_results_four_species_chromosome_u=2.4e-6/neutral_400_400_400'
    sample_name_dir = 'sample_name_files'
    population_file = 'sim_four_species_chromosome.population'
    tree_file = 'sim_four_species_chromosome_species_tree.tre'
    out_root = 'Dsuite_results_four_species_chromosome_u=2.4e-6/neutral_400_400_400'
    run_dsuite(result_root, sample_name_dir, population_file, tree_file, out_root)

    result_root = 'sim_results_four_species_chromosome_u=2.4e-6/neutral_2000_2000_80'
    sample_name_dir = 'sample_name_files'
    population_file = 'sim_four_species_chromosome.population'
    tree_file = 'sim_four_species_chromosome_species_tree.tre'
    out_root = 'Dsuite_results_four_species_chromosome_u=2.4e-6/neutral_2000_2000_80'
    run_dsuite(result_root, sample_name_dir, population_file, tree_file, out_root)

    result_root = 'sim_results_four_species_chromosome_u=2.4e-6/deleterious_400_400_400'
    sample_name_dir = 'sample_name_files'
    population_file = 'sim_four_species_chromosome.population'
    tree_file = 'sim_four_species_chromosome_species_tree.tre'
    out_root = 'Dsuite_results_four_species_chromosome_u=2.4e-6/deleterious_400_400_400'
    run_dsuite(result_root, sample_name_dir, population_file, tree_file, out_root)

    result_root = 'sim_results_four_species_chromosome_u=2.4e-6/deleterious_2000_2000_80'
    sample_name_dir = 'sample_name_files'
    population_file = 'sim_four_species_chromosome.population'
    tree_file = 'sim_four_species_chromosome_species_tree.tre'
    out_root = 'Dsuite_results_four_species_chromosome_u=2.4e-6/deleterious_2000_2000_80'
    run_dsuite(result_root, sample_name_dir, population_file, tree_file, out_root)


if __name__ == '__main__':
    run_all()