import os, multiprocessing, subprocess, time


def run_slim(slim_script, trees_file, out_dir, population_sizes):
    pop_size_1, pop_size_2, pop_size_3 = population_sizes
    slim_seed = (((int(time.time()) * os.getpid())**2) % 31324) * 1731 + int((time.time()*os.getpid() - int(time.time()*os.getpid()/799)*799))
    print('SLiM seed:', slim_seed)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    out_vcf_0 = os.path.join(out_dir, 'species_0.vcf')
    out_vcf_1 = os.path.join(out_dir, 'species_1.vcf')
    out_vcf_2 = os.path.join(out_dir, 'species_2.vcf')
    out_vcf_3 = os.path.join(out_dir, 'species_3.vcf')
    p = subprocess.Popen(['slim', '-s', str(slim_seed), '-d', 'O=\'{}\''.format(trees_file), \
        '-d', 'pop_size_1={}'.format(pop_size_1), '-d', 'pop_size_2={}'.format(pop_size_2), '-d', 'pop_size_3={}'.format(pop_size_3), \
        '-d', 'out_vcf_0=\'{}\''.format(out_vcf_0), '-d', 'out_vcf_1=\'{}\''.format(out_vcf_1), '-d', 'out_vcf_2=\'{}\''.format(out_vcf_2), '-d', 'out_vcf_3=\'{}\''.format(out_vcf_3), slim_script], stdout=subprocess.PIPE)
    p.communicate()


def run_all():
    result_root = 'sim_results_four_species_chromosome_u=2.4e-6'
    proc_num = 12
    repeat_num = 600
    sel_types = [
        'deleterious',
        'neutral', 
    ]
    population_size_combinations = [
        [400, 400, 400],
        [2000, 2000, 80]
    ]
    pool = multiprocessing.Pool(processes=proc_num)
    count_proc = 0
    for sel in sel_types:
        slim_script = 'sim_four_species_chromosome_' + sel
        for population_sizes in population_size_combinations:
            sel_result_root = os.path.join(result_root, '_'.join([sel] + [str(pop_size) for pop_size in population_sizes]))
            for i in range(repeat_num):
                out_dir = os.path.join(sel_result_root, 'repeat_{}'.format(i+1))
                trees_file = os.path.join(out_dir, 'sim_four_species_chromosome.trees')
                if not os.path.exists(trees_file):
                    pool.apply_async(run_slim, (slim_script, trees_file, out_dir, population_sizes))
    pool.close()
    pool.join()


if __name__ == '__main__':
    run_all()