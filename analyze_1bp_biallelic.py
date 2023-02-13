import os, time, random
import multiprocessing, subprocess
import tskit
import numpy as np


def sim_a_gene_tree(slim_script, trees_file, sel, pop_size_1, pop_size_2, pop_size_3):
    seed = (((int(time.time()) * os.getpid())**2) % 31324) * 1000 + int((time.time()*os.getpid() - int(time.time()*os.getpid()/1000)*1000))
    slim_seed = (((int(time.time()) * os.getpid())**2) % 31324) * 1731 + int((time.time()*os.getpid() - int(time.time()*os.getpid()/799)*799))
    print('Seed:', seed)
    print('SLiM seed:', slim_seed)
    random.seed(seed)
    p = subprocess.Popen(['slim', '-s', str(slim_seed), '-d', 'O=\'{}\''.format(trees_file), '-d', 'sel={}'.format(sel), '-d', 'pop_size_1={}'.format(pop_size_1), '-d', 'pop_size_2={}'.format(pop_size_2), '-d', 'pop_size_3={}'.format(pop_size_3), slim_script], stdout=subprocess.PIPE)
    p.communicate()
    ts = tskit.load(trees_file)
    sample_list = []
    for pop in ts.populations():
        a = ts.samples(pop.id)
        if len(a) > 0:
            sample = random.choice(a)
            sample_list.append(sample)
    # print(sample_list)
    sample_tree_record = ts.simplify(sample_list)
    sample_tree = sample_tree_record.first()
    outer_label = sample_tree.rank().label
    # outer_lineage = sample_list[outer_label]
    outer_population = str(outer_label+1)
    newick_str = str(sample_tree.as_newick(include_branch_lengths=False))
    try:
        var = next(sample_tree_record.variants())
        genotypes = var.genotypes
    except StopIteration:     # no variant
        genotypes = [0, 0, 0]   
    # print(genotypes)

    stem_mut_list = []
    leaf_mut_list = []
    first_coalescence_node = sample_tree_record.node(3)
    if first_coalescence_node.time > 3204:
        for mutation in sample_tree_record.mutations():
            mut_time = mutation.time
            # print(mutation.node, mutation.time)
            if mut_time > 3204:
                mut_node_id = mutation.node
                print('ancestral mutation:', mut_node_id, mutation.time)
                if mut_node_id < 3:
                    leaf_mut_list.append(mutation)
                elif mut_node_id == 3:
                    stem_mut_list.append(mutation)

    print('Out population:', outer_label + 1)
    out_mut_count = 0
    in_mut_count = 0
    print(leaf_mut_list)
    for mutation in leaf_mut_list:
        mut_node_id = mutation.node
        if mut_node_id == outer_label:
            out_mut_count += 1
        elif mut_node_id != outer_label and mut_node_id < 3:
            in_mut_count += 1
        else:
            raise Exception('impossible leaf_mut_node_id >= 3')

    return newick_str, genotypes, outer_population, out_mut_count, in_mut_count


def a_run(running, slim_script, trees_file, para_list):
    while True:
        try:
            newick_str, genotypes, outer_population, out_mut_count, in_mut_count = running(slim_script, trees_file, *para_list)
            if len(set(genotypes)) > 1:    # only when extant lineages are not all identical the gene tree is sampled
                print('Selected:', newick_str, genotypes, outer_population, out_mut_count, in_mut_count)
                break
        except ValueError:     # multiple roots
            pass
        except BrokenPipeError:
            pass
    return newick_str, outer_population, out_mut_count, in_mut_count


def a_scenario(running, slim_script, trees_file, para_list, repeat_num, tree_num_in_each_repeat, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    for i in range(repeat_num):
        out_file = os.path.join(out_dir, "repeat_"+str(i+1)+".txt")
        print('Analyzing:', out_file)
        if not os.path.isfile(out_file):
            with open(out_file, "w") as out_handle:
                for j in range(tree_num_in_each_repeat):
                    print('No.', j+1)
                    newick_str, outer_population, out_mut_count, in_mut_count = a_run(running, slim_script, trees_file, para_list)
                    out_handle.write('\t'.join([newick_str, 'n'+outer_population, str(out_mut_count), str(in_mut_count)]))
                    out_handle.write('\n')
                    os.sync()
            print('OK!', out_file)


def run_all():
    slim_script = 'sim_three_species_1bp_biallelic'
    result_root = 'sim_results_three_species_1bp_biallelic'
    if not os.path.isdir(result_root):
        os.makedirs(result_root)
    subset_num = 30
    repeat_num = 4
    tree_num_in_each_repeat = 1000
    sel_list = [-7.5e-3, 0.0]
    population_size_combinations = [
        [200, 200, 200], 
        [1000, 1000, 40], 
    ]
    
    for sel in sel_list:
        for population_sizes in population_size_combinations:
            pool = multiprocessing.Pool(processes=subset_num)
            sel_result_root = os.path.join(result_root, "_".join(["s={}".format(str(sel))]+[str(pop_size) for pop_size in population_sizes]))
            for i in range(subset_num):
                out_dir_base = "_".join([str(i+1)] + [str(pop_size) for pop_size in population_sizes])
                out_dir = os.path.join(sel_result_root, out_dir_base)
                trees_file = './tmp_output_1bp_biallelic_{}.trees'.format(i+1)
                para_list = [sel] + population_sizes
                time.sleep(1)
                # a_scenario(sim_a_gene_tree, slim_script, trees_file, para_list, repeat_num, tree_num_in_each_repeat, out_dir)
                pool.apply_async(a_scenario, (sim_a_gene_tree, slim_script, trees_file, para_list, repeat_num, tree_num_in_each_repeat, out_dir))
            pool.close()
            pool.join()


if __name__ == '__main__':
    run_all()