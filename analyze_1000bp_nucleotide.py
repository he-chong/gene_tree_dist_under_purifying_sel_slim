import os, time, subprocess, multiprocessing, random
import tskit, pyslim


def sim_a_gene_tree_nucleotide(slim_script, trees_file, pop_size_1, pop_size_2, pop_size_3):
    seed = (((int(time.time()) * os.getpid())**2) % 31324) * 1000 + int((time.time()*os.getpid() - int(time.time()*os.getpid()/1000)*1000))
    slim_seed = (((int(time.time()) * os.getpid())**2) % 31324) * 1731 + int((time.time()*os.getpid() - int(time.time()*os.getpid()/799)*799))
    print('Seed:', seed)
    print('SLiM seed:', slim_seed)
    random.seed(seed)
    p = subprocess.Popen(['slim', '-s', str(slim_seed), '-d', 'O=\'{}\''.format(trees_file), '-d', 'pop_size_1={}'.format(pop_size_1), '-d', 'pop_size_2={}'.format(pop_size_2), '-d', 'pop_size_3={}'.format(pop_size_3), slim_script], stdout=subprocess.PIPE)
    p.communicate()
    # time.sleep(1)
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
    outer_population = str(outer_label + 1)
    newick_str = str(sample_tree.as_newick(include_branch_lengths=False))
    try:
        var = next(sample_tree_record.variants())
        genotypes = var.genotypes
    except StopIteration:     # no variant
        genotypes = [0, 0, 0, 0]
    sample_tree_nucleotide = pyslim.convert_alleles(sample_tree_record)
    fasta_str = sample_tree_nucleotide.as_fasta()
    # print(fasta_str)

    stem_mut_list = []
    leaf_mut_list = []
    first_coalescence_node = sample_tree_record.node(4)
    if first_coalescence_node.time > 6416:
        for mutation in sample_tree_record.mutations():
            mut_time = mutation.time
            if mut_time > 6416:
                mut_node_id = mutation.node
                leaf = True if mut_node_id < 4 else False
                parent_of_leaves = True
                for child_id in sample_tree.children(mut_node_id):
                    if child_id >= 4:
                        parent_of_leaves = False
                        break
                if sample_tree.num_children(mut_node_id) == 0:
                    parent_of_leaves = False
                if leaf and parent_of_leaves:
                    print(sample_tree_record.mutation(mut_node_id))
                    print(sample_tree.children(mut_node_id))
                    raise Exception('Impossible situation')
                if leaf:
                    leaf_mut_list.append(mutation)
                elif parent_of_leaves:
                    stem_mut_list.append(mutation)

    print('Out population:', outer_label + 1)
    out_mut_count = 0
    in_mut_count = 0
    for mutation in leaf_mut_list:
        mut_node_id = mutation.node
        if mut_node_id > 0:
            if mut_node_id == outer_label + 1:
                out_mut_count += 1
            else:
                in_mut_count += 1

    return fasta_str, newick_str, genotypes, outer_population, out_mut_count, in_mut_count


def a_run(running, slim_script, trees_file, population_sizes):
    while True:
        try:
            fasta_str, newick_str, genotypes, outer_population, out_mut_count, in_mut_count = running(slim_script, trees_file, *population_sizes)
            print(newick_str, genotypes, outer_population, out_mut_count, in_mut_count)
            break
        except ValueError:     # multiple roots
            pass
        except BrokenPipeError:
            pass
    return fasta_str, newick_str, outer_population, out_mut_count, in_mut_count


def a_scenario(running, slim_script, trees_file, population_sizes, repeat_num, tree_num_in_each_repeat, out_root, filtering=False):
    true_topology_root = os.path.join(out_root, 'true_topology') 
    fasta_root = os.path.join(out_root, 'sequence_fasta')
    if not os.path.isdir(true_topology_root):
        os.makedirs(true_topology_root)

    for i in range(repeat_num):
        topology_file = os.path.join(true_topology_root, "repeat_"+str(i+1)+".txt")
        fasta_dir = os.path.join(fasta_root, "repeat_"+str(i+1))
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)
        if len(os.listdir(fasta_dir)) < tree_num_in_each_repeat:
            print('Analyzing:', fasta_dir)
            with open(topology_file, "w") as topology_handle:
                for j in range(tree_num_in_each_repeat):
                    print('No.', j+1)
                    fasta_str, newick_str, outer_population, out_mut_count, in_mut_count = a_run(running, slim_script, trees_file, population_sizes)
                    topology_handle.write('\t'.join([newick_str, 'n'+outer_population, str(out_mut_count), str(in_mut_count)]))
                    topology_handle.write('\n')
                    fasta_file = os.path.join(fasta_dir, 'sequences_{}.fasta'.format(j+1))
                    with open(fasta_file, 'w') as fasta_handle:
                        fasta_handle.write(fasta_str)
                print('OK!', topology_file)


def run_all():
    result_root = 'sim_results_four_species_1000bp_nucleotide_u=2.4e-6_count'
    if not os.path.isdir(result_root):
        os.makedirs(result_root)
    subset_num = 30
    repeat_num = 2
    tree_num_in_each_repeat = 1000
    sel_types = [
        'deleterious',
        'neutral',   
    ]
    population_size_combinations = [
        [400, 400, 400],
        [2000, 2000, 80]
    ]
    
    count_proc = 0
    for sel in sel_types:
        slim_script = 'sim_four_species_1000bp_nucleotide_' + sel
        for population_sizes in population_size_combinations:
            pool = multiprocessing.Pool(processes=subset_num)
            sel_result_root = os.path.join(result_root, '_'.join([sel] + [str(pop_size) for pop_size in population_sizes]))
            for i in range(subset_num):
                count_proc += 1
                out_dir_base = "_".join([str(i+1)] + [str(pop_size) for pop_size in population_sizes])
                out_dir = os.path.join(sel_result_root, out_dir_base)
                print(out_dir)
                trees_file = './tmp_output_nucleotide_{}.trees'.format(count_proc)
                time.sleep(1)
                pool.apply_async(a_scenario, (sim_a_gene_tree_nucleotide, slim_script, trees_file, population_sizes, repeat_num, tree_num_in_each_repeat, out_dir))
                # a_scenario(sim_a_gene_tree_nucleotide, slim_script, trees_file, population_sizes, repeat_num, tree_num_in_each_repeat, out_dir)
            pool.close()
            pool.join()


if __name__ == '__main__':
    run_all()