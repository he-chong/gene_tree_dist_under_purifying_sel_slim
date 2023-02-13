## gene_tree_dist_under_purifying_sel_slim

This folder contains Python scripts, SLiM scripts, and a Dsuite control file for the analyses of the paper 
'Phylogenetic reconstruction and introgression detection can be influenced by population sizes under 
purifying selection'

'analyze_1bp_biallelic.py' runs the SLiM simulation involving biallelic loci. 'sim_three_species_1bp_biallelic'
is the SLiM script used by 'analyze_1bp_biallelic.py'. The results generated by 'analyze_1bp_biallelic.py'
are the input of 'plot_biallelic_result.py' which generates Figure 2 of the paper. 

'analyze_1000bp_nucleotide.py' runs the SLiM simulation involving nucleotide sequences. 'sim_four_species_1000bp_nucleotide_neutral'
and 'sim_four_species_1000bp_nucleotide_deleterious' are the SLiM scripts used by 'analyze_1000bp_nucleotide.py'.
The nucleotide sequences generated by 'analyze_1000bp_nucleotide.py' were the input of 'run_iqtree.py'
which runs IQTREE to reconstruct gene trees. The true gene trees generated by 'analyze_1000bp_nucleotide.py'
and the inferred gene trees generated by 'run_iqtree.py' are the input of 'plot_iqtree_result.py' which
generates Figure 3 of the paper. 

'plot_iqtree_result.py' also generates the upper panels shown in Figure 5, which illustrate the influence
of population sizes on the value of delta under purifying selection.

'infer_species_tree.py' perform phylogenetic reconstruction using the nucleotide sequences generated by
'analyze_1000bp_nucleotide.py' and the inferred gene trees generated by 'run_iqtree.py'. The results of
'infer_species_tree.py' are the input of 'plot_species_tree.py', which generate Figure 4 of the paper.

'run_sim_four_species_chromosome.py' runs the SLiM simulation involving genetic variants. 'sim_four_species_chromosome_neutral'
and 'sim_four_species_chromosome_deleterious' are the SLiM scripts used by 'run_sim_four_species_chromosome.py'.
The VCF files generated by 'run_sim_four_species_chromosome.py' are the input of 'run_dsuite_for_four_species_chromosome.py',
which runs Dsuite to calculate the D-statistic. The results generated by 'run_dsuite_for_four_species_chromosome.py'
are the input of 'plot_dsuite_result.py', which generates the lower panels shown in Figure 5. Its results
illstrute the influence of population sizes on the value of D-statistic under purifying selection. Running Dsuite
needs information about the population/species where each sampled individual belongs to. This information is
written in 'sim_four_species_chromosome.population'

'plot_ancestral_mutation.py' generates the boxplots shown in Figure 7.

Using these files, one should be able to replicate the analyses described in the paper. Running these files
requres the installation of SLiM, tkits, pyslim, SciPy, NumPy, matplotlib