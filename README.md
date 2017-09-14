# Naqvi17-code
generate_tscan_kmers.py takes miRNA family seed sequences as input, and shuffles them, providing a specified number of shuffled sequences per miRNA family. In the manuscript, 6 shuffled 7mers per true seed were used

count_targetscan70_cons_output.py takes the output of targetscan_70.pl (available from the TargetScan website) as input, and counts gene-miRNA interactions with target sites conserved between human and a specified species

plot_transitions.py takes a basewise 3'UTR branch length file (provided as bls_big_outfile.txt) and a desired transcript, segments the transcript, and plots the output. Segmented UTRs for all 3'UTRs in bls_big_outfile.txt are provided in big_transitions.txt
