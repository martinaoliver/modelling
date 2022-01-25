#############################
#########README##############
#############################

This folder contains code to:
(1) Generate parameter sets using latin hypercube sampling (here parameter sets work with circuit 3954) in a loguniform distribution.
(2) Analyse the system with those parameter sets. This analysis consists of (2.1) finding steady states of the system
    and (2.2) analysing the stability of steady states using linear stability analysis.
    Recurring functions are stored a modules folder two directories above (6eq/modules)

To run the search in parallel (64 dataframes simultaneously), run in the linux command line the
file 'parallel_parameter_search_input_df.sh'. This will send a turing analysis of a whole dataframe to each of the 64 cores.
