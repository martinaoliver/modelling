#!/bin/bash

threshold=46
n_samples=10
n_cpus=1
continue_search=true
seed=0
while $continue_search; do
    
    # Query the database for counts and return whether to continue search 
    continue_search=$(python database_scripts/system_class_counts_fromDb.py $threshold) #this function outputs wheather to continue simulation or not based on the counts of the system classes
    echo "Continue simulation: $continue_search"


    # Call the Python function to create input/parameter files
    python create_input/parameterfiles_creator_turinghill_variant9.py $n_samples $seed


    # Call the Python function lsa_df
    python analytical/run_lsa_multithread_seed_dependent.py $n_cpus 9 $seed $n_samples 


    seed=$((seed + 1)) 
    # else
    #     echo "Counts are not below the threshold. No action needed."
    # fi

    # # Pause for a while before the next iteration
    # sleep 1  # Sleep for 30 seconds for db to recover when querying next
done
