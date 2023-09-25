#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=96:mem=96gb

module load anaconda3/personal
source activate env1
cd $PBS_O_WORKDIR



threshold=100
n_samples=1000000
n_cpus=96
continue_search=true
seed=0
while $continue_search; do
    
    # Query the database for counts and return whether to continue search 
    continue_search=$(python database_scripts/system_class_counts_fromDb.py $threshold) #this function outputs wheather to continue simulation or not based on the counts of the system classes
    echo "Continue simulation: $continue_search"


    # Call the Python function to create input/parameter files
    python create_input/parameterfiles_creator_turinghill_variant8.py $n_samples $seed


    # Call the Python function lsa_df
    python analytical/run_lsa_multithread_seed_dependent.py $n_cpus 8 $seed $n_samples 


    #insert instabilities into df 
    python database_scripts/insert_analytical_output_instabilities.py 8 $seed $n_samples

    seed=$((seed + 1)) 
    # else
    #     echo "Counts are not below the threshold. No action needed."
    # fi

    # # Pause for a while before the next iteration
    # sleep 1  # Sleep for 30 seconds for db to recover when querying next
done

#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=96:mem=96gb