


# threshold=100
# n_samples=10
# n_cpus=1
# continue_search=true
# seed=0


threshold=100
n_samples=1000000
n_cpus=96
continue_search=true
seed=0


while $continue_search; do
    


    #insert instabilities into df 
    python database_scripts/insert_analytical_output_instabilities.py 9 $seed $n_samples

    seed=$((seed + 1)) 
    # else
    #     echo "Counts are not below the threshold. No action needed."
    # fi

    # # Pause for a while before the next iteration
    # sleep 1  # Sleep for 30 seconds for db to recover when querying next
done

