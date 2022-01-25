#README
#This file sends a task to each of the 64 GPUS. The task is 'parameter_search_input_df.sh' which contains instructions to run
##a python code 'parameter_search_input_df.py' which executes a turing analysis of a dataframe.
# circuit_n = $1

cd /end/home/moliver/Documents/modelling/6eq/parameter_space_search/code
module load python
for ID in {0..63}
 do
   qsub -q long -e output/error.${ID} -o output/output.${ID} -t 1%64 parameter_search_input_df.sh -F "$1 $ID"
 done
