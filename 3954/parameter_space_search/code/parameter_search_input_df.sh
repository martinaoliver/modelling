#README
#This bash file contains instructions to run a python code 'parameter_search_input_df.py'
##which executes a turing analysis of a dataframe.cd /end/home/moliver/Documents/modelling/6eq/parameter_space_search/code
#$1 circuit_n, $2 ID

cd /end/home/moliver/Documents/modelling/6eq/parameter_space_search/code
module load python
python parameter_search_input_df.py $1 $2
