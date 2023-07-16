#%%
import sys
import os
pwd = os.getcwd()
modellingpath = pwd.rpartition("modelling")[0] + pwd.rpartition("modelling")[1] 
sys.path.append(modellingpath + '/lib')
from database.databaseFunctions import *
import pickle
#############
#############
###paths#####
#############

import sys
import os
import pickle
import psycopg2
import matplotlib.pyplot as plt

from numerical.cn_plot import plot1D, surfpattern


#%%




L=50; dx =0.1; J = int(L/dx)
T =2000; dt = 0.02; N = int(T/dt)
boundaryCoeff=1;rate=L/T
suggesteddt = float(dx*dx*2)
growth = 'edgegrowth2'
simulation_param_dict = {'L':L, 'dx':dx, 'J':J, 'T':T, 'dt':dt, 'N':N, 
            'boundaryCoeff':boundaryCoeff, 
            'growth':growth, 'growth rate': rate}

parID = 1544038
circuit_n='turinghill'
variant= '9'
n_samples=2000000
ssID = 0
folder = f'{circuit_n}_variant{variant}'
model_param_dict = {'parID':parID, 'circuit_n':circuit_n,'variant':variant, 'n_samples':n_samples}
#%%
credentials=f"postgresql://moliver:moliver@ld-rendres07.bc.ic.ac.uk/moliver"

def query_simulationOutput_from_sql(sim_param_dict,model_param_dict,query_column, ssID=0):
    with psycopg2.connect(credentials) as conn:
        with conn.cursor() as cursor:
        # Build the SQL query dynamically based on the provided dictionaries
            query = """
            SELECT *
            FROM simulation_output o
            JOIN model_param m on o.model_param_id = m.model_param_id
            JOIN simulation_param s ON o.simulation_param_id = s.simulation_param_id
            WHERE 1=1
            """
            
            # Add filters for model_params
            for key, value in model_param_dict.items():
                query += f"AND m.{key} = '{value}'\n"
            
            # Add filters for simulation_params
            for key, value in simulation_param_dict.items():
                query += f"AND s.{key} = '{value}'\n"
            
            try:
                # Execute the query
                cursor.execute(query)
                
                # Fetch all the rows
                rows = cursor.fetchall()
                
                # Print or process the rows as per your requirement
                for row in rows:
                    print(row)
                    
            except Exception as e:
                print("Error executing query:", e)
            
            # Close the cursor and the connection
            cursor.close()
            conn.close()



            return simulationOutput



import psycopg2  # Assuming you are using PostgreSQL

def query_simulation_output(model_param_dict, simulation_param_dict):
    # Connect to your database
    conn = psycopg2.connect(
        host="your_host",
        database="your_database",
        user="your_username",
        password="your_password"
    )
    
    # Create a cursor object
    cursor = conn.cursor()
    
    # Build the SQL query dynamically based on the provided dictionaries
    query = """
    SELECT *
    FROM simulation_output o
    JOIN model_param m on o.model_param_id = m.model_param_id
    JOIN simulation_param s ON o.simulation_param_id = s.simulation_param_id
    WHERE 1=1
    """
    
    # Add filters for model_params
    for key, value in model_param_dict.items():
        query += f"AND m.{key} = '{value}'\n"
    
    # Add filters for simulation_params
    for key, value in simulation_param_dict.items():
        query += f"AND s.{key} = '{value}'\n"
    
    try:
        # Execute the query
        cursor.execute(query)
        
        # Fetch all the rows
        rows = cursor.fetchall()
        
        # Print or process the rows as per your requirement
        for row in rows:
            print(row)
            
    except Exception as e:
        print("Error executing query:", e)
    
    # Close the cursor and the connection
    cursor.close()
    conn.close()



#%
simulationOutput = query_simulationOutput_from_sql(simulation_param_dict, model_param_dict,query_column = 'U_final_1D', ssID=ssID)
plot1D(simulationOutput, savefig=False,filename='')

# %%
