##########################
#########IMPORTS##########
##########################

import os
import os.path
import sys
sixeqpath ='/rds/general/user/mo2016/home/Documents/modelling/3954'
modulepath = sixeqpath + '/modules'
sys.path.append(modulepath)
while True:
    try:
        from linear_stability_analysis import *
        break
    except ImportError:
        sixeqpath ='/Volumes/mo2016/home/Documents/modelling/3954/'
        modulepath = sixeqpath + '/modules'
        sys.path.append(modulepath)
        from linear_stability_analysis import *
        # from randomfunctions import plot_all_dispersion

import sys
import time
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl

#######################
#########CODE##########
#######################


circuit_n = 6 #ID of circuit we want to analyse
n_species=2 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
batch_size = 10000
#obtain a dictionary with some parameters to use in our analysis

#turing1
turingI = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':5e-05, 'd_B':0.0025}

#noturing
nopattern = {'k1':100, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':5e-05, 'd_B':0.0025}
# df = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':1, 'd_B':0.0025}
# df = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':2.8,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':0.0025, 'd_B':0.0025}

# TuringI-hopf
turinghopf = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':0.3,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':5e-05, 'd_B':0.0025}

# hopf
hopf = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':0.3,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':0.1, 'd_B':0.0025}

set1 = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':0.5,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':0, 'd_B':0.0025}
set2 = {'k1':0.05, 'k2': 0.01, 'k3':2, 'n1':3, 'n2':3, 'n3':1, 'beta':0.01,'Va': 1, 'Vb':1, 'mua':1, 'mub':1, 'd_A':0.00005, 'd_B':0.0025}

parameterset_list = [turingI,nopattern,turinghopf,hopf,set1,set2]
parameterset_list_names = ['turingI','nopattern','turinghopf','hopf','set1','set2']

n=4
#Run analysis on 1M parameter sets
top_dispersion=6000
output_df = detailed_turing_analysis_dict(parameterset_list[n],circuit_n,n_species,top_dispersion=top_dispersion,calculate_unstable=True)
eigenvalues = output_df[4][0]
print(output_df)# print(np.shape(eigenvalues[0]))
print(np.shape(eigenvalues))
# mpl.rcParams['axes.spines.right']=False
# mpl.rcParams['axes.spines.top']=False
def plot_all_dispersion(eigenvalues, n_species=6, crop=5000, top=5000, L=100,name = 'general'):
    wvn_list = np.array(list(range(0, top + 1))) * np.pi / L
    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100
    color = ['limegreen','red']
    for n in range(n_species):

        plt.plot(wvn_list[:crop], eigenvalues.real[:crop,[n]],c=color[n], label=('Eigenvalue %s')%int(n+1))
        plt.plot(wvn_list[:crop], eigenvalues.imag[:crop,[n]], linestyle = '--',c=color[n])


    plt.xlabel('Wavenumber (k)', size=16)
    plt.ylabel('Eigenvalue', size=16)
    plt.ylim(-4,4)
    plt.xticks(size = '14')
    plt.yticks(size = '14')

    plt.axhline(y=0, color='green', linestyle='-', linewidth=0.1)
    plt.tight_layout()
    plt.savefig('dispersionrelation_%s.jpg'%name, dpi = 1500,bbox_inches='tight')
    plt.show()



plot_all_dispersion(eigenvalues,n_species=2,crop=top_dispersion,top=top_dispersion,name=parameterset_list_names[n])
