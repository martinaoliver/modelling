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


circuit_n = 7 #ID of circuit we want to analyse
n_species=2 #number of molecular species in circuit_n (#Circuit2 has 6 molecular species)
#obtain a dictionary with some parameters to use in our analysis

#r3 and r2 dont seem to affect dispersion relation.
# open parameter dictionaries
i = {'name':'i','alpha':0.85, 'beta': -0.91, 'D':0.45, 'r3':1,  'r2':1, 'gamma':6}#turing I
ii = {'name':'ii','alpha':0.99, 'beta': -0.85, 'D':0.45, 'r3':1,  'r2':1, 'gamma':6}#THB II
iii = {'name':'iii','alpha':0.899, 'beta': -0.7, 'D':0.45, 'r3':1,  'r2':1, 'gamma':6}#hopf III
iv = {'name':'iv','alpha':0.55, 'beta': -0.7, 'D':0.45, 'r3':1,  'r2':1, 'gamma':6}#stable IV
v = {'name':'v','alpha':0.1, 'beta': -0.7, 'D':0.45, 'r3':1,  'r2':1, 'gamma':6}#stable IV
i_r0 = {'name':'i_r0','alpha':0.85, 'beta': -0.91, 'D':0.45, 'r3':0.001,  'r2':0.001, 'gamma':6}#turing I
fig3a = {'name':'fig3a','alpha':0.88, 'beta': -0.91, 'D':0.39, 'r3':3.05,  'r2':2, 'gamma':4}
fig3b = {'name':'fig3b','alpha':0.799, 'beta': -0.91, 'D':0.45, 'r3':0.1,  'r2':0.278, 'gamma':6}
fig3c = {'name':'fig3c','alpha':0.899, 'beta': -0.91, 'D':0.516, 'r3':3.5,  'r2':2, 'gamma':4}
fig3d = {'name':'fig3d','alpha':0.899, 'beta': -0.91, 'D':0.45, 'r3':0.1,  'r2':0.296, 'gamma':6}
thb1 = {'name':'thb1','alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':1.5,  'r2':1, 'gamma':6}
thb2 = {'name':'thb1','alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':2,  'r2':1.5, 'gamma':6}
thb3 = {'name':'thb1','alpha':0.95, 'beta': -0.91, 'D':0.45, 'r3':2.5,  'r2':1.75, 'gamma':6}
thbII = {'name':'test','alpha':0.95, 'beta': -0.91, 'D':0.0001, 'r3':2.5,  'r2':1.75, 'gamma':6} #thbII turingII hopf
thbII1 = {'name':'test','alpha':0.95, 'beta': -0.91, 'D':0.01, 'r3':2.5,  'r2':1.75, 'gamma':6} #thbII turingII hopf
par_dict_list = [i,ii,iii,iv,v,i_r0,fig3a,fig3b,fig3c,fig3d,thb1,thb2,thb3,thbII,thbII1]

#Run analysis on 1M parameter sets
top_dispersion=500
n=int(sys.argv[1])-1
print(par_dict_list[n])
output_df = detailed_turing_analysis_dict(par_dict_list[n],circuit_n,n_species,top_dispersion=top_dispersion,calculate_unstable=True)
eigenvalues = output_df[4][0]
print(output_df[3])
# print(output_df)# print(np.shape(eigenvalues[0]))
# print(np.shape(eigenvalues))
# print(np.shape(output_df))
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
    plt.title(name,size=8)
    plt.axhline(y=0, color='green', linestyle='-', linewidth=0.1)
    plt.tight_layout()
    plt.savefig('dispersionrelation_%s.jpg'%name['name'], dpi = 1500,bbox_inches='tight')
    plt.show()



plot_all_dispersion(eigenvalues,n_species=2,crop=top_dispersion,top=top_dispersion,name=par_dict_list[n])
