# %config Completer.use_jedi = False
import sys
import os
pwd = os.getcwd()
root = pwd.split("home", 1)[0]
modelling_home = root + 'home/Documents/modelling'
modelling_ephemeral = root + 'ephemeral/Documents/modelling'
modulepath = modelling_home + '/3954/modules'
sys.path.append(modulepath)

from numerical_solvers_variableboundary import *
#execution parameters
circuit_n=7
# variant=0
# parametersets_n = 10000
save_figure = False
tqdm_disable = False #disable tqdm
dimension = str(sys.argv[1])
n_species=2
parID=int(sys.argv[2])-1

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
thbII1 = {'name':'test','alpha':0.95, 'beta': -0.91, 'D':0.0001, 'r3':2.5,  'r2':1.75, 'gamma':6} #thbII turingII hopf
par_dict_list = [i,ii,iii,iv,v,i_r0,fig3a,fig3b,fig3c,fig3d,thb1,thb2,thb3,thbII,thbII1]
# name_list = ['i','ii','iii','iv','v', 'i_r0',fig3a,fig3b,fig3c,fig3d]


print('parID = ' + str(parID))
mechanism = 'general'
boundary_coef = 0#1 is open boundary and 0 is closed boundary
shape = 'no_growth'
growth = False
variant=0

# boundary_coef = 1 #1 is open boundary and 0 is closed boundary
# shape = 'growing_colony'

x_gridpoints = int(sys.argv[3])
T =int(sys.argv[4])
# par_dict = general_df.loc[(parID,0)].to_dict()
par_dict =par_dict_list[parID]
print(par_dict)
L=int(sys.argv[5])
J = L *x_gridpoints  # number of equally spaced gridpoints in space domain (larger J means more spatial precision(tends towards continuum solution) )
t_gridpoints = t_gridpoints_stability(L, J, T)  # number of equally spaced gridpoints in domain (larger N means more temporal precision (tends towards continuum solution) )

N = T * t_gridpoints
initial_condition = [0.001]*n_species

filename = 'circuit%r_variant%r_boundary%r_%s_%sID%r_L%r_J%r_T%r_N%r'%(circuit_n,variant,boundary_coef, shape,mechanism,parID,L,J,T,N)

if dimension == '1D':
    # try:
    records, final_concentration, grids = crank_nicolson(par_dict, initial_condition, L, J, T, N, circuit_n, n_species=n_species,boundary_coef=boundary_coef,growth = growth,tqdm_disable=tqdm_disable)
    print(np.shape(records))
    plot_1D_final_concentration(final_concentration, grids,mechanism,shape,filename,n_species,parID,save_figure=save_figure,path='')
    surfpattern(records,grids,par_dict,morphogen = 0)
    if save_figure ==True:
        pickle.dump(final_concentration, open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/1Dfinal_%s.pkl'%filename, 'wb'))
        pickle.dump(records,open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/1Dtimeseries_%s.pkl'%filename, 'wb'))
        # timeseries_1D_rgb(records, L, J,T,filename,modelling_hpc,save_figure=save_figure)

    # except ValueError:
    #     print('!!!!!!!!!!!!!')
    #     print('ValueError --> unstable solution')
    #     print('!!!!!!!!!!!!!')
    #     print()

        # pass


if dimension == '2D':
 # Define 2D numerical parameters
    L_x = L
    L_y = L
    I = J

    # Run 2D simulation
    try:
        records,final_concentration,grids = adi_shape(par_dict,initial_condition,L_x,L_y,J,I,T,N, circuit_n,shape, boundary_coef=boundary_coef,tqdm_disable=tqdm_disable)#,tqdm_disable=tqdm_disable)
        plot_redgreen_contrast(final_concentration,L,mechanism,shape,filename,modelling_ephemeral,dimension=dimension,scale_factor=x_gridpoints,save_figure=save_figure)
        if save_figure ==True:
            pickle.dump(final_concentration, open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/2Dfinal_%s.pkl'%filename, 'wb'))
            pickle.dump(records,open(modelling_ephemeral + '/6eq/numerical_confocal/results/simulation/1M_growingcolony/2Dtimeseries_%s.pkl'%filename, 'wb'))

        # else:
        #     plt.show()

    except ValueError:
        print('!!!!!!!!!!!!!')
        print('ValueError --> unstable solution')
        print('!!!!!!!!!!!!!')
        print()

        pass
