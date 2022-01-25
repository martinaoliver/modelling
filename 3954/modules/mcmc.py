import pandas as pd
from randomfunctions import noise
from parametercombination_analysis import analysis_range, analysis_fromdict
import numpy as np

import pickle


def parameter_bounds(d):
    bounds = []
    Vm_range = (10,1000)
    b_range = (0.01,1)
    km_range = (0.1,250)
    mu_range = (0.001,50)
    for parameters,values in d.items():
        if str(parameters)[:1] == 'V':
            bounds.append(Vm_range)
        elif str(parameters)[:1] == 'b':
            bounds.append(b_range)
        elif str(parameters)[:1] == 'k':
            bounds.append(km_range)
        elif str(parameters)[:1] == 'm':
            bounds.append(mu_range)
        else:
            bounds.append([1,1])
    return bounds

def update_guess(xold,T,lower,upper):
    if lower == upper:
        return xold
    else:
        xnew=noise(xold,T)[0]
        return xnew

def propose_dict(d, T):
    dnew = {}
    for i, (parameters,xold) in enumerate(d.items()):
        bounds = parameter_bounds(d)
        lower = bounds[i][0]
        upper = bounds[i][1]
        xnew = update_guess(xold,T,lower,upper)
        dnew[parameters] = xnew
    return dnew

def metropolis(maxeig_list_new,maxeig_list_old):
    maxeig_new = max(maxeig_list_new,  default=0)
    maxeig_old = max(maxeig_list_old)

    if maxeig_new > maxeig_old:
        print ('Improved maxeig')
        return 1 #Accept update
    else:
        maxeig_ratio = maxeig_new/maxeig_old
        random_uniform_value = np.random.normal(1,0.01)
        print ('Random Value:' + str(random_uniform_value) + '   Maxeig Ratio:' +  str(maxeig_ratio))
        if maxeig_ratio > random_uniform_value:
            return 1 #Accept update
            print(1)
        else:
            return 0 #Reject update
            print(0)

def update_decision(maxeig_list_new, maxeig_list_old ,pattern_class_new, numbersteadystates_new):
    fails = ['unstablesteadystate', 'noinstability', 'turingII', 'travelling']
    fails_multiple_steadystates = ['noinstability', 'turingII', 'travelling']
    if numbersteadystates_new == 0:
        update = 0
        print('no steady states')
    elif numbersteadystates_new == 1:
        if fails[0] in pattern_class_new or fails[1] in pattern_class_new or fails[2] in pattern_class_new or fails[3] in pattern_class_new:
            update = 0
            print('No valid pattern type')
        else:
            update = metropolis(maxeig_list_new, maxeig_list_old)
            if update == 0:
                print ('move rejected')
            if update == 1:
                print ('move accepted')
            else:
                update=0
    elif numbersteadystates_new >1:
        if fails_multiple_steadystates[0] in pattern_class_new or fails_multiple_steadystates[1] in pattern_class_new or fails_multiple_steadystates[2] in pattern_class_new:
            update = 0
            print('No valid pattern type - Multiple steady states')
        else:
            update = metropolis(maxeig_list_new, maxeig_list_old)
            if update == 0:
                print ('move rejected')
            if update == 1:
                print ('move accepted')
            else:
                update=0
    return update

def dispersion_optimization(par_ID, T, breakafter,maxiter, path,  topic = ''):
    par_dict_first = pickle.load( open(path +  '/simulated_annealing/results/parameterfiles/%s.pkl'%par_ID, "rb" ) )
    steadystatelist_first, numbersteadystates_first ,maxeig_list_first, pattern_class_first, eigenvalues_first = analysis_fromdict(par_dict_first,path)

    if numbersteadystates_first == 1:
        optimization_path_df = pd.DataFrame()
        maxeig_list_iterations = []
        optimization_path_df = optimization_path_df.append(par_dict_first, ignore_index = True)
        maxeig_list_iterations.append(np.real(max(maxeig_list_first)))
        par_dict_old = par_dict_first
        maxeig_list_old = maxeig_list_first
        print ('Initial dispersion:' )
        print (maxeig_list_first)
        print ('--------')
        print ('')
        rejected_count = 0
        for n in range(maxiter):
            if rejected_count == breakafter:
                break

            else:
                print('Iteration: ' + str(n+1))
                par_dict_new = propose_dict(par_dict_old,T)
                steadystatelist_new, numbersteadystates_new ,maxeig_list_new, pattern_class_new , eigenvalues_new = analysis_fromdict(par_dict_new,path)
                update = update_decision(maxeig_list_new, maxeig_list_old, pattern_class_new, numbersteadystates_new)
                if update == 1:
                    par_dict_old = par_dict_new
                    maxeig_list_old = maxeig_list_new
                    optimization_path_df = optimization_path_df.append(par_dict_old, ignore_index = True)
                    print (np.real(maxeig_list_old))
                    rejected_count = 0
                else:
                    rejected_count +=1
                    print ('no update')
                maxeig_list_iterations.append(np.real(max(maxeig_list_old)))
                print ('--------')
                print ('')
    elif numbersteadystates_first > 1:
        print('Initially multiple steady states')

    # pickle.dump( optimization_path_df, open(path +  '/simulated_annealing/results/parameterfiles/optimization_path_df_%s_%s_T%r.pkl'%(par_ID, 'constant_temp' , T), "wb" ) )
    # pickle.dump( maxeig_list_iterations, open(path +  '/simulated_annealing/results/dispersion/maxeig_list_%s_%s_T%r.pkl'%(par_ID, 'constant_temp', T), "wb" ) )
    # pickle.dump( par_dict_old, open(path +  '/simulated_annealing/results/parameterfiles/%s_%s_T%r.pkl'%(par_ID, 'constant_temp', T), "wb" ) )


    pickle.dump( optimization_path_df, open('results/optimization_path_df_%s_%s_T%r.pkl'%(par_ID, 'constant_temp' , T), "wb" ) )
    pickle.dump( maxeig_list_iterations, open('results/maxeig_list_%s_%s_T%r.pkl'%(par_ID, 'constant_temp', T), "wb" ) )
    pickle.dump( par_dict_old, open('results/%s_%s_T%r.pkl'%(par_ID, 'constant_temp', T), "wb" ) )
