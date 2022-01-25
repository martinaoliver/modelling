from parametercombination_analysis import *
import pickle

def single_parameter_sensitivity(parameter, direction, par_dict, stepsize, par_ID,iter_ID,  workingpath, topic):
    print ('------%s------'%direction)
    par_dict_df= pickle.load(open(workingpath + '/sensitivity_analysis/parameterfiles/par%s_optimised_64.pkl'%(par_ID), 'rb'))
    par_dict_initial = par_dict_df.iloc[iter_ID].to_dict()

    par_dict_analysis = analysis_fromdict(par_dict_initial, par_ID, workingpath, topic)
    fails_multiple_steadystates = ['noinstability', 'turingII', 'travelling']

    if par_dict_analysis[1]>1:
        if fails_multiple_steadystates[0] in par_dict_analysis[3] or fails_multiple_steadystates[1] in par_dict_analysis[3] or fails_multiple_steadystates[2] in par_dict_analysis[3]:
            print('multi-steadystate no pattern')
        else:
            par_dict_analysis[3] = ['stableturing']
    while par_dict_analysis[3] == ['stableturing']:
        if direction == 'up':
            # par_dict[parameter]  =  par_dict[parameter] + par_dict_initial[parameter]*stepsize
            par_dict[parameter]  =  par_dict[parameter] + par_dict[parameter]*stepsize
            if par_dict[parameter] + par_dict_initial[parameter]*stepsize > 3000:
                break
        elif direction == 'down':
            par_dict[parameter]  =  par_dict[parameter] - par_dict_initial[parameter]*stepsize
            # par_dict[parameter]  =  par_dict[parameter] - par_dict[parameter]*stepsize
            if par_dict[parameter] + par_dict_initial[parameter]*stepsize < -3000:
                break
        print (par_dict[parameter] )
        par_dict_analysis = analysis_fromdict(par_dict, par_ID, workingpath, topic)
    return par_dict[parameter]

def allparameters_OFAT_sensitivity(stepsize,par_dict,par_ID,iter_ID,workingpath,topic):
    constants = ('d_A', 'd_B', 'n')
    # constants = ('Va', 'Vb', 'Vc', 'Vd', 'Ve', 'Vf', 'ba', 'bb', 'bc', 'bd', 'be', 'bf', 'd_A', 'd_B', 'kaa', 'kbd', 'kce', 'kda', 'keb', 'kee', 'kfe', 'mua', 'mub', 'mulva', 'n')
    limit_dict = {}
    for item in (i for i in par_dict.keys() if not i in constants):
        print('-----------------------')
        print(item)
        print('-----------------------')
        par_dict_df= pickle.load(open(workingpath + '/sensitivity_analysis/parameterfiles/par%s_optimised_64.pkl'%(par_ID), 'rb'))
        par_dict = par_dict_df.iloc[iter_ID].to_dict()
        limit_down = single_parameter_sensitivity(item, 'down', par_dict,stepsize,par_ID,iter_ID, workingpath,topic)

        par_dict_df= pickle.load(open(workingpath + '/sensitivity_analysis/parameterfiles/par%s_optimised_64.pkl'%(par_ID), 'rb'))
        par_dict = par_dict_df.iloc[iter_ID].to_dict()
        limit_up = single_parameter_sensitivity(item, 'up', par_dict,stepsize,par_ID, iter_ID, workingpath,topic)

        limit_dict[item] = (limit_down, limit_up)

    return limit_dict

def calculate_limits_range(iter_ID,par_ID,stepsize,workingpath,folder_ID='2020-06-15', iter_sa = 64):
    topic = 'sensitity_analysis_test'

    par_dict_df= pickle.load(open(workingpath + '/sensitivity_analysis/parameterfiles/%s_optimised_%s.pkl'%(par_ID,iter_sa), 'rb'))
    par_dict = par_dict_df.iloc[iter_ID].to_dict()

    limit_dict = pickle.load(open(workingpath + '/sensitivity_analysis/results/limits/limits_%s/limits_%s_iter%r_step%r.pkl'%(folder_ID,par_ID,iter_ID,stepsize),'rb'))
    constants = ('d_A', 'd_B', 'n')
    parameter_value = [value for key, value in par_dict.items() if key not in constants]
    limits_range = []
    for count,i in enumerate(limit_dict.values()):
        limits_range.append(i[1]-i[0])
    return  limits_range
#same as above but size of parameter is corrected for, by dividing the robustness range by the magnitude of the parameter
def calculate_limits_range_sizecorrection(iter_ID,par_ID,stepsize,workingpath,folder_ID, iter_sa = 64):
    topic = 'sensitity_analysis_test'

    par_dict_df= pickle.load(open(workingpath + '/sensitivity_analysis/parameterfiles/%s_optimised_%s.pkl'%(par_ID,iter_sa), 'rb'))
    par_dict = par_dict_df.iloc[iter_ID].to_dict()


    limit_dict = pickle.load(open(workingpath + '/sensitivity_analysis/results/limits/limits_%s/limits_%s_iter%r_step%r.pkl'%(folder_ID,par_ID,iter_ID,stepsize),'rb'))
    constants = ('d_A', 'd_B', 'n')
    parameter_value = [value for key, value in par_dict.items() if key not in constants]
    limits_range = []


    for count,i in enumerate(limit_dict.values()):
        limits_range.append((i[1]-i[0])/parameter_value[count])
#         limits_range.append(i[1]-i[0])


    return  limits_range
