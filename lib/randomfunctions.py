import pandas as pd
import scipy.io
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from sklearn.preprocessing import StandardScaler
import pickle
from numpy import trapz
from scipy import signal

#####manipulating dataframes

def mattodf():
    dictlist = []
    numberlist = []
    source_path = '/Users/mo2016/Documents/mres/6eqpipeline/lhs_initial_sampling'
    id = [f for f in os.listdir(source_path + '/results_2020-05-01/parameterfilesturing/parameterfilesturing_mat')]
    for file in id:

        number = file[8:file.find('.')]
        numberlist.append(number)
        parameters = scipy.io.loadmat(source_path + '/results_2020-05-01/parameterfilesturing/parameterfilesturing_mat/' + file, squeeze_me=True)
        parameters.pop('__header__');parameters.pop('__version__');parameters.pop('__globals__')
        dictlist.append(parameters)
    parameterdf = pd.DataFrame(dictlist)
    parameterdf = setindex(parameterdf,numberlist)
    parameterdf = parameterdf.drop(columns=['circuit_n'])


    return parameterdf, numberlist

def mat_to_parpickle(par_number,workingpath):
    Xturing, numberlist= mattodf()
    Xturing.to_pickle(workingpath + '/lhs_initial_sampling/parameterfiles/parameters-2020-05-01.pkl')


    par = Xturing.loc[str(par_number)].to_dict()
    pickle.dump(par, open(workingpath + '/lhs_initial_sampling/results_2020-05-01/parameterfilesturing/cir1_par%r.pkl'%par_number, "wb" ) )
    return par

def readdf_withoutconstant():
    source_path = '/Users/mo2016/Documents/mres/6eqpipeline/lhs_initial_sampling'
    df = pd.read_pickle(source_path + '/parameterfiles/parameters-2020-05-01.pkl')
    df=df.drop(['d_A','d_B','n'],axis=1)
    return df

def setindex(df, numberlist):
    df['index'] = numberlist
    df = df.set_index('index')
    return df
def addconstants(df,size):
    cooperativity = np.full((size, 1), 2)
    d_A = np.full((size, 1), 0.789)
    d_B = np.full((size, 1), 0.511)
    df['d_A'], df['d_B'], df['n'] = (d_A,d_B,cooperativity)
    return df

####Obtain ID for patterns of a certain type

def list_patterntype(typepattern):
    source_path = '/Users/mo2016/Documents/mres/6eqpipeline/lhs_initial_sampling'
    list_id = []
    for file in os.listdir(source_path + '/results_2020-05-01/lsa_results/'):
            if file.endswith(typepattern+'.mat'):
                number = file[19:]
                number = number[:number.find('_')]
                list_id.append(number)
    return list_id

def list_patterntype_perturbation(perturbation_value,typepattern,turing_value,date):
    list_id = []
    for file in os.listdir('/Users/mo2016/Documents/mres/6eqpipeline/lhs_initial_sampling/perturbations/perturbations_%r_%s/results/lsa_results'%(turing_value,date)):
        if file.endswith(perturbation_value +  '_100_' + typepattern + '.mat'):
            number = file[19:]
            number = number[:number.find('_')]
            number = (number + '_%s'%perturbation_value)
            list_id.append(number)
    return list_id

#######Creating noise to explore parameter space around point
def noise(value,deviation,size=1):
    fractional_deviation = 1/deviation
    noise= np.random.uniform(value+value/fractional_deviation,value-value/fractional_deviation, size)
    return noise

def point_exploration_df(point,deviation,size):
    noise_dict = {}
    df = readdf_withoutconstant()
    for parameter in list(df):
        noise_dict[parameter]=noise(df.loc[point,parameter],deviation,size)

    new_df = pd.DataFrame(noise_dict)
    new_df = setindex(new_df,np.arange(1,size+1))
    new_df = new_df.append(df.loc[point])
    new_df = addconstants(new_df,size+1)

    return new_df

def updatevar(par_n, path = 'lhs_initial_sampling/perturbations/perturbations_835723_2020-05-12', file='parameters_835723_0.01_100.pkl'):
    parameters_df = pd.read_pickle(path+'/results/parameterfiles/%s'%file)
    parameters_df_dict = parameters_df.to_dict('index')
    parameters_dict = parameters_df_dict[par_n]
    return parameters_dict

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

def wavelenght_from_dispersion(eigenvalues,n_species = 6, top = 5000, L=100):
    wvn_list = np.array(list(range(0,top+1)))*np.pi/L
    row_position = np.argmax(eigenvalues[:,n_species-1], axis=0)
    wavelenght = (2*np.pi)/wvn_list[row_position]
    return wavelenght

def wavelenght_from_numerical(records,grids,morphogen_number):
    peaks, _ = signal.find_peaks(records[morphogen_number][:,-10], height=0)

    x_grid = grids[0]
    i=0
    x_peaks = np.zeros(len(peaks))
    for n in peaks:
        x_peaks[i] = x_grid[n]
        i+=1


    distance_peaks  = np.diff(x_peaks)
    wavelenght = np.mean(distance_peaks)
    wavenumber = 2*np.pi/wavelenght
    return wavelenght, wavenumber
def average_wavelenght_from_numerical(records,grids):
    wavelenght_list = []
    for n in range(6):
        wavelenght,wavenumber = wavelenght_from_numerical(records,grids,n)
        wavelenght_list.append(wavelenght)
    average_wavelenght = np.mean(wavelenght_list)
    return average_wavelenght


def plot_highest_dispersion(eigenvalues,crop = 1000, top = 5000, L=100):
    wvn_list = np.array(list(range(0,top+1)))*np.pi/L
    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100

    plt.plot(wvn_list[:crop], eigenvalues.real[:crop,[-1]], label='Real highest eigenvalue', c='k')
    plt.plot(wvn_list[1:crop], eigenvalues.imag[1:crop,[-1]], linestyle = '--', label = 'Imaginary highest eigenvalue', c='k')

    plt.legend()
    plt.xlabel('Wavenumber')
    plt.ylabel('Eigenvalue')
    plt.axhline(y=0, color='k', linestyle='-', linewidth = 0.1)
    plt.grid()
    plt.tight_layout()


def plot_hopf_dispersion(eigenvalues,crop = 1000, top = 5000, L=100):
    wvn_list = np.array(list(range(0,top+1)))*np.pi/L
    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100

    plt.plot(wvn_list[:crop], eigenvalues.real[:crop,[-1]], label='Real highest eigenvalue', c='k')
    plt.plot(wvn_list[1:crop], eigenvalues.imag[1:crop,[-1]], linestyle = '--', label = 'Imaginary highest eigenvalue', c='k')
    plt.plot(wvn_list[1:crop], eigenvalues.imag[1:crop,[-2]], linestyle = '--', label = 'Imaginary highest eigenvalue', c='k')

    plt.legend()
    plt.xlabel('Wavenumber')
    plt.ylabel('Eigenvalue')
    plt.axhline(y=0, color='k', linestyle='-', linewidth = 0.1)
    plt.grid()
    plt.tight_layout()


def plot_all_dispersion(eigenvalues, n_species=6, crop=False, top=5000, L=100):
    wvn_list = np.array(list(range(0, top + 1))) * np.pi / L
    # wvn_list = np.array(list(range(0,5000+1)))*np.pi/100
    real_dominant_eig = eigenvalues.real[:,-1]
    indexZeros = np.where(np.diff(np.sign(real_dominant_eig)))[0]
    if crop==False:
        crop=indexZeros[-1]+5
    for n in range(n_species):
        plt.plot(wvn_list[:crop], eigenvalues.real[:crop,[n]])
        plt.plot(wvn_list[:crop], eigenvalues.imag[:crop,[n]], linestyle = '--',c='k')
        # plt.plot(wvn_list[:crop], np.real(eigenvalues[n][:crop]))
        # plt.plot(wvn_list[:crop], np.imag(eigenvalues[n][:crop]), linestyle = '--',c='k')


    plt.xlabel('Wavenumber')
    plt.ylabel('Eigenvalue')
    plt.axhline(y=0, color='green', linestyle='-', linewidth=0.1)
    plt.grid()
    plt.tight_layout()
    # plt.show()

def standardise_df(df):
    standardised_array = StandardScaler().fit_transform(df.values)
    standardised_df = pd.DataFrame(standardised_array, index=df.index, columns=df.columns)
    return standardised_df

def calculate_eucledian_distance(df, n_runs = 64):
    initial_par = df.iloc[0]
    eucledian_distance = []
    for i in range (n_runs):
        final_par = df.iloc[i+1]
        dist = np.linalg.norm(initial_par-final_par)
        eucledian_distance.append(dist)
    return eucledian_distance

def eucledian_dict(X,Y):
    eucledian = math.sqrt(sum((X.get(d,0) - Y.get(d,0))**2 for d in set(X) | set(Y)))
    return eucledian

def mat_to_dict(parameterfile,source_path):
    par_dict = scipy.io.loadmat(source_path + '/cir1_%s.mat'%parameterfile , squeeze_me=True)
    par_dict.pop('__header__');par_dict.pop('__version__');par_dict.pop('__globals__')
    par_dict.pop('circuit_n')
    return par_dict




def calculate_max_eigenvalue(iter_ID,par_dict_df,par_ID,workingpath, dx=0.03141592653589792):
        topic = ''
        par_dict = par_dict_df.iloc[iter_ID].to_dict()
        par_dict_analysis = analysis_fromdict(par_dict, par_ID, workingpath, topic)
        max_eigenvalue = np.real(np.nanmax(par_dict_analysis[2][:]))

        # dispersion = list(par_dict_analysis[4][0][:,5])
        # realdispersion = np.real(dispersion)
        # positivedispersion = list(i for i in realdispersion if i>0)
        # instability_area = trapz(positivedispersion, dx=dx)

        # positive_dispersion_list = [num for num in dispersion if np.real(num) > 0]
        # width_instability = len(positive_dispersion_list)
        return max_eigenvalue#, instability_area, width_instability

def constant_allparameters_perturbation(par_dict,perturbation_scale=0.01):
    constants = ('d_A', 'd_B', 'n')
    par_dict_constant_down = dict(par_dict)
    par_dict_constant_up = dict(par_dict)


    par_dict_constant_down.update({key:value -value*perturbation_scale for (key,value) in par_dict.items() if key not in constants})
    par_dict_constant_up.update({key:value + value*perturbation_scale for (key,value) in par_dict.items() if key not in constants})

    return par_dict_constant_down, par_dict_constant_up


def delta_peak(par_dict, perturbation_scale = 0.01, workingpath = ''):
    par_dict_constant_down, par_dict_constant_up = constant_allparameters_perturbation(par_dict, perturbation_scale)

    turing_analysis_par_dict = analysis_fromdict(par_dict, workingpath)
    print('down')
    turing_analysis_par_dict_constant_down = analysis_fromdict( par_dict_constant_down, workingpath)
    print('up')
    turing_analysis_par_dict_constant_up = analysis_fromdict( par_dict_constant_up, workingpath)

    max_eigenvalue_par_dict = np.real(np.nanmax(turing_analysis_par_dict[2]))
    max_eigenvalue_par_dict_constant_down = np.real(np.nanmax(turing_analysis_par_dict_constant_down[2]))
    max_eigenvalue_par_dict_constant_up = np.real(np.nanmax(turing_analysis_par_dict_constant_up[2]))

    relative_delta_peak_down = abs(max_eigenvalue_par_dict - max_eigenvalue_par_dict_constant_down)/max_eigenvalue_par_dict
    relative_delta_peak_up = abs(max_eigenvalue_par_dict - max_eigenvalue_par_dict_constant_up)/max_eigenvalue_par_dict
    average_delta_peak = (relative_delta_peak_down + relative_delta_peak_up)/2


    return average_delta_peak
