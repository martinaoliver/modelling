import pickle as pkl
import numpy as np

def load_dose_response(filename,HSL_transform =0.14*10**3):
    data = pkl.load(open('input/liquid_culture/curatedData/%s'%filename,'rb'))
    OC14_list= np.array(data['OC14'])*HSL_transform; gfpExp_list = list(data['mean_gfp']); rfpExp_list = list(data['mean_rfp'])
    semGreen = list(data['std_gfp']); semRed = list(data['std_rfp'])
    return OC14_list, gfpExp_list, rfpExp_list, semGreen, semRed

