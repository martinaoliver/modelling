#############################
#IMPORTS#
#############################
import sys
import os

from numpy import searchsorted

pwd = os.getcwd()
root = pwd.rpartition("mo2016")[0] + pwd.rpartition("mo2016")[1] #/Volumes/mo2016/ or '/Users/mo2016/' or '/rds/general/mo2016/'
if root == '/Users/mo2016':
    modelling_ephemeral = '/Volumes/mo2016/ephemeral/Documents/modelling'
    modelling_home = '/Volumes/mo2016/home/Documents/modelling'
    modelling_local = root + '/Documents/modelling'


if root == '/Volumes/mo2016' or root=='/rds/general/user/mo2016': #'/rds/general' or root=='/Volumes':
        modelling_ephemeral = root + '/ephemeral/Documents/modelling'
        modelling_home = root  + '/home/Documents/modelling'
        modelling_local = modelling_home

if root == '/Users/mo2016' or  root == '/Volumes/mo2016':
    import matplotlib as mpl
    mpl.use('tkagg')

modulepath = modelling_local + '/3954/modules/new_CN'

sys.path.append(modulepath)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


#image
def plot_2D_final_concentration(final_concentration,L,J,filename,path, n_species=6, save_figure=False ):

    # A,B,C,D,E,F = final_concentration
    dx = float(L)/float(J-1)
    grid = np.array([j*dx for j in range(J)])
    labels=['A','B','C','D','E','F']

    fig, axs = plt.subplots(2,int(n_species/2))#,figsize=(7.5,4))
    ax = axs.flatten()
    black_yellow= [(45/255,45/255,45/255),(255/255,255/255,0/255)]
    black_red = [(45/255,45/255,45/255),(255/255,0/255,0/255)]
    black_green = [(45/255,45/255,45/255),(50/255,205/255,50/255)]
    yellow_cmap = LinearSegmentedColormap.from_list('black_yellow',black_yellow)
    red_cmap = LinearSegmentedColormap.from_list('black_red',black_red)
    green_cmap = LinearSegmentedColormap.from_list('black_green',black_green)
    cmap_list = [yellow_cmap,yellow_cmap,yellow_cmap,green_cmap,red_cmap,green_cmap]
    ims = [0]*n_species
    for n in range(n_species):
        specie = final_concentration[n]
        ims[n] = ax[n].pcolormesh(grid, grid, specie, shading='auto', label = labels[n],cmap=cmap_list[n])
    for ax in axs.flat:
        ax.label_outer()
    #
    count1=0
    morphogens = ('A','B','C','D','E','F')
    for ax in axs.flat:
        ax.set(title=morphogens[count1])
        fig.colorbar(ims[count1], ax=ax)

        count1+=1

    fig.tight_layout()
    if save_figure == True:
        plt.savefig(path + '/%s_%s.jpeg' % ('2D', filename),dpi=2000)
        plt.close()
    else:
        plt.show()


def matrix_rgb_normalisation(matrix):
    row_n = 0
    NewMatrix = np.zeros(matrix.shape)

    OldMin = np.min(matrix[np.nonzero(matrix)])
    OldMax = np.amax(matrix[np.nonzero(matrix)])+0.0001 #Add 0.0001 so that patterns with no var dont give errors
    if OldMin < 0:
        print('WARNING: Negative numbers!!!!!')
    # NewMin = 1 #make newmin 1 instead of zero so 1 can represent cells
    NewMin = 0
    NewMax = 255
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)

    for row in matrix:
        column_n = 0
        for value in row:
            if value!=0:
                NewMatrix[column_n, row_n] = int((((value- OldMin) * NewRange) / OldRange) + NewMin)
            column_n += 1
        row_n += 1
    return NewMatrix, OldMin, OldMax

def plot_redgreen_contrast(final_concentration, mm, mechanism, shape, filename, path, parID=0, scale_factor=10, save_figure=False, dimension='2D'):
    green = final_concentration[-1]
    red = final_concentration[-2]
    x_grid = np.linspace(0, mm, len(green))
    normalised_red, redmin, redmax = matrix_rgb_normalisation(red)
    normalised_green, greenmin, greenmax = matrix_rgb_normalisation(green)

    zeros = np.zeros(normalised_green.shape)
    rgb = np.dstack((normalised_red, normalised_green, zeros))
    rgb = np.rot90(rgb)
    if save_figure != 'LargeImage':
        plt.imshow(rgb.astype('uint8'), origin='lower')
        tick_positions = np.arange(0, len(normalised_green), len(normalised_green) / 4)
        tick_labels = np.arange(0, len(normalised_green) / scale_factor,
                                len(normalised_green) / scale_factor / 4).round(decimals=2)
        plt.xticks(tick_positions, tick_labels)
        plt.yticks(tick_positions, tick_labels)
        plt.ylabel('y axis (mm)', size=16)
        plt.xlabel('x axis (mm)', size=16)
        plt.yticks(size=15)
        plt.xticks(size=15)
        plt.title('parID=' + str(parID), size=14)
        np.set_printoptions(precision=2)
        plt.text(1,1,'mCherry = [%r-%r]'%(np.around(redmin,2),np.around(redmax,2)),c='r')
        plt.text(1,5,'GPF = [%r-%r]'%(np.around(greenmin,2),np.around(greenmax,2)),c='g')
        plt.tight_layout()

        if save_figure == True:
            plt.savefig(path + '/%s_%s.jpeg' % (dimension, filename),dpi=2000)
            plt.close()
        else:
            plt.show()





    return rgb

def plot_redgreenblue_contrast(final_concentration, mm, mechanism, shape, filename, path, mask, parID=0, scale_factor=10, save_figure=False, dimension='2D'):
    green = final_concentration[-1]
    red = final_concentration[-2]
    blue = final_concentration[0]
    if mask != None:
        blue = blue*mask[:,:,-1]
    x_grid = np.linspace(0, mm, len(green))
    normalised_red, redmin, redmax = matrix_rgb_normalisation(red)
    normalised_green, greenmin, greenmax = matrix_rgb_normalisation(green)
    normalised_blue, bluemin, bluemax = matrix_rgb_normalisation(blue)

    zeros = np.zeros(normalised_green.shape)
    rgb = np.dstack((normalised_red, normalised_green, normalised_blue))
    rgb = np.rot90(rgb)
    if save_figure != 'LargeImage':
        plt.imshow(rgb.astype('uint8'), origin='lower')
        tick_positions = np.arange(0, len(normalised_green), len(normalised_green) / 4)
        tick_labels = np.arange(0, len(normalised_green) / scale_factor,
                                len(normalised_green) / scale_factor / 4).round(decimals=2)
        plt.xticks(tick_positions, tick_labels)
        plt.yticks(tick_positions, tick_labels)
        plt.ylabel('y axis (mm)', size=16)
        plt.xlabel('x axis (mm)', size=16)
        plt.yticks(size=15)
        plt.xticks(size=15)
        plt.title('parID=' + str(parID), size=14)
        np.set_printoptions(precision=2)
        plt.text(1,1,'mCherry = [%r-%r]'%(np.around(redmin,2),np.around(redmax,2)),c='r')
        plt.text(1,5,'GPF = [%r-%r]'%(np.around(greenmin,2),np.around(greenmax,2)),c='g')
        plt.text(1,10,'GPF = [%r-%r]'%(np.around(bluemin,2),np.around(bluemax,2)),c='g')
        plt.tight_layout()

        if save_figure == True:
            plt.savefig(path + '/%s_%s.jpeg' % (dimension, filename),dpi=2000)
            plt.close()
        else:
            plt.show()


    return rgb



#video
def redgreen_contrast_timeseries(records):
    rgb_timeseries = []
    simulation_time = len(records[0][0][0])
    for time in range (simulation_time):
        red_timeseries,green_timeseries = records[-2],records[-1]
        red = red_timeseries[:,:,time]
        green = green_timeseries[:,:,time]
        normalised_red = matrix_rgb_normalisation(red)[0]
        normalised_green = matrix_rgb_normalisation(green)[0]
        zeros = np.zeros(red.shape)
        rgb = np.dstack((normalised_red,normalised_green,zeros))
        rgb = np.rot90(rgb)
        rgb_timeseries.append(rgb)
    return rgb_timeseries

def show_rgbvideo(timeseries_unstacked,parID):
    time=0

    fig = plt.plot()
    rgb_timeseries=timeseries_unstacked # Read the numpy matrix with images in the rows
    im=plt.imshow(rgb_timeseries[0].astype('uint8'), origin= 'lower')
    for time in range(len(rgb_timeseries)):
        im.set_data(rgb_timeseries[time].astype('uint8'))
        plt.title(parID)
        plt.xlabel(time)
        plt.pause(0.01)
    plt.show()
