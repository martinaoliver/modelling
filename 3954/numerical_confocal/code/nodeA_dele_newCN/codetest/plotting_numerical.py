import numpy as np
# import matplotlib
# matplotlib.use('TKAgg')
import matplotlib.pyplot as plt


def matrix_rgb_normalisation(matrix):
    row_n = 0
    NewMatrix = np.zeros(matrix.shape)

    OldMin = np.min(matrix)
    OldMax = np.amax(matrix)
    NewMin = 0
    NewMax = 255
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)

    for row in matrix:
        column_n = 0
        for value in row:
            NewMatrix[column_n, row_n] = int((((value - OldMin) * NewRange) / OldRange) + NewMin)
            column_n += 1
        row_n += 1
    return NewMatrix, OldMin, OldMax


def plot_redgreen_contrast(final_concentration, mm, mechanism, shape, filename, path, parID=0, scale_factor=10, save_figure=False, dimension='2D'):
    green = final_concentration[-1]
    print(green)
    red = final_concentration[-2]
    var_red,var_green = [(np.amax(red) - np.amin(red)), (np.amax(green) - np.amin(green))]
    print(var_green)
    x_grid = np.linspace(0, mm, len(green))

    normalised_red, redmin, redmax = matrix_rgb_normalisation(red)
    normalised_green, greenmin, greenmax = matrix_rgb_normalisation(green)
    zeros = np.zeros(normalised_green.shape)
    rgb = np.dstack((normalised_red, normalised_green, zeros))
    rgb = np.rot90(rgb)
    # if save_figure != 'results':
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
