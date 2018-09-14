#-*-coding:utf-8-*-
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import CirclePolygon, Circle
import pickle
import argparse
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# get_ipython().magic(u'matplotlib inline')
def setBeamPlotShift(length=1.7, bias=0.2):
    #length, bias = 1.7, 0.2 # best 1.7, 0.2
    #figure_size = 5 # 
    shift_x = np.array([length, 0])
    bias_y = np.array([0, bias])
    shift_60 = np.array([length*np.cos(np.pi/3), length*np.sin(np.pi/3)]) - bias_y
    beam_shift = {1:[0, 0],
                  2:-shift_x,
                  3:-shift_60,
                  4:-shift_60 + shift_x,
                  5:shift_x,
                  6:shift_60,
                  7:shift_60 - shift_x,
                  8:-shift_x * 2,
                  9:-shift_60 - shift_x,
                  10:-shift_60 * 2,
                  11:-shift_60 * 2 + shift_x,
                  12:-shift_60 * 2 + shift_x * 2,
                  13:-shift_60 + shift_x*2,
                  14:shift_x * 2,
                  15:shift_60 + shift_x,
                  16:shift_60 * 2,
                  17:shift_60 * 2 - shift_x,
                  18:shift_60 * 2 - shift_x * 2,
                  19:shift_60-shift_x*2
                 }
    return beam_shift

def hexagonPlot(axin, x, y, centerxy=(0,0), radius=90, title='beam_0'): # best: radius=90
    new_ax = inset_axes(axins,
                        width="100%",
                        height="100%",
                        loc=10,
                        bbox_to_anchor=(centerxy[0], centerxy[1], 1, 1),
                        bbox_transform=axins.transAxes,
                        borderpad=0,
                       )
    hexagon = CirclePolygon((0, 0), radius, resolution=6, 
                            ec="black", fc='blue', alpha=0.1,
                            lw=2, label='hexagon {}'.format(centerxy))
    ada = AnchoredDrawingArea(0.0, 0.0, 0, 0, loc=10, pad=0, frameon=False)
    ada.da.add_artist(hexagon)
    new_ax.add_artist(ada)
    new_ax.plot(x, y, 'r-', linewidth=2)
    new_ax.set_title(title)
    return new_ax

def hexagonPlotInit(**kw):
    fig, ax = plt.subplots(**kw)
    axins = zoomed_inset_axes(ax, 1/figure_size, loc=10) # best 0.2
    return axins, fig, ax

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot 19 beams in hexagons")
    parser.add_argument('-f',type=str, help='directory where to load data from ... ')
    parser.add_argument('-png', type=str, help='filepath where to save the png')
    args = parser.parse_args()
    filein = args.f
    png_name = args.png

    filelist = os.listdir(filein)
    assert os.path.isdir('hexplotcache')==False,"hexplotcache exist."
    os.makedirs('hexplotcache') 
    assert len(filelist)==38, "{} should have exact 38 fits files.".format(filein)
    #save data to a cache directory
    to_plot = dict()
    for filename in filelist:        
        ffrom = os.path.join(filein, filename)
        fto = os.path.join('hexplotcache', filename[:-4]+'pkl')
        os.system("python load_file_to_pkl.py -f {} -t {}".format(ffrom, fto))
        # print some info
        print("Load: {}".format(filename))
        
        prefix, mxx, cache, beam_id, beam_pole, date, index = filename.split('_')
        assert cache=='multibeam'         # "__"
        assert 1<=int(beam_id)<=19 # beam id
        assert 0<=int(beam_pole)<=1  # pole
        to_plot.setdefault(int(beam_id), []).append((filename, int(beam_pole)))

    # Adapto python 2.7
    radius = 110 # hexagon radius e.g. 60
    fsize = 20 # figure size e.g. 12
    beam_shift = setBeamPlotShift(length=1.8, bias=-0.14) # adjust the distance between hexagons e.g 1.65, -0.05
    coef = 8 # Circle coef. e.g 8
    
    new_ax = dict()
    # Plot lines in hexagon
    fig, ax = plt.subplots(figsize=(fsize,fsize))
    axins = zoomed_inset_axes(ax, 0.1, loc=10)
    for beam_id, filenames in to_plot.items(): 
        
        for name in filenames:
            fullpath = os.path.join('hexplotcache', name[0][:-4]+'pkl')
            file = open(fullpath, 'rb')
            data = pickle.load(file)                      
            x = data['rfreq']
            y = data['mean_power']
            file.close()
            
            beam_pole = name[1]
            c = 'r-' if beam_pole else 'b-'
            centerxy=beam_shift[beam_id]
            
            if not new_ax.has_key(beam_id):
                new_ax[beam_id] = hexagonPlot(axins, x, y, centerxy=tuple(beam_shift[beam_id]),
                                              radius=radius, title='beam_{}'.format(beam_id))
            new_ax[beam_id].semilogy(x, y, c, linewidth=2, label='pole_{}'.format(beam_pole)) #  log view
        new_ax[beam_id].legend(loc="upper left",fontsize="x-small")
        
    os.system("rm -rf hexplotcache")
    # Big cirle
    hexagon = Circle((0, 0), radius*5, 
                     ec="black", fc='blue', alpha=0.1,
                     lw=2)
    ada = AnchoredDrawingArea(radius*coef, radius*coef, radius*coef/2, radius*coef/2, loc=10, pad=0, frameon=False)
    ada.da.add_artist(hexagon)
    ax.add_artist(ada)
    ax.axis('off')
    axins.axis('off')
    ax.set_title('Load files from dir:{}'.format(filein))
    
    plt.savefig(png_name)
