#~/opt/local/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import MaxNLocator

# Update the matplotlib configuration parameters:
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams.update({'font.size': 12,
                    'font.family': 'STIXGeneral',
                    'mathtext.fontset': 'stix'})

def setup_axes( j, label ):
    # definitions for plot size
    width, height = 0.75, 0.145     ## size of figure
    width_hist = 0.10              ## width of histogram
    wx, wy = 0.005, 0.0            ## space between plots
    left   = 0.10                  ## left margin
    bottom = 0.10 + j*(height + wy)## bottom margin
    # plot coordinates
    pos_axScatter = [left, bottom, width, height]
    pos_axHist    = [left + width + wx, bottom, width_hist, height]

    # start with a rectangular Figure
    plt.figure(2, figsize=(12, 10))

    axScatter = plt.axes(pos_axScatter) # walker plot
    axHist = plt.axes(pos_axHist)     # histogram of walker plot

    # set labels
    axScatter.set_xlabel('Iteration #' if j == 0 else ' ')
    axScatter.set_ylabel(label)

    axHist.set_xlabel('N' if j == 0 else ' ')
    axHist.yaxis.set_major_formatter(NullFormatter())         # no labels

    axScatter.set_xticks(np.arange(0,11000,2000))
    axScatter.yaxis.set_major_locator(MaxNLocator(3))

    # no tick labels for stacked plots
    if j > 0:
        axScatter.set_xticklabels([])
        axHist.set_xticklabels([])

    return (axScatter, axHist)

def plot_data(y, axScatter, axHist):
    # draw the scatter plot
    axScatter.scatter( np.arange(1,np.size(y)+1,1), y, s=3, alpha=0.5 )
    # xlimit for axScatter -- could make # of steps completed, but i like total steps better
    axScatter.set_xlim( [0, 11000] )

    # draw the histogram -- set nbins = 20 -- should be fine for walkers
    axHist.hist(y, bins=20, histtype='step', log=True, orientation='horizontal',
                color='black')
    # set limits for axHhist
    axHist.set_ylim( axScatter.get_ylim() )
    axHist.set_xlim( 1, axHist.get_xlim()[1] )

def main():
    # read in the .res log file -- should be grabbed from Mcmc call
    fin = 'log.res'
    data = np.genfromtxt(fin, names=True)

    fout = 'mcmcWalkers.png' #  output files for walkers

    for j in range( np.size(data.dtype.names) ):
        # set up the axes for the walker and the hist based on the .res file
        axScatter, axHist = setup_axes( j, data.dtype.names[j] )
        plot_data(data[data.dtype.names[j]], axScatter, axHist)

        plt.savefig(fout, format='png')

if __name__ == "__main__":
    main()
