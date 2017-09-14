# Given a UTR branch length feel and the name of a transcript, segments the UTR
# into smaller regions of distinct conservation and plots output 
# USAGE: plot_transitions.py INFILE TRANSCRIPT OUTFILE
# INFILE is a file of base wise 3` UTR branch length scores, provided in the repo 
# TRANSCRIPT is a transcript name from a line in INFILE
# OUTFIELD is the .pdf filename for the output segmentation plot

import sys

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.stats import ttest_ind

COLOR1 = '#1b9e77'
COLOR2 = '#d95f02'
#COLOR1 = colors.to_hex('m')
#COLOR2 = colors.to_hex('y')
COLOR1= 'm'
COLOR2= 'y'


INFILE, TRANSNAME, OUTFILE = sys.argv[1:]


def find_transitions(obs, tol=-15, window=50, min_peak_space=100):
    logpvals = []
    for i in range(2,len(obs)-2):
        # obtain right and left windows
        left = np.array(obs[max(0,i-window):i])
        right = np.array(obs[i:min(i+window, len(obs))])

        # calculate t-statistic and p val
        t, p = ttest_ind(left, right, equal_var=False)

        # record log p-value, scale by window size for truncated windows
        logpvals.append(np.log(p)*(min(len(left),len(right))/(window*1.0)))

    # obtain the highest peaks, as defined by the tolerance
    locs = [i+2 for (i,x) in enumerate(logpvals) if x <= tol]
    peaks = [x for x in logpvals if x <= tol]
    i = 0
    while i < (len(locs)-1):
        if ((locs[i+1] - locs[i]) < min_peak_space):
            if peaks[i+1] < peaks[i]:
                peaks.pop(i)
                locs.pop(i)
            else:
                peaks.pop(i+1)
                locs.pop(i+1)
        else:
            i += 1

    # return peak locations and probs
    return locs + [len(obs)], np.array(logpvals)

found = False
with open(INFILE, 'r') as infile:
    for line in infile:
        line = line.split(',')

        # extract the gene we want
        if line[0] == TRANSNAME:
            obs = [float(x) for x in line[1:]]
            peaks, logp = find_transitions(obs)

            # plot results
            #fig = plt.figure(figsize=(10,5))
            #ax = fig.add_subplot(1, 1, 1)
            fig, (ax1, ax2) = plt.subplots(2,sharex=True,figsize=(10,5))

            ax2.plot(np.arange(len(logp)) + 2, logp, color='grey')
            ax2.axhline(y=-15,ls='--',color='black')
            # plot segments, alternating between the two colors
            prev = 0
            for i,peak in enumerate(peaks):
                if i % 2 == 0:
                    color = COLOR1
                else:
                    color = COLOR2
                ax1.scatter(np.arange(prev, peak), obs[prev: peak], color=color, alpha=0.5)

                prev = peak

            # Eliminate upper and right axes
            ax1.spines['right'].set_color('none')
            ax1.spines['top'].set_color('none')
            ax2.spines['right'].set_color('none')
            ax2.spines['bottom'].set_color('none')
            # Show ticks in the left and lower axes only
            ax1.xaxis.set_ticks_position('bottom')
            ax1.yaxis.set_ticks_position('left')
            ax2.xaxis.set_ticks_position('top')
            ax2.yaxis.set_ticks_position('left')

            ax1.set_ylabel('Branch length score')
            ax2.set_ylabel('log(transition p-value)')
            ax1.set_ylim(bottom=0)
            plt.xlim(-30, len(obs) + 30)
            plt.xlabel('Position along 3\'UTR')
            plt.savefig(OUTFILE)
            plt.close()

            found = True
            break

if not found:
    print('{} not found'.format(TRANSNAME))
