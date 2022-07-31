#!/usr/bin/env python3
from pyparsing import *
import numpy as np
import sys

#
# Parse textual output and make a plot
#

to_int = pyparsing_common.integer

to_float = pyparsing_common.fnumber

env_pattern = Suppress(Literal("env") + ":") + to_float

def tally_pattern(label):
    return Suppress(label + ":") + OneOrMore(Suppress(to_int + ":") + to_int)

envs = []
line = next(sys.stdin)
envs += [env_pattern.parseString(line,parseAll=True)[0]]
line = next(sys.stdin)
col = tally_pattern("no_sex").parseString(line,parseAll=True).asList()
a_arr = np.array(col)
line = next(sys.stdin)
col = tally_pattern("sex").parseString(line,parseAll=True).asList()
s_arr = np.array(col)

for line in sys.stdin:
    envs += [env_pattern.parseString(line,parseAll=True)[0]]
    line = next(sys.stdin)
    col = tally_pattern("no_sex").parseString(line,parseAll=True).asList()
    a_arr = np.vstack((a_arr,col))
    line = next(sys.stdin)
    col = tally_pattern("sex").parseString(line,parseAll=True).asList()
    s_arr = np.vstack((s_arr,col))

# Binning
nloci = a_arr.shape[1] - 1
a_arr_old = a_arr
s_arr_old = s_arr
nbins = nloci + 1
if nloci > 9:
    nbins = 10
    lendpts = np.linspace(0,nloci+1,nbins+1)
    bins = []
    i = 0
    for j,x in enumerate(lendpts):
        while i < x:
            i += 1
        bins += [i]
    a_arr_old = a_arr
    s_arr_old = s_arr
    a_arr = np.empty((len(a_arr),len(bins)-1))
    s_arr = np.empty((len(s_arr),len(bins)-1))
    for i in range(len(bins)-1):
        a_arr[:,i] = a_arr_old[:,bins[i]:bins[i+1]].sum(1)
        s_arr[:,i] = s_arr_old[:,bins[i]:bins[i+1]].sum(1)

    
    labels = ["{:d}-{:d}a".format(bins[i],bins[i+1]-1) if bins[i] != bins[i+1]-1 else "{:d}a".format(bins[i]) for i in range(a_arr.shape[1])] + \
            ["{:d}-{:d}s".format(bins[i],bins[i+1]-1) if bins[i] != bins[i+1]-1 else "{:d}s".format(bins[i]) for i in range(s_arr.shape[1])] 
else:
    labels = ["{:2d}a".format(i) for i in range(a_arr.shape[1])] + \
            ["{:2d}s".format(i) for i in range(s_arr.shape[1])]

import matplotlib as mpl
import matplotlib.pyplot as plt
import getopt

ugly = False
if "--ugly" in sys.argv[1:]:
    ugly = True

size = (a_arr[0] + s_arr[0]).sum()

plt.rcParams.update({'font.family':'serif',
                     'font.sans-serif':'Palatino',
                     'text.usetex':True,
                     'axes.prop_cycle':plt.rcParams['axes.prop_cycle'][:nbins]})

plt.figure(figsize=(10,5),tight_layout=True)
ax1 = plt.subplot(1,6,(1,5))
ax2 = plt.subplot(1,6,6,frame_on=False,xticks=[],yticks=[])

lines = ax1.plot(range(len(a_arr)),a_arr)
ax1.set_prop_cycle(plt.rcParams['axes.prop_cycle'])
lines += ax1.plot(range(len(s_arr)),s_arr,ls='dashed' if ugly else 'dotted')
# Autoscale does padding, but set_ylim doesn't.  So I add invisible points.
ax1.plot([0,0],[0,size],'y')[0].set_visible(False)

leg =ax2.legend(lines,labels,loc="center")
if ugly:
    for line in leg.get_lines():
        line.set_linewidth(1.5*line.get_linewidth())

ax1_r = ax1.twinx()
ax1_r.plot(range(len(s_arr_old)),envs,'go',mfc='white')
# See above comment
ax1_r.plot([0,0],[0,s_arr_old.shape[1]-1])[0].set_visible(False)
ax1_r.grid()
ax1.grid(axis='x')

plt.savefig("out.pdf")
print("Plotted in out.pdf",file=sys.stderr)
