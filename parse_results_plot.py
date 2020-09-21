#!/usr/bin/env python3
from pyparsing import *
import numpy as np
import sys

to_int = pyparsing_common.integer

to_float = pyparsing_common.fnumber

env_pattern = Suppress(Literal("env") + ":") + to_int

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

import matplotlib as mpl
import matplotlib.pyplot as plt

size = (a_arr[0] + s_arr[0]).sum()

plt.rcParams.update({'font.family':'serif',
                     'font.sans-serif':'Palatino',
                     'text.usetex':True})

plt.figure(figsize=(10,5))
ax1 = plt.gca()

plt.plot(range(len(a_arr)),a_arr)
ax1.set_prop_cycle(mpl.rcParams['axes.prop_cycle'])
plt.plot(range(len(s_arr)),s_arr,ls='dotted')
# Autoscale pads, but set_ylim doesn't.  So I add invisible points.
plt.plot([0,0],[0,size],'y')[0].set_visible(False)

labels = ["{:2d}a".format(i) for i in range(a_arr.shape[1])] + \
        ["{:2d}s".format(i) for i in range(s_arr.shape[1])]
plt.legend(labels)

ax2 = plt.twinx()
ax2.plot(range(len(s_arr)),envs,'go',mfc='white')
# See above comment
ax2.plot([0,0],[0,s_arr.shape[1]-1])[0].set_visible(False)
ax2.grid()
ax1.grid(axis='x')

plt.savefig("out.pdf")
