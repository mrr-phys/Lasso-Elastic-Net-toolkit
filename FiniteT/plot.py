#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rc('text', usetex=True)

with open("/Users/mohammad/Dropbox/convex\ optimization/BackgroundLiterature/Code/AlphaTau2.txt") as f:
    dmath = f.read()

dmath = dmath.split('\n')

xmath = [row.split(' ')[0] for row in dmath]
ymath = [row.split(' ')[1] for row in dmath]

with open("/Users/mohammad/Dropbox/convex\ optimization/BackgroundLiterature/Code/AlphaTau2.txt") as f:
    drep = f.read()

drep = drep.split('\n')

xrep = [row.split(' ')[0] for row in drep]
yrep = [row.split(' ')[1] for row in drep]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title(r'$\rho = 0.1$')    
ax.set_xlabel(r'\tau')
ax.set_ylabel(r'\alpha')

ax.plot(xmath, ymath, 'r', label ='CVX(ADMM)')
ax.plot( xrep, yrep, 'b',label = 'Code')

leg = ax.legend()

# display the plot
plt.show()

dpi = 300
fig.savefig('plot.pdf', dpi=dpi)

