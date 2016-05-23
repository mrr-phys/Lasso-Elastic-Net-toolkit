#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rc('text', usetex=True)



with open("c:\\users\mohammad\desktop\zerotemp\datROC600.txt") as f:
    d0 = f.read()

d0 = d0.split('\n')

x0 = [row.split(' ')[0] for row in d0]
y0 = [row.split(' ')[1] for row in d0]

with open("c:\\users\mohammad\desktop\zerotemp\datROC601.txt") as f:
    d1 = f.read()

d1 = d1.split('\n')

x1 = [row.split(' ')[0] for row in d1]
y1 = [row.split(' ')[1] for row in d1]

with open("c:\\users\mohammad\desktop\zerotemp\datROC602.txt") as f:
    d2 = f.read()

d2 = d2.split('\n')

x2 = [row.split(' ')[0] for row in d2]
y2 = [row.split(' ')[1] for row in d2]



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title(r'$\alpha/\rho = 6.0, \sigma^{2} = 10^{-5}, \lambda_{1} = 1.0$')    
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')


ax.plot( x0, y0, marker = 'o', color='r',markersize=7, label=r'$\lambda_{2} = 0.0$')
ax.plot( x1, y1, marker = 's', color='b',markersize=7, label=r'$\lambda_{2} = 0.5$')
ax.plot( x2, y2, marker = '^', color='g',markersize=7, label=r'$\lambda_{2} = 1.0$')

leg = ax.legend(loc = 4)

# display the plot
plt.show()

dpi = 300
fig.savefig('p6.pdf', dpi=dpi)

