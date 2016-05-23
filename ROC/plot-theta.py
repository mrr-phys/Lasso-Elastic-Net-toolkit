#!/usr/bin/python

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rc('text', usetex=True)

with open("c:\\users\mohammad\desktop\zerotemp\TPR600.txt") as f:
    dq = f.read()

dq = dq.split('\n')

xq = [row.split(' ')[0] for row in dq]
yq = [row.split(' ')[1] for row in dq]

with open("c:\\users\mohammad\desktop\zerotemp\TNR600.txt") as f:
    dl0 = f.read()

dl0 = dl0.split('\n')

xl0 = [row.split(' ')[0] for row in dl0]
yl0 = [row.split(' ')[1] for row in dl0]



if 1:

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()


        
    #par2.axis["right"].toggle(all=True)

    host.set_xlim(0, 1.5)
    host.set_ylim(0, .9)

    host.set_title(r'$\alpha/\rho = 6.0, (\lambda_{1} = 1.0, \lambda_{2} = 0.0), (\mu = 0.0, \sigma^{2} = 10^{-5}), $')

    host.set_xlabel('criterion value')
    host.set_ylabel("True Positive Rate")
    par1.set_ylabel("True Negative Rate")


    p1, = host.plot(xq, yq, marker = 'o', color='b')
    p2, = par1.plot(xl0, yl0, marker = 's',color='r')

    par1.set_ylim(.7, 1.)

    host.legend()

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())

    plt.draw()
    plt.show()
  
