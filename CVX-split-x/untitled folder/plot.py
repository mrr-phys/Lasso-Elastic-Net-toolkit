#!/usr/bin/python

import numpy as np
from pylab import *


rc('text', usetex=True)
rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]




with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_suscep35.txt') as f:
    dmath2 = f.read()

dmath2 = dmath2.split('\n')

xmath2 = [row.split(' ')[0] for row in dmath2]
ymath2 = [row.split(' ')[1] for row in dmath2]

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_suscep40.txt') as f:
    dmath3 = f.read()

dmath3 = dmath3.split('\n')

xmath3 = [row.split(' ')[0] for row in dmath3]
ymath3 = [row.split(' ')[1] for row in dmath3]

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_suscep45.txt') as f:
    dmath4 = f.read()

dmath4 = dmath4.split('\n')

xmath4 = [row.split(' ')[0] for row in dmath4]
ymath4 = [row.split(' ')[1] for row in dmath4]

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_suscep50.txt') as f:
    dmath5 = f.read()

dmath5 = dmath5.split('\n')

xmath5 = [row.split(' ')[0] for row in dmath5]
ymath5 = [row.split(' ')[1] for row in dmath5]




##########################

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_susc45.txt') as f:
    d3 = f.read()

d3 = d3.split('\n')

x3 = [row.split(' ')[0] for row in d3]
y3 = [row.split(' ')[1] for row in d3]

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_susc40.txt') as f:
    d4 = f.read()

d4 = d4.split('\n')

x4 = [row.split(' ')[0] for row in d4]
y4 = [row.split(' ')[1] for row in d4]

with open('C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\cvx-split-x\dat_susc35.txt') as f:
    d5 = f.read()

d5 = d5.split('\n')

x5 = [row.split(' ')[0] for row in d5]
y5 = [row.split(' ')[1] for row in d5]

#xlim(0,1.5)
plot( xmath2,ymath2, 'r', label=r'$\boldsymbol{\alpha}=.175$')
plot( xmath3,ymath3, 'b', label=r'$\boldsymbol{\alpha}=.2$')
plot( xmath4,ymath4, 'g', label=r'$\boldsymbol{\alpha}=.225$')
plot( xmath5,ymath5, 'k', label=r'$\boldsymbol{\alpha}=.25$')
title(r'$\rho = 0.05$')    
xlabel(r'$p\,\,(\mathrm{external \, field \,at\, each\, site})$')
ylabel(r'$\boldsymbol{\bar{x}}$')
xlim(0,1.501)
ylim(0,.25)
a = axes([.15, .5, .2, .3], axisbg='y')
(m3,b3)=np.polyfit(x3,y3,1)
yp3=polyval([m3,b3],x3)

plot( x5,y5, 'r', label=r'$\boldsymbol{\alpha}=.175$')
plot( x4,y4, 'b', label=r'$\boldsymbol{\alpha}=.2$')
plot( x3,y3, 'g',label=r'$\boldsymbol{\alpha}=.225$')
plot( x3,yp3, '--k')
xlim(0,.1)
ylim(0,.025)
setp(a, xticks=[], yticks=[])

# display the plot

show()

#plt.savefig('plot.pdf')
#plt.close(fig)

