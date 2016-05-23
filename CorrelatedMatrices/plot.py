# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:59:59 2016

@author: mohammad
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:09:08 2015

@author: mohammad
"""

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#rc('text', usetex=True)
#rcParams['text.latex.preamble']=[r"\usepackage{amsmath,bm}"]
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)
labels = []



xl=np.zeros(81)
yl=np.zeros(81)


with open('/Users/mohammad/Dropbox/Anirvan-Mohammad/Code-CD/phase.txt') as myfile:
    
    data = myfile.read()
        #data = data[:-1]
    data = data.split('\n')
        
    x = [row.split(' ')[0] for row in data]
    y = [row.split(' ')[2] for row in data]
        
   # x.pop()
   # y.pop()
        
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)


p1=plt.plot(x,y,'rs',linewidth=2.)
#p1=plt.plot(x,x/(4*.02*np.log(1/x)),'bs',linewidth=2.)

plt.xlim(0.16, .2)
plt.ylim(1, 2   )
#plt.xticks([0,0.04,0.08,0.12,0.16,0.2,0.24])
#plt.yticks([0,0.01,0.02,0.03,0.04])
#plt.xticks(x, labels, rotation='vertical')
#print j
#a, b=np.polyfit(x[:20]*3., yl[:20],1)
#print a
#plt.legend(loc=2, prop={'size':8})
plt.title(r'$N=2000,\;\rho = 0.02$')    
plt.xlabel(r'$\alpha$')
#plt.ylabel(r'$\overline{\Delta u(f)}=\frac{1}{N}\sum\limits_a^N\big(u_(f)-u(0)\big)$')
plt.ylabel(r'$Tr(D)Tr(Dinv)$')

     
plt.legend(loc=2)
plt.show()

plt.savefig('test.pdf')

