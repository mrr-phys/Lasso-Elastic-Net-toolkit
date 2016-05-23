# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt import solvers
import numpy as np
from scipy.linalg import block_diag
import pylab
from time import time



def solver(m,n,s,lamda2):
    print "Number of observations:"+str(m)
    print "Signal size:"+str(n)
    print "Sparsity:"+str(s)
    print "lambda2:"+str(lamda2)

    hor = np.hstack
    ver = np.vstack
    cbm = cvx.base.matrix
    
    #Construction of s-sparse vector x:
    x = np.zeros(n)
    x[:s] = 10*np.random.normal(0,1,s)
    np.random.shuffle(x)
    xsparse = x.T

    #Construction of matrix H and observation vector y :

    H1 = 1/np.sqrt(m/2)*np.random.normal(0,1,[m/2,n/2])
    H2 = 1/np.sqrt(m/2)*np.random.normal(0,1,[m/2,n/2])
    H = block_diag(H1,H2)
    y = np.dot(H,xsparse)
    
    HdH=np.dot(H.T,H)
    Hd= H.T
   

    #x-split l1  
    c = cbm(hor([np.ones(n),np.ones(n)]).T)
    
    G1 = hor([HdH+lamda2,-(HdH+lamda2)])
    G2 = hor([-(HdH+lamda2),HdH+lamda2])
    G3 = hor([-np.eye(2*n)])
    G = cbm(ver([G1,G2,G3]))    
    
    hh = cbm(hor([np.dot(Hd,y),-np.dot(Hd,y),np.zeros(2*n).T]))
 

    ini_time = time()
    sol = cvx.solvers.lp(c, G, hh)
    fin_time = time()-ini_time
    print "Computation time :"+str(fin_time)
    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(cbm(hor([np.eye(n),-np.eye(n)])),doublex)

    #Figures
    pylab.subplot(211)
    pylab.plot(xsparse)
    pylab.subplot(212)
    pylab.plot(xestimate)



    pylab.show()

#Launching :
solver(60,200,10,0.0)

