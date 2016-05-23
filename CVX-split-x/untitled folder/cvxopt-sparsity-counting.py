
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt import solvers
import numpy as np
import pylab
from time import time


def cvx_l1(H, y, p, lamda2):
    
    # p is the purturbing field
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T
   

    #x-split l1 approach  
    c = u(h([np.ones(n)-p,np.ones(n)+p]).T)
    
    G1 = h([HdH+lamda2,-(HdH+lamda2)])
    G2 = h([-(HdH+lamda2),HdH+lamda2])
    G3 = h([-np.eye(2*n)])
    G = u(v([G1,G2,G3]))    
    
    hh = u(h([np.dot(Hd,y),-np.dot(Hd,y),np.zeros(2*n).T]))
 


    sol = cvx.solvers.lp(c, G, hh)

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T
    
def sq_norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)
    
def mse(a):
    # return mean square error
    return (1./np.size(a))*sq_norm2(a)

def avg(a):
    # return mean square error
    return (1./np.size(a))*np.sum(a)
    
def counter(a, th):

    a = np.array(a)
    x = np.zeros(np.size(a))

    x = np.abs(a) >= th
    a[x.nonzero()] = 1

    x = np.abs(a) < th
    a[x.nonzero()] = 0

    return np.sum(a)
#Launching :

if __name__=='__main__':

    n = 200 # number of predictors
    k = 10 # number of sparsity
    h_iter = 1 # h-iteration over different random sensing matrix 
     
    lamda2=0
    #sigma2_n = 1./np.sqrt(100)
    
    m = 100  # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    #x[:k] = np.sign(np.random.normal(0,1,k))
    
    #np.random.shuffle(x)
    xsparse = x.T
    
    
    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
    y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
    m=47       

    while (m>=47):
        
        m-=1
        mse_list=[]
        
        for i in range(1):
            
            subsamp= np.random.randint(H.shape[0],size=m)
            ysub= y[subsamp]
            Hsub= H[subsamp,:]
            
            xest = cvx_l1(Hsub, ysub, p, lamda2) 
            
            mse_list.append(mse(xest-xsparse))
            
        #cnt_nonzero=counter((xest[:k]/h_iter-xsparse[:k])/xsparse[:k],1e-1)  #count number of non-zero component change
       # cnt_zero=counter(xest[k::]/h_iter-xsparse[k::],1e-1) #count number of zero component change
        #error=mse(xest/30-xsparse)
                
                
        with open("dat_countp.txt", "a") as myfile:
            myfile.write(str(m/float(n))+' ' +str(np.sum(mse_list)/30)+'\n')
            myfile.close()
            
            
        
            
       # with open("dat_countpp.txt", "a") as myfile:
      #      myfile.write(str(m/float(n))+' ' +str(xest[:k]/h_iter-xsparse[:k])+'\n')
      #      myfile.close()
    

    print "run is finished"

