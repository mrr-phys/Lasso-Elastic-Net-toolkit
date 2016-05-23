

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
from sklearn.decomposition import PCA


def cvx_l1(H, y, f, i, lamda2):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T
    

    stzl=np.zeros((2*n,2*n))
    stzr=np.zeros(2*n)
    if f==0:
        stzl=np.zeros((2*n,2*n))
        stzr=np.zeros(2*n)
    if f>0:
        stzl[i,i]=2
        stzl[i+n,i+n]=2
        stzr[i]=f 
        
    if f<0:
        stzl[i,i]=2
        stzl[i+n,i+n]=2
        stzr[i+n]=f 


    #x-split l1 approach  
    c = u(h([np.ones(n),np.ones(n)]).T)
    
    G1 = h([HdH+lamda2,-(HdH+lamda2)])
    G2 = h([-(HdH+lamda2),HdH+lamda2])
    G3 = h([-np.eye(2*n)])
    G4 = h([-np.eye(2*n)+stzl])
    G = u(v([G1,G2,G3,G4]))    
    
    hh = u(h([np.dot(Hd,y),-np.dot(Hd,y),(np.zeros(2*n)-stzr).T,(np.zeros(2*n)+stzr).T]))
 
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

def l1(a):
    # return mean square error
    return np.sum(np.abs(a))

#Launching :

if __name__=='__main__':

    n = 200 # number of predictors
    k = 10 # number of sparsity
    h_iter = 1 # h-iteration over different random sensing matrix 
     
    lamda2=0
    #sigma2_n = 1./np.sqrt(100)
    
    m =50 # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    #x[:k] = np.sign(np.random.normal(0,1,k))
    
    #np.random.shuffle(x)
    xsparse = x.T
    
    

    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
    y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
      

    while (m>=30):
        
        
        x_f_list=np.zeros(20)
        
        for l in range(20):
            
            f=0
            
            subsamp= np.random.randint(H.shape[0],size=m)
            ysub= y[subsamp]
            Hsub= H[subsamp,:]
        
            x = cvx_l1(Hsub, ysub, 0, 0, lamda2)
 
            for j in range(20):
                
                
                xest = cvx_l1(Hsub, ysub, f, 0, lamda2) 
                x_f_list[j]+=l1(xest)-l1(x)
                f+=.01
       
        f=0
        for h in range(20):
            
            with open("dat_l1{}.txt".format(m), "a") as myfile:
                myfile.write(str(f)+' ' +str(x_f_list[h]/20.)+'\n')
                myfile.close()
            f+=.01
        
            
        m-=5
    

    print "run is finished"

