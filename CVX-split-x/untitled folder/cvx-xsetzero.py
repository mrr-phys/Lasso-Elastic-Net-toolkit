# 

# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""

import matplotlib
import matplotlib.pyplot as plt
import cvxopt as cvx
from cvxopt import solvers
import numpy as np
import pylab
from time import time
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import *
from scipy.stats.stats import pearsonr
from numpy import linalg as la
#import xlrdmatplotlib.use('qt')




def svd(data, S=2):
     
    #calculate SVD
    U, s, Vt = linalg.svd( data )
    rowl=U.shape[0]
    V = Vt.T
    #take out columns we don't need
 #   Sig = np.eye(S)*s[:S]
    
    newU = U[:,:S]
  #  newV = V[:,:S]
    
    
    
  #  fig = plt.figure()
   # ax = fig.add_subplot(1,1,1)
    #fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True) 

    plt.xlabel('SVD1')
    plt.ylabel('SVD2')
    labels = ['point{0}'.format(i) for i in range(rowl)]
        
    for i in xrange(rowl):
        
        plt.scatter(newU[i,0],newU[i,1], color= 'blue')
    
    for label, x, y in zip(labels, newU[:,0], newU[:,1]):
        plt.annotate(label, xy = (x, y), xytext = (0, 0), textcoords = 'offset points')


       
    plt.show() 


def cvx_l1(H, y, i, lamda2):
    
    
    m,n=H.shape
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T

    stz=np.zeros((2*n,2*n))
    if i==n+1:
        stz=np.zeros((2*n,2*n))
    else:
        stz[i,i]=2
        stz[i+n,i+n]=2


    #x-split l1 approach  
    c = u(h([np.ones(n),np.ones(n)]).T)
    
    G1 = h([HdH+lamda2,-(HdH+lamda2)])
    G2 = h([-(HdH+lamda2),HdH+lamda2])
    G3 = h([-np.eye(2*n)])
    G4 = h([-np.eye(2*n)+stz])
    G = u(v([G1,G2,G3,G4]))    
    
    hh = u(h([np.dot(Hd,y),-np.dot(Hd,y),np.zeros(2*n).T,np.zeros(2*n).T]))
 


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

#Launching :

if __name__=='__main__':
    

    n = 200 # number of predictors
    k = 10 # number of sparsity
    h_iter = 1 # h-iteration over different random sensing matrix 
     
    lamda2=0
    #sigma2_n = 1./np.sqrt(100)
    
    m =50  # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    #x[:k] = np.sign(np.random.normal(0,1,k))
    
    #np.random.shuffle(x)
    xsparse = x.T

    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
    y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
      

    while (m>=30):
        
        subsamp= np.random.randint(H.shape[0],size=m)
        ysub= y[subsamp]
        Hsub= H[subsamp,:]
        
        x = cvx_l1(Hsub, ysub, n+1, lamda2)
        xcenter=x
        
        
        for i in range(k):
        
            xest = cvx_l1(Hsub, ysub,  i, lamda2)
            
            x = np.vstack((x,xest))
        
        
           
        #mean = x.mean(axis=0)
        #mean.shape = (1,mean.shape[0])
        #x = x-mean #takeout column mean
        #mean = x.mean(axis=1)
        #mean.shape = (mean.shape[0], 1)
        #x = x-mean  #takeout raw mean  
        
        x=x-xcenter
        svd(x,2)
        


            
       # with open("dat_pca.txt", "a") as myfile:
       #     myfile.write(str(float(m/200.))+' ' +str(pca.explained_variance_ratio_)+'\n')
        #    myfile.close()

        
            
        m-=2
    

    print "run is finished"