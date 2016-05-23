


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


def cvx_l1(H, y, p, i, lamda2):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T
    
    prtrb=np.zeros(n)
    prtrb[i]=p


    #x-split l1 approach  
    c = u(h([np.ones(n)-prtrb,np.ones(n)+prtrb]).T)
    
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

#Launching :

if __name__=='__main__':

    n = 200 # number of predictors
    k = 10 # number of sparsity
    h_iter = 1 # h-iteration over different random sensing matrix 
     
    lamda2=0
    #sigma2_n = 1./np.sqrt(100)
    
    m =150  # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    #x[:k] = np.sign(np.random.normal(0,1,k))
    
    #np.random.shuffle(x)
    xsparse = x.T
    
    

    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
    y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
      

    while (m>=150):
        
        p = -100
        subsamp= np.random.randint(H.shape[0],size=m)
        ysub= y[subsamp]
        Hsub= H[subsamp,:]
        
        x = cvx_l1(Hsub, ysub, 0, 0, lamda2)
        
        
        for i in range(k):
        
            xest = cvx_l1(Hsub, ysub, p, i, lamda2)
            
            x = np.vstack((x[:k],xest[:k]))
        
        
           
        mean = x.mean(axis=0)
        mean.shape = (1,mean.shape[0])
        x = x-mean #takeout column mean
        mean = x.mean(axis=1)
        mean.shape = (mean.shape[0], 1)
        x = x-mean  #takeout raw mean      
        pca = PCA(n_components=2)
        pca.fit(x)
        PCA(copy=True, n_components=2, whiten=False)
        pca_score = pca.explained_variance_ratio_
        V = pca.components_
        
        x_pca_axis, y_pca_axis = V.T * pca_score / pca_score.min()
        
        for i in xrange(11):
        ax.scatter(x_pca_axis,y_pca_axis, color= 'gray')
 
        show()
            
        with open("dat_pca.txt", "a") as myfile:
            myfile.write(str(float(m/200.))+' ' +str(pca.explained_variance_ratio_)+'\n')
            myfile.close()

        
            
        m-=5
    

    print "run is finished"

