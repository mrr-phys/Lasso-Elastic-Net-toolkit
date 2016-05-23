# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 14:46:40 2016

@author: mohammad
"""
from numpy import linalg as LA
import numpy as np
from time import time
import matplotlib.pylab as plt
import l1ls as L
import cvxopt as cvx
from cvxopt import solvers
from scipy.linalg import inv


def cvx_lq(H, y, lmd_1,lmd_2):

    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    
    ydH = np.dot(y.T,H)/2.#/(2.*sig2)
    Hdy = np.dot(H.T,y)/2.#/(2.*sig2)
    HdH= np.dot(H.T,H)#/sig2
    q = u(h([-ydH-Hdy+lmd_1*np.ones(n),ydH+Hdy+lmd_1*np.ones(n)]).T)
    #P=u(lamda2*np.eye(2*n))
    hor1 = h([HdH+lmd_2*np.eye(n),-HdH-lmd_2*np.eye(n)])
   # hor2 = h([lamda2*np.eye(n),lamda2*np.eye(n)])
    P=u(v([hor1,-hor1]))
    

    G3 = h([-np.eye(2*n)])
    G = u(v([G3]))    
    
    hh = u(h([np.zeros(2*n).T]))
 
    sol = cvx.solvers.qp(P,q, G, hh)

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T


def cvx_lp(H, y):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
   

    H2 = u(h([H,-H]))
   

    #x-split l1 approach  
    c = u(h([np.ones(n),np.ones(n)]).T)
    G = u(-np.eye(2*n).T)    
    
    hh = u(np.zeros(2*n).T)

    sol = cvx.solvers.lp(c, G, hh, H2,u(y))

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T

def sq_norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)
    
"The least square function"
def LS(x, H, y):

    return .5*sq_norm2(y - H.dot(x))
    
    
"Tthe gradient of the least square function"    
def LS_grad(x):
    
    return (HtH.dot(x)-Hty)
    
    #return -H.T.dot(y - H.dot(beta))/m
    
#"Proximal operator of the l1 norm"       
#def prox_l1(beta, lmd=1.):

#    return np.sign(beta) * (np.abs(beta) - lmd) * (np.abs(beta) > lmd)

#" Proximal operator of the l2 norm"
#def prox_l2(beta, lmd=1.):

#    return 1. / (1 + lmd) * beta

    
"Proximal operator for Elastic Net"
def prox_enet(beta,  t=1.):

   
    prox_l1 = np.sign(beta) * (np.abs(beta) -  t) * (np.abs(beta) > t)
    return prox_l1 / (1. + lmd_2) 

"Proximal Gradiant Descent algorithm"
def phd(x_init, lmbda=1., n_iter=500):

    x = x_init.copy()
    for _ in range(n_iter):
        xt=x
        x = prox_enet(xt - LS_grad(xt),lmbda)
        lmbda = .92*lmbda
        
        if LA.norm(x,0)>2*k:
            xhat=xt
            break

    return xhat,lmbda
   

def mse(a):
    # return mean square error
    return (1./np.size(a))*sq_norm2(a)

    
if __name__=='__main__':
    
    n=2000
    k=40
    


    s_iter=1
    h_iter=30
    x_init = np.zeros(n)
    n_iter = 800

    m = 32  # initialize number of observant to sparsity

    while (m>=28):
        
        mselist_cvx = []
        mselist_phd = []
        
        for i in range(s_iter):
            
            x = np.zeros(n)
            x[:k] = np.random.normal(0,1,k)
            np.random.shuffle(x)
            xsparse = x.T
            ##################


            for j in range(h_iter):
            #producing sample sparse array

            
                mean = np.zeros(n)
                cov = np.random.exponential(scale=.1, size=[n,n])
                #cov = np.cov(xi)
                print LA.inv(cov)
                print cov.dot(LA.inv(cov))

                break
                
            
                
                #Hx = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
                H = 1/np.sqrt(m)*np.random.multivariate_normal(mean,cov,m)

            
                #print H.shape
                y = np.dot(H,xsparse)# + np.random.normal(0,sigma2_n,m) 
                HtH = H.T.dot(H)
                Hty = H.T.dot(y)

                lmd_1 = LA.norm(Hty , np.inf)
                lmd_2 = 0
                
            
                xest_phd,lmbda = phd(x_init, lmbda=lmd_1, n_iter=n_iter) 
              #  lmd_1=lmbda                
             #   xest_cvx = cvx_lq(H, y,lmd_1,lmd_2)
                
                            
               # plt.figure(figsize=(18, 5))
               # plt.suptitle("", fontsize=18)
               # plt.subplot(1, 2, 1)
               # plt.title("Original Signal")
               # plt.stem(-np.sort(-np.abs(xsparse))[:300])
               # plt.subplot(1, 2, 2)
               # plt.title("Homotopy_Proximal Gradiant Descent")
               # plt.stem(-np.sort(-np.abs(xest))[:300], color='red')  
                
                mselist_phd.append(mse(xsparse-xest_phd))
              #  mselist_cvx.append(mse(xsparse-xest_cvx))

            
        with open("datMSE.txt", "a") as myfile:
            myfile.write(str(m/float(n))+' ' +str(np.sum(mselist_phd)/len(mselist_phd))+'\n')
            myfile.close()
            
        m-=5 # increment of undersampling step


    print "run is finished"