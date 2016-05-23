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
from scipy.linalg import toeplitz
from scipy.linalg import inv
from scipy.linalg import sqrtm



    
"The least square function"
def LS(x, H, y):

    return .5*sq_norm2(y - H.dot(x))
    
    
"Tthe gradient of the least square function"    
def LS_grad(x):
    
    return (HtH.dot(x)-Hty)

    
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
   

def sq_norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)
    
def mse(a):
    # return mean square error
    return (1./np.size(a))*sq_norm2(a)

    
if __name__=='__main__':
    
    n=20000
    k=50
    


    s_iter=4
    h_iter=8
    
    x_init = np.zeros(n)

    n_iter = 800

      # initialize number of observant to sparsity
    
    c0=0.
    for i in range(9):
        c_list = np.zeros(n)
        c0+=.1
        
        
        for i in range(n):
            c_list[i]=c0**i
            
        
        mean = np.zeros(n)
        cov = toeplitz(c_list)
        DDinv = (1+c0**2)/(1-c0**2)#np.trace(cov)*(np.trace(LA.inv(cov)))/n**2
        cov_half = sqrtm(cov)
        
        m = 0.02*DDinv*n
        print m, c0, DDinv


        while (m>=0.01*DDinv*n):               
 
       # mselist_cvx = []
            mselist_phd = []
        
            for i in range(s_iter):
            
                x = np.zeros(n)
                x[:k] = np.random.normal(0,1,k)
                np.random.shuffle(x)
                xsparse = x.T
                ##################


                for j in range(h_iter):
                    #producing sample sparse array
                
                    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n]).dot(cov_half)
                   # H = 1/np.sqrt(m)*np.random.multivariate_normal(mean,cov,m)

                    y = np.dot(H,xsparse)# + np.random.normal(0,sigma2_n,m) 
                    HtH = H.T.dot(H)
                    Hty = H.T.dot(y)

                    lmd_1 = LA.norm(Hty , np.inf)
                    lmd_2 = 0
                
            
                    xest_phd,lmbda = phd(x_init, lmbda=lmd_1, n_iter=n_iter) 
                  #  xest_phd = np.real(xest_phd)
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

            
            with open("datMSE1.txt", "a") as myfile:
                myfile.write(str(m/float(n))+' ' +str(np.sum(mselist_phd)/len(mselist_phd))+'\n')
                myfile.close()
            
            m-=30 # increment of undersampling step
        with open("datMSE1.txt", "a") as myfile:
            myfile.write(str(c0)+' '+str(DDinv)+'\n')
            myfile.close()


    print "run is finished"