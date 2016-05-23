# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 14:46:40 2016

@author: mohammad
"""
from numpy import linalg as LA
import numpy as np
from time import time
import matplotlib.pylab as plt



"The least square function"
def LS(x, H, y):

    return .5*sq_norm2(y - H.dot(x))
    
    
"Tthe gradient of the least square function"    
def LS_grad(x):
    
    return (HtH.dot(x)-Hty)

    
"Proximal operator for Elastic Net"
def prox_enet(beta,  t1=1.,t2=1.):

   
    prox_l1 = np.sign(beta) * (np.abs(beta) -  t1) * (np.abs(beta) > t1)
    return prox_l1 / (1. + t2) 

"Proximal Gradiant Descent algorithm"
def phd(x_init, lmbda1=1.,lmbda2=1.,sigma2n=1., n_iter=500):

    x = x_init.copy()
    for _ in range(n_iter):
        xt=x
        x = prox_enet(xt - LS_grad(xt),lmbda1/sigma2n,lmbda2/sigma2n)
        lmbda1 = .92*lmbda1
       # lmbda2 = .82*lmbda2
        y = mse(xsparse-x)
        if lmbda1<1e-5:
            break
       #     break
            
        
        
        #if LA.norm(x,0)>2*k:
        #    xhat=xt
        #    break
    print x, lmbda1

    return x,lmbda1,lmbda2
    

def sq_norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)
    
def mse(a):
    # return mean square error
    return (1./np.size(a))*sq_norm2(a)




if __name__=='__main__':
    
    n=1000
    k=20
    m = 150
    sig2n=.01
    


    s_iter=1
    h_iter=1
    
    x_init = np.zeros(n)

    n_iter = 800
    mselist = []
    
    
    for i in range(s_iter):
        
        x = np.zeros(n)
        x[:k] = np.random.normal(0,1,k)
        np.random.shuffle(x)
        xsparse = x.T
        
        for j in range(h_iter):
            
            #H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n]).dot(cov_half)
            H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
            
            y = np.dot(H,xsparse)# + np.random.normal(0,sig2n,m) 
            HtH = H.T.dot(H)
            Hty = H.T.dot(y)
            
            lmd_1 = LA.norm(Hty , np.inf)
            lmd_2 = 0#10000#*lmd_1
                
            
            xest_phd,lmbda1,lmbda2 = phd(x_init, lmbda1=lmd_1,lmbda2=lmd_2, sigma2n = sig2n, n_iter=n_iter)
            
           # plt.figure(figsize=(18, 5))
          #  plt.suptitle("", fontsize=18)
          #  plt.subplot(1, 2, 1)
          #  plt.title("Original Signal")
           # plt.stem(-np.sort(-np.abs(xsparse))[:m])
           ## plt.subplot(1, 2, 2)
           # plt.title("Estimated Signal")

           # plt.stem(-np.sort(-np.abs(xest_phd))[:m])
            
            
            
            
            mselist.append(mse(xsparse-xest_phd))
              #  mselist_cvx.append(mse(xsparse-xest_cvx))

    print lmbda1, lmbda2        
    with open("datMSE.txt", "a") as myfile:
        myfile.write(str(m/float(n))+' ' +str(np.sum(mselist)/len(mselist))+'\n')
        myfile.close()
            



    print "run is finished"



