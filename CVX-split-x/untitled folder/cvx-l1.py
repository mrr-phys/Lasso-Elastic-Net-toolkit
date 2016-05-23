# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt import solvers
import numpy as np



def cvx_l1(H, y, lamda2):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T
   

    #x-split l1 approach  
    c = u(h([np.ones(n),np.ones(n)]).T)
    
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



#Launching :

if __name__=='__main__':

    n = 200 # number of predictors
    k = 10 # number of sparsity
    h_iter = 10 # g-iteration over different random sensing matrix 
    s_iter = 30
     
    lamda2=0
    sigma2_n = 1./np.sqrt(100)
    
    m = k  # initialize number of observant to sparsity

    while (m<n):
        
        mse_list = []
        
        for i in range(h_iter):
                
            H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
        
            for j in range(s_iter):
            #producing sample sparse array
                x = np.zeros(n)
                x[:k] = np.random.normal(0,1,k)
                np.random.shuffle(x)
                xsparse = x.T
                ##################
                #additive noise
                y = np.dot(H,xsparse) + np.random.normal(0,sigma2_n,m) 
            
                xest = cvx_l1(H, y, lamda2)  
                
                mse_list.append(mse(x-xest))
            


            # print mse_list

        with open("datn.txt", "a") as myfile:
            myfile.write(str(m/float(n))+' ' +str(np.log(np.sum(mse_list)/len(mse_list))/np.log(sigma2_n))+'\n')
            myfile.close()
            
        with open("datq.txt", "a") as myfile:
            myfile.write(str(m/float(n))+' ' +str(np.sum(mse_list)/len(mse_list))+'\n')
            myfile.close()
            
        m+=5 # increment of undersampling step


    print "run is finished"

