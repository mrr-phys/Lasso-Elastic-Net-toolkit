
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt.solvers import qp
import numpy as np
import pylab
from time import time
import scipy
from scipy.linalg import inv
from scipy.linalg import block_diag



def cvx_l1(H, y, lamda1,lamda2, sig2):

    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix



#Linear programming#############################

    #x-split l1 approach  
  #  c = u(h([np.ones(n),np.ones(n)]).T)
    
  #  G1 = h([H,-H])
  #  G2 = h([-H,H])
  #  G3 = h([-np.eye(2*n)])
   # G = u(v([G1,G2,G3]))    
    
   # hh = u(h([y,-y,np.zeros(2*n).T]))
 


   # sol = cvx.solvers.lp(c, G, hh)

   # doublex = np.array(sol['x'][:2*n]).T[0]
   # xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
  #  return xestimate.T
    
    
#Quadratic programming#############################
    
    ydH = np.dot(y.T,H)/(2.*sig2)
    Hdy = np.dot(H.T,y)/(2.*sig2)
    HdH= np.dot(H.T,H)/sig2
    q = u(h([-ydH-Hdy+lamda1*np.ones(n),ydH+Hdy+lamda1*np.ones(n)]).T)
    #P=u(lamda2*np.eye(2*n))
    hor1 = h([HdH+lamda2*np.eye(n),-HdH-lamda2*np.eye(n)])
   # hor2 = h([lamda2*np.eye(n),lamda2*np.eye(n)])
    P=u(v([hor1,-hor1]))
    
    
    

    G3 = h([-np.eye(2*n)])
    G = u(v([G3]))    
    
    hh = u(h([np.zeros(2*n).T]))
 


    sol = qp(P,q, G, hh)

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T
    
    
def norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)

def mse(a):
    # return mean square error
    return (1./np.size(a))*norm2(a)



#Launching:

if __name__=='__main__':

    n = 1000 # number of predictor
    k = 1 # number of sparsity
    h_iter = 10 # g-iteration over different random sensing matrix 
    s_iter = 40
     
    lamda1=1e-4
    lamda2=0#.4e-4
    sig2=1.
    sigma2_n = 0#1./np.sqrt(100)

    m = 18  # initialize number of observant to sparsity

    
    Dhalf=(1-.4)*np.eye(n,n)+(.4/np.sqrt(1000))*np.ones((n,n))
    D = np.dot(Dhalf.T,Dhalf)

   

    while (m>=6):
        
        mse_hlist = []
        
        for i in range(h_iter):
                
            H = cvxopt.matrix(np.dot((1./np.sqrt(m))*np.random.randn(m,n),Dhalf))
            #H = cvxopt.matrix((1./np.sqrt(m))*np.random.randn(m,n))
            
            mse_slist = []
        
            for j in range(s_iter):
            #producing sample sparse array
                x = np.zeros(n)
                x[:k] = np.random.normal(0,1,k)
                np.random.shuffle(x)
                xsparse = x.T
                ##################
                #additive noise
                y = np.dot(H,xsparse)# + np.random.normal(0,sigma2_n,m) 
            
                xest = cvx_l1(H, y,lamda1,lamda2,sig2)  
                
                mse_slist.append(mse(x-xest))
                
            mse_hlist.append(np.sum(mse_slist)/s_iter)
            


            # print mse_list

        #with open("datn.txt", "a") as myfile:
       #     myfile.write(str(m/float(n))+' ' +str(np.log(np.sum(mse_list)/len(mse_list)))+'\n')
       #     myfile.close()
            
        with open("datMSE-CD4.txt", "a") as myfile:
            myfile.write(str(m/float(n))+' ' +str(np.sum(mse_hlist)/h_iter)+'\n')
            myfile.close()
            
        m-=2 # increment of undersampling step
        
        print np.trace(D)*np.trace(inv(D))/n**2.
        print np.trace(D)
        print np.trace(inv(D))


    print "run is finished"

