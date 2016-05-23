# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 00:00:50 2015

@author: mohammad
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt import solvers
from cvxopt.solvers import qp
import numpy as np
import pylab
from time import time


def cvx_l1(H, y, p, i, lamda1,lamda2, sig2):

    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
    
    prtrb=np.zeros(n)
    prtrb[i]=p
    
    ydH = np.dot(y.T,H)/(2.*sig2)
    Hdy = np.dot(H.T,y)/(2.*sig2)
    HdH= np.dot(H.T,H)/sig2
    q = u(h([-ydH-Hdy+lamda1*np.ones(n)-prtrb,ydH+Hdy+lamda1*np.ones(n)+prtrb]).T)
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
    k = 30 # number of sparsity
    
    
    lamda1=1.
    lamda2=0#.4e-4
    sig2=1e-4
    #sigma2_n = 1./np.sqrt(100)
    
    m =40 # initialize number of observant to sparsity
    
    x = np.zeros(n)
    #x[:k] = np.random.normal(0,1,k)
    x[:k] = np.random.normal(0,1,k)
    
    np.random.shuffle(x)
    xsparse = x.T
    
   # print xsparse

    H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
    y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
      

    while (m>=40):
        
        
        avg_list=np.zeros(10)
        for l in range(6):
            
            p=0
            i=0
  
            subsamp= np.random.randint(H.shape[0],size=m)
            ysub= y[subsamp]
            Hsub= H[subsamp,:]
        
            xpzro = cvx_l1(Hsub, ysub, p, i, lamda1, lamda2, sig2)
            
            
           
            du_list=np.zeros(10)
            
            for i in range(n):
                
                p=0
                
                for j in range(10):
                    
                    xp = cvx_l1(Hsub, ysub, p, i, lamda1, lamda2, sig2)
                    du_list[j]+=xp[i]-xpzro[i]#avg(xp-xpzro)
                    #p_list[j]=p 
                    p+=.01
            avg_list+=du_list/200.
            
        p=0            
        for r in range(10):
            
                    
            with open("susc{}.txt".format(m), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
                myfile.write(str(p)+' ' +str(avg_list[r]/6.)+'\n')
                myfile.close()
            p+=.01
                        
               # with open("susc{}.txt".format(i), "a") as myfile:
            
               #     myfile.write(str(xsparse[i])+' ' +str(xsparse[i]))
               #     myfile.close()                        
                   
      
        m-=5
    

    print "run is finished"
