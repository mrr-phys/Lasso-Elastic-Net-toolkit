# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:20:35 2015

@author: mohammad
"""

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
import time


def cvx_l1(H, y, xinit, f_a, a):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
   
   
    
    f_prtrb=np.zeros(n)
    f_prtrb[a]=f_a
    H2 = u(h([H,-H]))
   

    #x-split l1 approach  
    c = u(h([(np.ones(n)-f_prtrb),(np.ones(n)+f_prtrb)]).T)
    G = u(-np.eye(2*n).T)    
    
    hh = u(np.zeros(2*n).T)
    
  #  y = u(h([np.dot(Hd,y),-np.dot(Hd,y)
 


    sol = cvx.solvers.lp(c, G, hh, H2,u(y))
    
    return sol


    
    

    
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

    n = 500 # number of predictors
    k = 100 # number of sparsity
    

    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
    
    m =125 # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    
    np.random.shuffle(x)
    xsparse = x.T
    

      

    while (m>=125):
        
        
        for l in range(20):
            
            p=0
            i=0
            
            np.random.seed()
            H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
            y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
            
            xp = np.zeros(n)
            xm = np.zeros(n)
            ap=xsparse>0
            xp[ap]=xsparse[ap]
            am=xsparse<0
            xm[am]=xsparse[am]
            xinit=v([xp,xm])
      
            
  
          #  subsamp= np.random.randint(H.shape[0],size=m)
          #  ysub= y[subsamp]
          #  Hsub= H[subsamp,:]
        
            sol = cvx_l1(H, y, xinit, p, i)
            doublex_zr = np.array(sol['x'][:2*n]).T[0]
            xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex_zr)
            xpzro= xestimate.T
            
           # print mse(xpzro-xsparse)
            
          #  np.savetxt('H.txt', Hsub) 
          #  np.savetxt('b.txt', ysub) 
           
            du_list=np.zeros(20)
            
            for i in range(n):
                
                p=0
                doublex=doublex_zr
                
                for j in range(20):
                        
                        xinit = doublex
                        sol = cvx_l1(H, y, xinit, p, i)
                        doublex = np.array(sol['x'][:2*n]).T[0]
                        xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
                        xp = xestimate.T
                        
                        du_list[j]=xp[i]-xpzro[i]#avg(xp-xpzro)
                        p+=.1

                    
                p=0            
                for r in range(20):
            
                    
                    with open("sucf{}_{}.txt".format(l,i), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
                        myfile.write(str(p)+' ' +str(du_list[r])+'\n')
                        myfile.close()
                        p+=.1
                        
                with open("sucf{}_{}.txt".format(l,i), "a") as myfile:
                    myfile.write(str(xsparse[i])+' ' +str(xsparse[i]))
                    myfile.close()                        


                   
      
        m-=40
        


    

    print "run is finished"
