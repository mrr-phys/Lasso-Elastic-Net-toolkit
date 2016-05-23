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


# Optimizing \lambda*\sigma^2(x_1+x_2) s.t. y=Hx_0 with CVX LP Solver 


def cvx_lp1(H, y, f_a, a):
    
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

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T
    
    
def cvx_lp2(H, y, f_a, a):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix

    HdH=np.dot(H.T,H)
    Hd= H.T
    
    f_prtrb=np.zeros(n)
    f_prtrb[a]=f_a
   

    #x-split l1 approach  
    c = u(h([np.ones(n)-f_prtrb,np.ones(n)+f_prtrb]).T)
    
    G1 = h([HdH,-HdH])
    G2 = h([-HdH,HdH])
    G3 = h([-np.eye(2*n)])
    G = u(v([G1,G2,G3]))    
    
    hh = u(h([np.dot(Hd,y),-np.dot(Hd,y),np.zeros(2*n).T]))
 


    sol = cvx.solvers.lp(c, G, hh)

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T
    
    
# Optimizing \frac{1}{\sigma^2}(y-Hx)^2+ \lambda(x_1+x_2)- f_a (x_1-x_2)_a
#    +\frac{1}{2}\lambda2(x_1+x_2)^2 s.t. x_1,x_2>=0 with CVX QP Solver 
def cvx_lq(H, y, f_a, a, lamda1,lamda2, sig2):

    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
    
    f_prtrb=np.zeros(n)
    f_prtrb[a]=f_a
    
    ydH = np.dot(y.T,H)/(2.*sig2)
    Hdy = np.dot(H.T,y)/(2.*sig2)
    HdH= np.dot(H.T,H)/sig2
    q = u(h([-ydH-Hdy+lamda1*np.ones(n)-f_prtrb,ydH+Hdy+lamda1*np.ones(n)+f_prtrb]).T)
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

    n = 300 # number of predictors
    k = 60 # number of sparsity
    
    
    lamda1=1.
    lamda2=0#.4e-4
    sigma2=1.
    #sigma2_n = 1./np.sqrt(100)
    
    m = 210 # initialize number of observant 
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    np.random.shuffle(x)
    xsparse = x.T
    
  #  H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
  #  y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
      

    while (m>=90):
        
        
        for l in range(1):  # \ell is the number of different H instances
            

            
            H = 1./np.sqrt(m)*np.random.normal(0,1,[m,n])
            y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
            
  
         #   subsamp= np.random.randint(H.shape[0],size=m)
         #   ysub= y[subsamp]
         #   Hsub= H[subsamp,:]
  
            f_a=0   # initial size of the perturbation f_a set to zero
            a=0    # initialize the starting site a
        
           # x_f_zro2 = cvx_lp2(H, y, f_a, a)
            x_f_zro1 = cvx_lp1(H, y, f_a, a)
            
          #  np.savetxt('H.txt', Hsub) 
          #  np.savetxt('b.txt', ysub) 
           
            deltau_list1=np.zeros(21)  #makes an array of size 10 for 10 different value of f_a
          #  deltau_list2=np.zeros(10) 
        

         #   start = time.time()


  

            for a in range(n):
                f_a = 0
                  
                for j in range(21):
                        
                 #       x_f_nzro2 = cvx_lp2(H, y, f_a, a)
                        x_f_nzro1 = cvx_lp1(H, y, f_a, a)
                        deltau_list1[j]=x_f_nzro1[a]-x_f_zro1[a]#avg(xp-xpzro)
                 #       deltau_list2[j]=x_f_nzro2[a]-x_f_zro2[a]
                        f_a+=.1

                f_a=0
                for r in range(21):

                    with open("suscf{}_{}.txt".format(m,a), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
                        myfile.write(str(f_a)+' ' +str(deltau_list1[r])+'\n')
                        myfile.close()
                        f_a+=.1
                        
                with open("suscf{}_{}.txt".format(m,a), "a") as myfile:   # Added to recognize the contribution from zero or non-zero signal
                        myfile.write(str(xsparse[a])+' ' +str(xsparse[a]))
                        myfile.close()      
                        
                        
                        
         #       f_a=0
         #       for r in range(10):

          #          with open("suscc{}_{}.txt".format(l,a), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
             #           myfile.write(str(f_a)+' ' +str(deltau_list2[r])+'\n')
             #           myfile.close()
             #           f_a+=.01
                        
             #   with open("suscc{}_{}.txt".format(l,a), "a") as myfile:   # Added to recognize the contribution from zero or non-zero signal
               #         myfile.write(str(xsparse[a])+' ' +str(xsparse[a]))
               #         myfile.close()  
      
        m-=120
        


    

    print "run is finished"
