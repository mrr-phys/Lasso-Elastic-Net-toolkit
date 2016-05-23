
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:25:13 2015

@author: mohammad
"""


import cvxopt as cvx
from cvxopt import solvers
import numpy as np
import pylab
import time

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

def cvx_l1(H, y, f_a, a):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
    
    c_a = np.zeros(n)
    c_a[a]=1.
   
    f_prtrb=np.zeros(4*n)
    P=np.zeros((2*n,2*n))
    P[a,a]=1
    P[a+n,a+n]=1
    
    if f_a>=0:
        f_prtrb[a]=-f_a
        f_prtrb[a+2*n]=f_a
    else:
        f_prtrb[a+n]=f_a
        f_prtrb[a+3*n]=-f_a        
        
        
    H2 = h([H,-H])
    
    A = u(H2) 
    b = u(y)
   

    #x-split l1 approach  
    c = u(h([np.ones(n)-c_a,np.ones(n)-c_a]).T)
    G1 = u(-np.eye(2*n).T) 
    G2 = u(P.T)
    
    G = u(v([G1,G2]))      
    
    
    hh = u(f_prtrb.T)
    
  #  y = u(h([np.dot(Hd,y),-np.dot(Hd,y)
 


    sol = cvx.solvers.lp(c, G, hh, A,b)

    doublex = np.array(sol['x'][:2*n]).T[0]
    xestimate = np.dot(u(h([np.eye(n),-np.eye(n)])),doublex)
    
    return xestimate.T

def cvx_lp2(H, y):
    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
       
        
        
    H2 = h([H,-H])
    
    A = u(H2) 
    b = u(y)
   

    #x-split l1 approach  
    c = u(h([np.ones(n),np.ones(n)]).T)
    G = u(-np.eye(2*n).T) 

    
    G = u(G)      
    
    
    hh = u(np.zeros(2*n).T)
    
  #  y = u(h([np.dot(Hd,y),-np.dot(Hd,y)
 


    sol = cvx.solvers.lp(c, G, hh, A,b)

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
    
def F(xest,xsparse,a):
    # return mean square error
    x=np.zeros(n)
    x[a]=xest[a]
  
    return np.sum(np.abs(xest-x))#np.sqrt(sq_norm2(xsparse))

#Launching :

if __name__=='__main__':

    n = 300 # number of predictors
    k = 60 # number of sparsity
    

    
    h = np.hstack
    v = np.vstack
    u = cvx.base.matrix
    
    m =210 # initialize number of observant to sparsity
    
    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    
    np.random.shuffle(x)
    xsparse = x.T
    

      

    while (m>=90):
        
        
        for l in range(1):
            
            p=0
            i=0
            
            np.random.seed()
            
            H = 1/np.sqrt(m)*np.random.normal(0,1,[m,n])
            y = np.dot(H,xsparse) #+ np.random.normal(0,sigma2_n,m) 
            #print y.dot(y)
            #print xsparse.dot(xsparse)
            print H.shape

            x_recover = cvx_lp2(H, y)
      
            
  
          #  subsamp= np.random.randint(H.shape[0],size=m)
          #  ysub= y[subsamp]
          #  Hsub= H[subsamp,:]
        
          #  xpzro = cvx_l1(H, y, p, i)
          #  print xpzro
 
            
           # print mse(xpzro-xsparse)
            
          #  np.savetxt('H.txt', Hsub) 
          #  np.savetxt('b.txt', ysub) 
           
            df_list=np.zeros(81)
            
            for a in range(10):
                
                f_a=-4.
               
                
                for j in range(81):
                        
    
                        xest = cvx_l1(H, y, f_a, a)
                        

                        #du_list[j]=xp[i]-xpzro[i]#avg(xp-xpzro)
                        df_list[j]=F(xest,xsparse,a)
                        f_a+=.1

                    
                p=-4.            
                for r in range(81):
            
                    
                    with open("cost{}_{}.txt".format(m,a), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
                        myfile.write(str(p)+' ' +str(df_list[r])+'\n')
                        myfile.close()
                        p+=.1
                        
                with open("cost{}_{}.txt".format(m,a), "a") as myfile:
                    myfile.write(str(x_recover[a])+' ' +str(x_recover[a]))
                    myfile.close()                        


            f_a=0   # initial size of the perturbation f_a set to zero
            a=0    # initialize the starting site a
        
           # x_f_zro2 = cvx_lp2(H, y, f_a, a)
            x_f_zro1 = cvx_lp1(H, y, f_a, a)
            
          #  np.savetxt('H.txt', Hsub) 
          #  np.savetxt('b.txt', ysub) 
           
            deltau_list1=np.zeros(41)  #makes an array of size 10 for 10 different value of f_a
          #  deltau_list2=np.zeros(10) 
        

         #   start = time.time()


  

            for a in range(n):
                f_a = 0
                  
                for j in range(41):
                        
                 #       x_f_nzro2 = cvx_lp2(H, y, f_a, a)
                        x_f_nzro1 = cvx_lp1(H, y, f_a, a)
                        deltau_list1[j]=x_f_nzro1[a]-x_f_zro1[a]#avg(xp-xpzro)
                 #       deltau_list2[j]=x_f_nzro2[a]-x_f_zro2[a]
                        f_a+=.05

                f_a=0
                for r in range(41):

                    with open("suscf{}_{}.txt".format(m,a), "a") as myfile:
                    #with open("susc.txt", "a") as myfile:
            
                        myfile.write(str(f_a)+' ' +str(deltau_list1[r])+'\n')
                        myfile.close()
                        f_a+=.05
                        
                with open("suscf{}_{}.txt".format(m,a), "a") as myfile:   # Added to recognize the contribution from zero or non-zero signal
                        myfile.write(str(xsparse[a])+' ' +str(xsparse[a]))
                        myfile.close()                    
      
        m-=120
        


    

    print "run is finished"
