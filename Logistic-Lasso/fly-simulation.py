#!/usr/bin/env python

import numpy as np
from scipy import sparse
import itertools


def Energy(obs,J,K):
    "Energy of each observation at particular configuration"
    Ene = obs.dot(J) +.5*np.diag(obs.dot(K.dot(obs.T)))
    return Ene

def Sample(Ene):
    #
    p = 1./(1.+np.exp(-Ene))
    th = np.random.rand()#np.random.normal(p.mean(axis=0),.01,m)

    x  = np.where(p>th)
    return x

def Combs(seq1,seq2):

    out = []
    num = np.random.randint(5,6)

    
    seq1=[seq1[i::num] for i in range(num)]
    seq2=[seq2[i::num] for i in range(num)]
    for i in range(num):
        a = np.random.randint(2)
          
        if a==0:
            out = np.append(out,seq1[i])
        else:
            out = np.append(out,seq2[i])

    return out          


if __name__ == '__main__':

    gen = 150      # total number of generations 
    m = 2000      # number of flies in each generation
    n = 10
    

    J = np.zeros(n)#.001*np.random.rand(n)

    K = np.zeros((n,n))
    K[0,2]=1.
    K[2,0]=1.


    obs = np.zeros((m,n))

   
    #obs = np.array(np.round(sparse.rand(m,n,density=.5).todense()).astype('int'))
    i = 0
    while (i<m):
        col = np.random.randint(n)
        obs[i,col]=1
        i+=1

    it = 0                  

    while (it<gen):

        Ene = Energy(obs,J,K)
        Smpl = Sample(Ene)[0]

        j = 0
        seq0 = obs[0,:]

        T=seq0
        
        while (j<m):
            
            a = np.random.randint(len(Smpl))
            b = np.random.randint(len(Smpl))
            if a != b:
                seq1=obs[a,:]
                seq2=obs[b,:]
                new_seq = Combs(seq1,seq2)
                
                T = np.vstack((T,new_seq))            
            j+=1
        
        T = np.delete(T, (0), axis=0)
        obs = T.astype('int')

        m,n = T.shape

        it+=1
    print obs.mean(axis=0)
 

    np.savetxt('fly_test0.txt', obs, delimiter=',') 


        








