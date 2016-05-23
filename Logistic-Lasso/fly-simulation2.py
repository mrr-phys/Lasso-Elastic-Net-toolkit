#!/usr/bin/env python

import numpy as np
from scipy import sparse
import itertools


def Energy(seq,J,K):
    "Energy of each observation at particular configuration"
    Ene = seq.dot(J) +.5*seq.dot(K.dot(seq.T))
    return Ene

def Combs(seq1,seq2):

    out = []
    num = 4#np.random.randint(2,n/2)
    
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

    gen = 15      # total number of generations 
    m = 1000      # number of flies in each generation
    n = 9
    

    J = .01*np.random.rand(n)

    K = -.01*np.ones((n,n))
    K[6,0]=1.
    K[0,6]=1.
    K[6,6]=1.
    K[0,0]=1.
   
  

    obs = np.zeros((m,n))

   
    #obs = np.array(np.round(sparse.rand(m,n,density=.5).todense()).astype('int'))
    i = 0
    while (i<m):
        col = np.random.randint(n)
        obs[i,col]=1
        i+=1

    it = 0

    while (it<gen):

        j = 0
        seq0 = obs[0,:]

        T=seq0
        
        while (j<m):
            
            a = np.random.randint(1,m)
            b = np.random.randint(1,m)
            if a != b:
                seq1=obs[a,:]
                seq2=obs[b,:]

                Ene_ini = max(Energy(seq1,J,K), Energy(seq2,J,K))
               # print seq1
               # print seq2
               # print Ene_ini
                new_seq = Combs(seq1,seq2)
                Ene_fin = Energy(new_seq,J,K)
                if Ene_fin>Ene_ini:
                    T = np.vstack((T,new_seq))
                    j+=1
                else:
                   # if np.random.randint(2)==0:
                       # T = np.vstack((T,new_seq))
                       # j+=1
                    
                    p = 1./(1.+np.exp(-(Ene_fin-Ene_ini)))
                    th = np.random.rand()
                    if p>th:
                        T = np.vstack((T,new_seq))
                        j+=1   

        T = np.delete(T, (0), axis=0)
        obs = T.astype('int')

        m,n = T.shape

        it+=1
        
        if it%5==0:
            print obs.mean(axis=0)
 

    np.savetxt('fly_test.txt', obs, delimiter=',') 


        








