#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Oct. 7, 2014. Code based on Deanna Needell and Joel Tropp CoSaMP paper

import numpy as np
from scipy.linalg import inv
import random


def norm2(a):
    #return ||a||_2
    return np.sqrt(a.transpose().dot(a))

def supp(a,N):
    return np.argsort(a)[::-1][:N]


#Main

if __name__=='__main__':

    n = 200 # number of predictors
    k = 10 # sparsity level
    m = 100

    ##############
    #producing sample sparse array

    x = np.zeros(n)
    x[:k] = np.random.normal(0,1,k)
    #np.random.shuffle(x)
    xsparse = x.T
    
 #   x0 = np.zeros(n)
 #   loc = random.sample(range(n),k)
   # for i in range(k):
   #     x0[loc[i]] = np.random.randn()
    ##################

    print xsparse[:k]

    H = np.random.normal(0,1./np.sqrt(m),(m,n))   # sensing  matrix
    y = H.dot(xsparse)  # samples

    

    
    r = y   #Current residual
    j = 1   #jth iteration
    precision = 1e-8
    max_itr = 350
    dist = 1e-5

    T = []

    while(j <= max_itr and norm2(r)/float(norm2(y)) > dist):
        
        err = np.abs(H.transpose().dot(r))  #Compute current error

        w = supp(err,2*k)        #Support identiﬁcation (the best 2k support set)

        Omega = np.where( err[w] > precision)
        Omega = w[Omega]


        T = np.union1d(Omega,T).astype('int')  #Strongest support merging


        HT = H[:,T]
        Ht = inv(HT.transpose().dot(HT)).dot(HT.transpose())
        bT = Ht.dot(y)   #Least-square signal estimation
        

        #Pruning

        w = supp(np.abs(bT),k)
        
        
        updated_indices = np.where( np.abs(bT[w]) > precision)
        updated_indices = w[updated_indices]
        
        T = T[updated_indices]
        x = np.zeros(H.shape[1])
        x[T] = bT[w]

        
        r = y - H.dot(x)

        
        j +=1

    print x
    print np.sum((x-xsparse)**2)

