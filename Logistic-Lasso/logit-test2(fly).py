#PATHWISE COORDINATE OPTIMIZATION For LOGISTIC_ENET

import numpy as np
from scipy.linalg import inv
import random
from pylab import scatter, show, legend, xlabel, ylabel, contour, title
from matplotlib import rc
import matplotlib.pyplot as plt
import itertools
from scipy.misc import comb


def eta_soft(t, theta):

    if np.abs(t) <= theta:
        return 0

    if t > theta:
        return t-theta

    if t < -theta:
        return t+theta
    
def norm2(a):
    #return ||a||_2
    return np.sqrt(a.transpose().dot(a))

def mse(a):
    # return mean square error
    return (1./np.size(a))*norm2(a)**2

def pi(x0, Hdotx_est, th):

    a = np.zeros(shape=(m, 1))
    x = np.zeros(shape=(m, 1))

    x  = x0+Hdotx_est>th
    a[x.nonzero()] = 1.

    x  = x0+Hdotx_est<-th
    a[x.nonzero()] = 0

    x  = np.abs(x0+Hdotx_est)<th
    a[x.nonzero()] = 1./(1.+np.exp(-(x0+Hdotx_est)[x.nonzero()]))

    return a

def map_variables(h0, h1,h2,h3, h4,h5,h6,h7,h8):

    #Kernel chosen to be the map of all Permutation (returns a new array of variables h1, h2, h1*h2, h1*h2*h3, etc...)
    
    h0.shape = (h0.size, 1)  
    h1.shape = (h1.size, 1)
    h2.shape = (h2.size, 1)
    h3.shape = (h3.size, 1)
    h4.shape = (h4.size, 1)
    h5.shape = (h5.size, 1)
    h6.shape = (h6.size, 1)
    h7.shape = (h7.size, 1)
    h8.shape = (h8.size, 1)
    #h9.shape = (h9.size, 1)    
  

    out = np.zeros(shape=(h1[:,0].size, 1))

    #m, n = out.shape

    degree = n
    perm = np.zeros(n,dtype=int)

    it=0

    for j in range(degree):
        perm[j]+=1
        k = int(round(comb(degree,j+1)))
        for i in range(k):
            config = list(set(itertools.permutations(perm,degree)))[i]
            print config[0], config[1], config[2], config[3], config[4], config[5],config[6], config[7], config[8], it+1
            r = h0**config[0] * h1**config[1] * h2**config[2] * h3**config[3] * h4**config[4] * h5**config[5] * h6**config[6] * h7**config[7] * h8**config[8] 
            it+=1
            out = np.append(out, r, axis=1)
            if j==degree-1:
                break
    out = np.delete(out, (0), axis=1)
    return out


def logistic_test(x0, x_est, H, y):

    p = np.zeros(shape=(m, 1))

    h = 1./(1.+np.exp(-x0-H.dot(x_est)))

    for it in range(0, h.shape[0]):
        if h[it] > 0.5:
            p[it, 0] = 1
        else:
            p[it, 0] = 0

    return (y[np.where(p == y)].size / float(y.size))


#Main
  
if __name__=='__main__':

    #load the dataset
    H = np.loadtxt('fly_test.txt', delimiter=',')
    m, n = H.shape

    print m, n

    y = np.ones((m,1))

    H = map_variables( H[:,0], H[:,1], H[:,2], H[:,3], H[:,4], H[:,5],H[:,6], H[:,7], H[:,8] )

    epsilon = 0.2 #lambda1 = lambda*(1-epsilon), lambda2 = lambda*epsilon
    logit_accuracy = .85  # accuracy of the classification more than 85%  

    lambda_max = max(np.abs(H.transpose().dot(y-.5)))/(1.-epsilon)  #max and min of lambda
    lambda_min = .001*lambda_max


    step = np.logspace(np.log10(lambda_max),np.log10(lambda_min) , num=100)
            
    x_est = np.zeros(shape=(H.shape[1], 1))
    x0 = 0

    logit_test = logistic_test(x0, x_est,H,y)
    
    stat = False
    j=0

    while(j<100):# and logit_test < logit_accuracy):

        lambda1 =step[j]*(1.-epsilon)
        lambda2 = step[j]*epsilon

        if stat == False:

            Hdotx_est = H.dot(x_est)
            pi_est = pi(x0, Hdotx_est,100)

            x0 = x0 + (4./m)*np.sum(y-pi_est)

            for a in xrange(H.shape[1]):
                
                Hdotx_est = H.dot(x_est)
                pi_est = pi(x0, Hdotx_est,100)

                
                z = x0+ Hdotx_est + 4.*(y-pi_est)
                x_est[a] = eta_soft(.25*H[:,a].dot(z-x0-Hdotx_est)+.25*H[:,a].dot(H[:,a])*x_est[a],lambda1)/(.25*H[:,a].dot(H[:,a])+lambda2)
                             
  
        logit_test = logistic_test(x0, x_est,H,y)
        j+=1
    print logit_test
    print x0
    print x_est
    np.savetxt('xest.txt', x_est, delimiter=',') 


print "run is finished"





