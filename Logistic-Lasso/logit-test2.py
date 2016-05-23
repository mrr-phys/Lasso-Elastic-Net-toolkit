#PATHWISE COORDINATE OPTIMIZATION For LOGISTIC_ENET
#Example is for M=118,N=2 variables with higher correlated (here the assumpstion is two body correlation, i.e., degree = 2)
import numpy as np
from scipy.linalg import inv
import random
from pylab import scatter, show, legend, xlabel, ylabel, contour, title
from matplotlib import rc
import matplotlib.pyplot as plt


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

def pi(Hdotx_est, th):

    a = np.zeros(shape=(m, 1))
    x = np.zeros(shape=(m, 1))

    x  = Hdotx_est>th
    a[x.nonzero()] = 1.

    x  = Hdotx_est<-th
    a[x.nonzero()] = 0

    x  = np.abs(Hdotx_est)<th
    a[x.nonzero()] = 1./(1.+np.exp(-Hdotx_est[x.nonzero()]))

    return a

def map_variables(h1, h2):

    #Returns a new array of variables h1, h2, h1**2, h2**2, h1*h2, h1*h2**2, etc...

    h1.shape = (h1.size, 1)
    h2.shape = (h2.size, 1)

    degree = 2
    out = np.ones(shape=(h1[:, 0].size, 1))

    m, n = out.shape

    for i in range(1, degree + 1):
        for j in range(i + 1):
            r = (h1 ** (i - j)) * (h2 ** j)
            out = np.append(out, r, axis=1)

    return out


def logistic_test(x_est, H, y):

    m, n = H.shape
    p = np.zeros(shape=(m, 1))

    h = 1./(1.+np.exp(-H.dot(x_est)))

    for it in range(0, h.shape[0]):
        if h[it] > 0.5:
            p[it, 0] = 1
        else:
            p[it, 0] = 0

    return (y[np.where(p == y)].size / float(y.size))


#Main
  
if __name__=='__main__':

    #load the dataset
    data = np.loadtxt('data.txt', delimiter=',')

    H = data[:, 0:2]
    y = data[:, 2]

    pos = np.where(y == 1)
    neg = np.where(y == 0)
    scatter(H[pos, 0], H[pos, 1], marker='o', c='b')
    scatter(H[neg, 0], H[neg, 1], marker='x', c='r')
    xlabel('h1')
    ylabel('h2')
    legend(['y = 1', 'y = 0'])
    #show()

    m, n = H.shape

    y.shape = (m, 1)

    H = map_variables(H[:, 0], H[:, 1])

    epsilon = .1 #lambda1 = lambda*(1-epsilon), lambda2 = lambda*epsilon
    logit_accuracy = .85    

    lambda_max = max(np.abs(H.transpose().dot(y-.5)))/(1.-epsilon)  #max and min of lambda
    lambda_min = .001*lambda_max


    step = np.logspace(np.log10(lambda_max),np.log10(lambda_min) , num=100)
            
    x_est = np.zeros(shape=(H.shape[1], 1))
    j=0
    r = y
    a_list = []

    logit_test = logistic_test(x_est,H,y)
    
    stat = False

    while(j<100 and logit_test < logit_accuracy):

        lambda1 =step[j]*(1-epsilon)
        lambda2 = step[j]*epsilon

        if stat == False:

            for a in xrange(H.shape[1]):
                
                Hdotx_est = H.dot(x_est)
                pi_est = pi(Hdotx_est,100)

                
                z = Hdotx_est + 4.*(y-pi_est)
                x_est[a] = eta_soft(.25*H[:,a].dot(z-Hdotx_est)+.25*H[:,a].dot(H[:,a])*x_est[a],lambda1)/(.25*H[:,a].dot(H[:,a])+lambda2)

          #      if x_est[a] != 0:
           #         if a not in a_list:
            #            a_list.append(a)
             #       elif logit_test > .7:
                        
              #          stat=True
              #          break

        #if stat == True:
           # for a in a_list:
             #   print a_list

              #  Hdotx_est = H.dot(x_est)
             #   pi_est = pi_est(Hdotx_est,100)

                                
               # z = Hdotx_est + 4.*(y-pi_est)                       
              #  x_est[a] = eta_soft(.25*H[:,a].dot(z-Hdotx_est)+.25*H[:,a].dot(H[:,a])*x_est[a],lambda1)/(.25*H[:,a].dot(H[:,a])+lambda2)                                
  
        logit_test = logistic_test(x_est,H,y)
        print logit_test
        j+=1
    print x_est
    #Plot Boundary
    u = np.linspace(-1, 1.5, 50)
    v = np.linspace(-1, 1.5, 50)
    z = np.zeros(shape=(len(u), len(v)))
    for i in range(len(u)):
        for j in range(len(v)):
            z[i, j] = (map_variables(np.array(u[i]), np.array(v[j])).dot(np.array(x_est)))

    z = z.transpose()

    contour(u, v, z)
    plt.rc('text', usetex=True)
    title(r'$\lambda_1 = %f,\, \lambda_2 = %f$' %(lambda1,lambda2))
    xlabel('y1')
    ylabel('y2')
    legend(['y = 1', 'y = 0', 'Decision boundary'])
    show()



print "run is finished"





