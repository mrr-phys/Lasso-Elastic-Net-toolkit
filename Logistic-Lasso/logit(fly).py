#PATHWISE COORDINATE OPTIMIZATION For LOGISTIC_ENET

import numpy as np
from scipy.linalg import inv
import random
from pylab import scatter, show, legend, xlabel, ylabel, contour, title
from matplotlib import rc
import matplotlib.pyplot as plt
import itertools
from scipy.misc import comb
import xlrd




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
"""
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
"""

def choose1_iter(length):

    out =[]
    for i in xrange(length):
        perm = np.zeros(length,dtype=int)
        perm[i]=1
        out.append(perm)
    return out

def choose2_iter(length):

    out =[]
    for i in xrange(length):
        for j in xrange(i):
            perm = np.zeros(length,dtype=int)
            perm[i]=1
            perm[j]=1
            out.append(perm)
    return out

def choose3_iter(length):

    out =[]
    for i in xrange(length):
        for j in xrange(i):
            for k in xrange(j):
                perm = np.zeros(length,dtype=int)
                perm[i]=1
                perm[j]=1
                perm[k]=1
                out.append(perm)
    return out


def polynomial_iter(length):

    out =[]
    for i in xrange(length):
        perm = np.zeros(length,dtype=int)
        perm[i]=1
        out.append(perm)
        
    for i in xrange(length):
        perm = np.zeros(length,dtype=int)
        perm[i]=2
        out.append(perm)

    for i in xrange(length):
        for j in xrange(i):
            perm = np.zeros(length,dtype=int)
            perm[i]=1
            perm[j]=1
            out.append(perm)
            
    for i in xrange(length):
        for j in xrange(i):
            perm = np.zeros(length,dtype=int)
            perm[i]=1
            perm[j]=2
            out.append(perm)

    for i in xrange(length):
        for j in xrange(i):
            perm = np.zeros(length,dtype=int)
            perm[i]=2
            perm[j]=2
            out.append(perm)                   
        
    return out

def support(a,N):
    return np.argsort(a)[::-1][:N]


def map_variables(h0, h1,h2,h3, h4,h5,h6,h7,h8,h9, h10,h11,h12,h13, h14,h15,h16,h17,h18,h19, h20,h21,h22):

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
    h9.shape = (h9.size, 1)
    h10.shape = (h10.size, 1)  
    h11.shape = (h11.size, 1)
    h12.shape = (h12.size, 1)
    h13.shape = (h13.size, 1)
    h14.shape = (h14.size, 1)
    h15.shape = (h15.size, 1)
    h16.shape = (h16.size, 1)
    h17.shape = (h17.size, 1)
    h18.shape = (h18.size, 1)
    h19.shape = (h19.size, 1)
    h20.shape = (h20.size, 1)  
    h21.shape = (h21.size, 1)
    h22.shape = (h22.size, 1)



    out = np.zeros(shape=(h1[:,0].size, 1))

    #m, n = out.shape

    degree = n
    config = polynomial_iter(degree)

    it=0
    for it in xrange(len(config)):
        #print config[it][0], config[it][1], config[it][2], config[it][3], config[it][4], config[it][5],config[it][6], config[it][7], config[it][8],\
         #     config[it][9], config[it][10], config[it][11], config[it][12], config[it][13], config[it][14],config[it][15], config[it][16], \
          #    config[it][17], config[it][18], config[it][19], config[it][20], config[it][21], config[it][22], it+1
        
        r = h0**config[it][0] * h1**config[it][1] * h2**config[it][2] * h3**config[it][3] * h4**config[it][4] * h5**config[it][5] * h6**config[it][6] * h7**config[it][7] * h8**config[it][8] * \
            h9**config[it][9] * h10**config[it][10] * h11**config[it][11] * h12**config[it][12] * h13**config[it][13] * h14**config[it][14] * h15**config[it][15] * h16**config[it][16] * \
            h17**config[it][17] * h18**config[it][19] * h19**config[it][19] * h20**config[it][20] * h21**config[it][21] * h22**config[it][22]
        it+=1
        out = np.append(out,r,axis=1)

    out = np.delete(out, (0), axis=1)
    return out, config





def logistic_test(x0_est, x_est, H, y):

    p = np.zeros(shape=(m, 1))

    h = 1./(1.+np.exp(-x0_est-H.dot(x_est)))

    for it in range(0, h.shape[0]):
        if h[it] > 0.5:
            p[it, 0] = 1
        else:
            p[it, 0] = 0

    return (y[np.where(p == y)].size / float(y.size))


#Main
  
if __name__=='__main__':

    #load the dataset
    #H = np.loadtxt('fly_test.txt', delimiter=',')

    gen11 = xlrd.open_workbook("C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\logit\uniqueGenotypes_gen25_allImputations.xls")
    MOR = gen11.sheet_by_name("MOR")
    MUL = gen11.sheet_by_name("MUL")
    LEW = gen11.sheet_by_name("LEW")
    BRI = gen11.sheet_by_name("BRI")
    DOB = gen11.sheet_by_name("DOB")
    STU = gen11.sheet_by_name("STU")
    Selected = gen11.sheet_by_name("Selected")
    Unselected = gen11.sheet_by_name("Unselected")

    H = []   
    for r in range(MOR.nrows-1):
        H.append(MOR.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(MUL.nrows-1):
        H.append(MUL.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(LEW.nrows-1):
        H.append(LEW.row_values(r+1,start_colx=0,end_colx=23))
         

   # num_Sel=MOR.nrows+MUL.nrows+LEW.nrows-3
    #num_UnS=BRI.nrows+DOB.nrows+STU.nrows-3


    H = np.array(H,dtype=float)
    m, n = H.shape
    #a = np.zeros(shape=(m, n))
    #a  = H>1.
    #H[a.nonzero()] = 1.

    print H.mean(axis=0)

    mirror = -np.eye(n)*H.mean(axis=0)

    #mean.shape = (1,mean.shape[0])
   # H = H-mean #takeout column mean
    #mean = H.mean(axis=1)
    #mean.shape = (mean.shape[0], 1)
    #H = H-mean  #takeout raw mean     

    
    y = np.vstack((np.zeros((m, 1),dtype=H.dtype), np.ones((n, 1),dtype=H.dtype)))
    H = np.vstack((H,mirror))

    m, n = H.shape

    #H=H[np.random.randint(H.shape[0],size=600),:]  #sub-sampling

    #y = np.ones((m,1))

    #obs=H

    H, cfig = map_variables( H[:,0], H[:,1], H[:,2], H[:,3], H[:,4], H[:,5],H[:,6], H[:,7], H[:,8],H[:,9], H[:,10], H[:,11], H[:,12], H[:,13], H[:,14], \
                             H[:,15],H[:,16], H[:,17], H[:,18],H[:,19], H[:,20], H[:,21], H[:,22] )


    epsilon = 0.1 #lambda1 = lambda*(1-epsilon), lambda2 = lambda*epsilon
    logit_accuracy = .85  # accuracy of the classification more than 85%  

    lambda_max = max(np.abs(H.transpose().dot(y-.5)))/(1.-epsilon)  #max and min of lambda
    lambda_min = .001*lambda_max


    step = np.logspace(np.log10(lambda_max),np.log10(lambda_min) , num=100)
            
    x0_current = np.log(np.mean(y)/(1.-np.mean(y)))
    x_current = np.zeros(shape=(H.shape[1], 1),dtype=H.dtype)
    x0_est = 0.
    x_est = np.zeros(shape=(H.shape[1], 1))

    #logit_test = logistic_test(x0_current, x_current,H,y)

    j=0
    while(j<100):# and logit_test < logit_accuracy):


        lambda1 =step[j]*(1.-epsilon)
        lambda2 = step[j]*epsilon

        Hdotx_current = H.dot(x_current)
        pi_current = pi(x0_current, Hdotx_current,100)

        x0_est = x0_current + (4./m)*np.sum(y-pi_current)
        x0_current = x0_est

        for a in xrange(H.shape[1]):

            Hdotx_current = H.dot(x_current)
            pi_current = pi(x0_current, Hdotx_current,100)

            z = x0_current+ Hdotx_current + 4.*(y-pi_current)

            partial_res = .25*H[:,a].dot(z-x0_current-Hdotx_current)+.25*H[:,a].dot(H[:,a])*x_current[a]
            #print partial_res

            x_est[a] = eta_soft(partial_res,lambda1)/(.25*H[:,a].dot(H[:,a])+lambda2)

            x_current = x_est
  
        #logit_test = logistic_test(x0_est, x_est,H,y)
        j+=1
        print x_est[:23]
    print x0_est
    print x_est
    #print np.exp(x_est)
"""
        precision=1e-4
        w = support(x_est[:H.shape[1],0],H.shape[1])
        Omega =  x_est[w] > precision
        Omega = w[Omega[:,0]]
        x_list = []
        loc =[]

        if j==100:
            for i in Omega:
                loc=[item for item in xrange(len(cfig[i])) if cfig[i][item] >= 1.]
                a=np.column_stack((x_est[i],loc))
                #print a
                #np.savetxt('xest.txt', a, delimiter=',')
            #x_list=np.append(x_list,a,axis=1)                

"""
    



"""
    precision=.001

    
    w = support(x_est[:n,0],n)
    Omega =  x_est[w] > precision
    Omega = w[Omega[:,0]]
    z = np.zeros((n,1))
    for i in Omega:
        z[i]=x_est[i]
    x_est=z

    
    Omega=Omega[::-1]
    for l in Omega:
        x_est[l]=0.0
        print x_est
        print logistic_test(x0_est, x_est,obs,y)
    np.savetxt('xest.txt', x_est, delimiter=',') 
"""

print "run is finished"





