#PATHWISE COORDINATE OPTIMIZATION
import numpy as np
from scipy.linalg import inv
import random



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


    
#Main
  
if __name__=='__main__':

    n = 200 # number of predictors
    k = 30 # number of sparsity
    m = 100
    
    s_itr = 100      
    lambda2 = 0
    sigma2 = 1.          # sigma2
    dist = 1e-5

    while (m <=100):

        H = np.random.normal(0,1./np.sqrt(m),(m,n))   # sensing  matrix

        s_list = []

        for i in range(s_itr):

            #producing sample sparse array

            s = np.zeros(n)

            loc = random.sample(range(n),k)

            for i in range(k):
                s[loc[i]] = np.random.randn()
                    
         ##################

            y = H.dot(s)
            lambda_max = max(np.abs(H.transpose().dot(y)))/sigma2
            lambda_min = .0001*lambda_max

            step = np.logspace(np.log10(lambda_max),np.log10(lambda_min) , num=100)

            
            x_estimate = np.zeros(n)
            j=0
            r = y
            F_test = norm2(r)/float(norm2(y))

            a_list = []
            stat = False

            while(j<100 and F_test > dist):

                F_test = norm2(r)/float(norm2(y))

                lambda1 =step[j]

                if stat == False:

                    for a in xrange(n):
                        x_estimate[a] = eta_soft((1./sigma2)*H[:,a].dot(y-H.dot(x_estimate)+H[:,a]*x_estimate[a]),lambda1)/((1./sigma2)*H[:,a].dot(H[:,a])+lambda2)


                        if x_estimate[a] != 0:
                            if a not in a_list:
                                a_list.append(a)

                            elif F_test <1e-2:
                                stat=True
                                break

                if stat == True:
                    for a in a_list:
                        x_estimate[a] = eta_soft((1./sigma2)*H[:,a].dot(y-H.dot(x_estimate)+H[:,a]*x_estimate[a]),lambda1)/((1./sigma2)*H[:,a].dot(H[:,a])+lambda2)                                  
  
                r = y - H.dot(x_estimate)
                j+=1

                
            s_list.append(mse(x_estimate-s))


        with open("MSE_exp.txt", "a") as myfile:
            # write error data for each different \delta = M/N
            myfile.write(str(m/float(n))+' ' +str(np.sum(s_list)/s_itr)+'\n')
            myfile.close()
        m+=10 # increment of undersampling step


    print "run is finished"





