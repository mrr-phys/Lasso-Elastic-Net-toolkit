"""
Created on Apr 20 2014

"""

# Zero Temperature Code for Elastic_Net \lambda_{1}|x|+\frac{1}{2}\lambda_{2}|x|^2
# Computes and stores \sigmaeta2 and sigmaeff2 for given \rho and \alpha
# Computes and stores MSE(error), L0, L1 and L2 norm for given \rho and \alpha


import numpy as np
from mpmath import mp, exp, erf, erfc
from scipy.integrate import quad
import random
import math

def fq(x0, seta2, seff2):

    coeff = 1./np.sqrt(np.pi)
    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    psi = lamda2*seff2*x0/np.sqrt(2.*seta2)
    
    intgr1 = coeff*(-tau+2*psi+xm)*exp(-(tau+xm)**2) +.5*(1.+2.*(tau-psi)**2)*erfc(tau+xm)
    intgr2 = coeff*(-tau-2*psi-xm)*exp(-(tau-xm)**2) +.5*(1.+2.*(tau+psi)**2)*erfc(tau-xm)
    intgr3 = .5*x0**2*(erf(tau+xm)+erf(tau-xm))

    return (seta2/(1.+lamda2*seff2)**2)*(intgr1+intgr2)+intgr3

def fdq(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)  

    intgr1 = erfc(tau-xm)
    intgr2 = erfc(tau+xm)

    return .5*(seff2/(1.+lamda2*seff2))*(intgr1+intgr2)


def pxz(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    ki = theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2) 

    intgr1 = erf(tau+ki-xm)
    intgr2 = erf(tau+ki+xm)

    return .5*(intgr1+intgr2)

def pxnz(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    ki = theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2) 

    intgr1 = erfc(tau+ki-xm)
    intgr2 = erfc(tau+ki+xm)

    return .5*(intgr1+intgr2)


def fl1(x0, seta2, seff2):

    coeff = 1./np.sqrt(np.pi)
    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)

    intgr1 = coeff*exp(-(tau-xm)**2) -(tau-xm)*erfc(tau-xm)
    intgr2 = coeff*exp(-(tau+xm)**2) -(tau+xm)*erfc(tau+xm)

    return .5*(np.sqrt(2.*seta2)/(1.+lamda2*seff2))*(intgr1+intgr2)

def fl2(x0, seta2, seff2):

    coeff = 1./np.sqrt(np.pi)
    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)

    intgr1 = -coeff*(tau-xm)*exp(-(tau-xm)**2) +(.5+(tau-xm)**2)*erfc(tau-xm)
    intgr2 = -coeff*(tau+xm)*exp(-(tau+xm)**2) +(.5+(tau+xm)**2)*erfc(tau+xm)

    return (seta2/(1.+lamda2*seff2)**2)*(intgr1+intgr2)  

    

# sigma_{eta} computed at each run
# Computes for Gaussian dist. with zero mean and variance \mu
# Computes for Bernoulli distribution with input \pm 1

def sigma_eta2(seta2, seff2):
 
    q = (1.-rho)*fq(0,seta2,seff2) + rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fq(x0,seta2,seff2), -np.inf, np.inf)[0] 
    #q = (1.-rho)*fq(0,seta2,seff2) + rho*(1./2.)*(fq(+1.0,seta2,seff2)+fq(-1.,seta2,seff2))
    return q/alpha

# sigma_{effective} computed at each run

def sigma_eff2(seta2, seff2):

    dq =(1.-rho)*fdq(0,seta2,seff2) + rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fdq(x0,seta2,seff2), -np.inf, np.inf)[0] 
    #dq =(1.-rho)*fdq(0,seta2,seff2) + rho*(1./2.)*(fdq(+1.,seta2,seff2)+fdq(-1.,seta2,seff2))
    return sigma2 + dq/alpha

# Computes L0, L1, L2 norms

def l0(seta2, seff2):

    #a =(1.-rho)*pxnz(0,seta2,seff2)
    #b = rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*pxnz(x0,seta2,seff2), -np.inf, np.inf)[0] 
    a =(1.-rho)*pxnz(0,seta2,seff2)
    b = rho*(1./2.)*(pxnz(+1.,seta2,seff2)+pxnz(-1.,seta2,seff2))

    return a+b
    

def l1(seta2, seff2):

    #a = (1.-rho)*fl1(0,seta2,seff2)
    #b = rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fl1(x0,seta2,seff2), -np.inf, np.inf)[0]    
    a = (1.-rho)*fl1(0,seta2,seff2)
    b = rho*(1./2)*(fl1(+1,seta2,seff2)+fl1(-1,seta2,seff2))

    return lamda1*(a+b)

def l2(seta2, seff2):

    #a =(1.-rho)*fl2(0,seta2,seff2)
    #b = rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fl2(x0,seta2,seff2), -np.inf, np.inf)[0]    
    a =(1.-rho)*fl2(0,seta2,seff2)
    b = rho*(1./2)*(fl2(+1,seta2,seff2)+fl2(-1,seta2,seff2))

    return lamda2*(a+b)
    
    

   
#Main
  
if __name__=='__main__':

    N = 200 # number of predictors
    K = 20  # number of sparsity
    M = 30 # initial number of measurements

    xi = 1.  # Converging Stabilizer Parameters
    omega = 1.

    # initializing parameters

    mu = 0
    sigma2 = 1.#1e-6#1e-10
    sigmax02 = 1.0

    lamda1 = 1e-4
    lamda2 = 0#8e-5

    theta = 0 # the hard cutoff parameter over the P(x|x_0)

    s_cutoff = 1e-9 # Sigma Cutoff

    while (M>=30):

            alpha = M/float(N)
            rho = K/float(N)

            itr=0

            #sigmaeff2 = sigma2*s_cutoff   # initialize \sigmaeff2
            #sigmaeta2 = sigma2*s_cutoff   # initialize \sigmaeta2
            #seff2_ini = 0
            #seta2_ini = 0   
            sigmaeff2 = sigma2  # initialize \sigmaeff2
            sigmaeta2 = 1e-7   # initialize \sigmaeta2        

            while (itr<=120):#while (np.abs((seta2_ini-sigmaeta2)/sigmaeta2)>s_cutoff):
                #print np.abs(seta2_ini-sigmaeta2)

                itr+=1

                # update initial \sigma_eff and \sigma_eta
                seff2_ini = sigmaeff2
                seta2_ini = sigmaeta2
                
                # compute new \sigma_eta^{2} and \sigma_eff^{2}
                sigmaeta2 = xi/(1.+xi)*seta2_ini + sigma_eta2(seta2_ini,seff2_ini)/(1.+xi)
                sigmaeff2 = omega/(1.+omega)*seff2_ini + sigma_eff2(seta2_ini,seff2_ini)/(1.+omega)

                if itr%20==0:
                    print ""
                    print np.abs((seta2_ini-sigmaeta2)/sigmaeta2)
                    print""
                    print itr, rho, alpha, alpha*sigmaeta2, lamda1*sigmaeff2/np.sqrt(sigmaeta2)    # keeping track of the variables on the screen

            #L0 = l0(sigmaeta2, sigmaeff2) # Computes L0 norm(non-zero elements over number of predictors)
           # L1 = l1(sigmaeta2, sigmaeff2) # Computes L1 norm(magnitude of the potential term)


            with open("MSE.8.txt", "a") as myfile:
                # write error data for each different \delta = M/N
                myfile.write(str(alpha)+' ' +str(sigmaeta2*alpha)+'\n')
                myfile.close()
                
      #      with open("L0.txt", "a") as myfile:
                # write L0 data for each different \delta = M/N
       #         myfile.write(str(alpha)+' ' +str(L0)+'\n')
       #         myfile.close()
                
        #    with open("L1.txt", "a") as myfile:
                # write L1 data for each different \delta = M/N
        #        myfile.write(str(alpha)+' ' +str(L1)+'\n')
        #        myfile.close()                 
                
            M-=5 # decrement of undersampling measurements
       


    print "run is finished"
            
