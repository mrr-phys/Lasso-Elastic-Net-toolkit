import numpy as np
from mpmath import mp, exp, erf, erfc
from scipy.integrate import quad
import random

def fq(x0, seta2, seff2):

    coeff = 1./np.sqrt(np.pi)
    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    psi = lamda2*seff2*x0/np.sqrt(2.*seta2)
    ki = 0#theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2)
    
    intgr1 = coeff*(-tau+2*psi+ki+xm)*exp(-(tau+ki+xm)**2) +.5*(1.+2.*(tau-psi)**2)*erfc(tau+ki+xm)
    intgr2 = coeff*(-tau-2*psi+ki-xm)*exp(-(tau+ki-xm)**2) +.5*(1.+2.*(tau+psi)**2)*erfc(tau+ki-xm)
    intgr3 = .5*x0**2*(erf(tau+ki+xm)+erf(tau+ki-xm))

    return (seta2/(1.+lamda2*seff2)**2)*(intgr1+intgr2)+intgr3

def fdq(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    ki = 0#theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2)  

    intgr1 = erfc(tau+ki-xm)
    intgr2 = erfc(tau+ki+xm)

    return .5*(seff2/(1.+lamda2*seff2))*(intgr1+intgr2)


def fl0xz(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    ki = theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2) 

    intgr1 = erf(tau+ki-xm)
    intgr2 = erf(tau+ki+xm)

    return .5*(intgr1+intgr2)

def fl0xnz(x0, seta2, seff2):

    xm = x0/np.sqrt(2.*seta2)
    tau = lamda1*seff2/np.sqrt(2.*seta2)
    ki = theta*(1.+lamda2*seff2)/np.sqrt(2.*seta2) 

    intgr1 = erfc(tau+ki-xm)
    intgr2 = erfc(tau+ki+xm)

    return .5*(intgr1+intgr2)


def fl1(x0, seta2, seff2):

    return None

def fl2(x0, seta2, seff2):

    return None  

    

# sigma_{eta} computed at each run

def sigma_eta2(seta2, seff2):

    #q = (1.-rho)*fq(0,seta2,seff2) + rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fq(x0,seta2,seff2), -np.inf, np.inf)[0] 
    q = (1.-rho)*fq(0,seta2,seff2) + rho*(1./2.)*(fq(+1.0,seta2,seff2)+fq(-1.,seta2,seff2))
    return q/alpha

# sigma_{effective} computed at each run

def sigma_eff2(seta2, seff2):

    #dq =(1.-rho)*fdq(0,seta2,seff2) + rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fdq(x0,seta2,seff2), -np.inf, np.inf)[0] 
    dq =(1.-rho)*fdq(0,seta2,seff2) + rho*(1./2.)*(fdq(+1.,seta2,seff2)+fdq(-1.,seta2,seff2))
    return sigma2 + dq/alpha



def l0xz(seta2, seff2):

    #a =(1.-rho)*fl0xz(0,seta2,seff2)
    #b = rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fl0xz(x0,seta2,seff2), -np.inf, np.inf)[0] 
    a =(1.-rho)*fl0xz(0,seta2,seff2)
    b = rho*(1./2.)*(fl0xz(+1.,seta2,seff2)+fl0xz(-1.,seta2,seff2))

    return a, b

def l0xnz(seta2, seff2):

    #a =(1.-rho)*fl0xnz(0,seta2,seff2)
    #b = rho*(1./np.sqrt(2.*np.pi*sigmax02))*quad(lambda x0:exp(-(x0-mu)**2/(2.*sigmax02))*fl0xnz(x0,seta2,seff2), -np.inf, np.inf)[0] 
    a =(1.-rho)*fl0xnz(0,seta2,seff2)
    b = rho*(1./2.)*(fl0xnz(+1.,seta2,seff2)+fl0xnz(-1.,seta2,seff2))

    return a, b
    

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

    N = 100000 # number of predictors
    K = 100 # number of sparsity
    M = 2000

    xi = 1000.  # Stabilizer
    omega = 1000.

    # initializing parameters

    mu = 0
    sigma2 = 1e-8
    sigmax02 = 1.0
    
    while (M>=K):

        lamda1 = 1.
        lamda2 = 0
        
        alpha = M/float(N)
        rho = K/float(N)        

        for l in range(1):

            theta = 0.

            while (theta <=1.5):

                sigmaeff2 = theta*(1.+lamda2*sigma2)/lamda1+sigma2   # initial \sigma_eff
                sigmaeta2 = theta*(1.+lamda2*sigma2)/lamda1+sigma2  # initial \sigma_eta
                seta2_ini=100
                seff2_ini=100

                print theta               


                itr=0   
                while (np.abs((seta2_ini-sigmaeta2)/sigmaeta2)>1e-5):

                    #print np.abs(seta2_ini-sigmaeta2)

                    itr+=1

                    # update initial \sigma_eff and \sigma_eta]
                    seff2_ini = sigmaeff2
                    seta2_ini = sigmaeta2
                    # compute new \sigma_eta^{2} and \sigma_eff^{2}

                    sigmaeta2 = xi/(1.+xi)*seta2_ini + sigma_eta2(seta2_ini,seff2_ini)/(1.+xi)
                    sigmaeff2 = omega/(1.+omega)*seff2_ini + sigma_eff2(seta2_ini,seff2_ini)/(1.+omega)
               

                    if itr%50==0:
                        print ""
                        print np.abs((seta2_ini-sigmaeta2)/sigmaeta2)
                        print""
                        print M, itr, alpha*sigmaeta2, lamda1*sigmaeff2/np.sqrt(sigmaeta2)    # keeping track of the variables on the screen


                FP, TP = l0xnz(sigmaeta2, sigmaeff2)
                TN, FN = l0xz(sigmaeta2, sigmaeff2)

                TPR = TP/(TP+FN) # True positive rate (sensitivity)
                TNR = TN/(TN+FP) # True negative rate (specificity)


                with open("datROC{}{}.txt".format(M, l), "a") as myfile:

                    # write error data for each different \delta = M/N
                    myfile.write(str(1.-TNR)+' ' + str(TPR)+'\n')
                    myfile.close()
                    
                with open("TPR{}{}.txt".format(M, l), "a") as myfile:

                    # write error data for each different \delta = M/N
                    myfile.write(str(theta)+' ' + str(TPR)+'\n')
                    myfile.close()

                with open("TNR{}{}.txt".format(M, l), "a") as myfile:

                    # write error data for each different \delta = M/N
                    myfile.write(str(theta)+' ' + str(TNR)+'\n')
                    myfile.close()
                    
                theta +=.1                   
                print ""

            #lamda2 += 1.0
            
        M -=20



    print "run is finished"


