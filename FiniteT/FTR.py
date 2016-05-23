#!/usr/bin/python
# algorithm to solve arg min_{x} \frac{1}{2*sigma^{2})*||a.x - b||_{2}^{2} + lambda*||x||_{1} via replica method
# packages to be installed: mpmath, scipy polynomials legendre and laguerre


import matplotlib.pyplot as plt
from numpy.polynomial import legendre as leg  
from numpy.polynomial import laguerre as lag
import numpy as np
from mpmath import mp, exp, erf, erfc
import random

# function to solve <x>_{thermal} problem

def avg_x(eta):
    

    iT = beta/(2.*sigmaeff2)
    
    bm = -beta*lamda*eta
    bp = beta*lamda*eta
    
    tp = np.sqrt(iT)*(-eta+lamda1)
    tm = np.sqrt(iT)*(eta+lamda1)

    nom1 = -(lamda1/2.)*np.sqrt(np.pi/iT)*erfc(tp) + (.5/iT)*exp(-tp**2)
    nom2 = (lamda1/2.)*np.sqrt(np.pi/iT)*erfc(tm) - (.5/iT)*exp(-tm**2)
    dnom1 = .5*np.sqrt(np.pi/iT)*erfc(tp)
    dnom2 = .5*np.sqrt(np.pi/iT)*erfc(tm)

    integrand = (exp(bm)*nom1+exp(bp)*nom2)/(exp(bm)*dnom1+exp(bp)*dnom2)
    
    return integrand

def fq(x0):

    c0=0
    c1=0
    c2=0
    

    for i in range(len(roots_leg)):

        c0 += weights_leg[i]*coeff*(lstd/qstd)*exp(-(lstd*roots_leg[i]-x0)**2/qstd**2)*(avx_leg0[i]+lstd*roots_leg[i]-x0)**2        
    
    for i in range(len(roots_lag)):

        c1 += weights_lag[i]*coeff*exp(roots_lag[i] - (qstd*roots_lag[i]+ lstd+x0)**2/qstd**2)*(avx_lag1[i] -qstd*roots_lag[i]-lstd-x0)**2
        
        c2 += weights_lag[i]*coeff*exp(roots_lag[i] - (qstd*roots_lag[i]+ lstd-x0)**2/qstd**2)*(avx_lag2[i] + qstd*roots_lag[i]+lstd-x0)**2


    return c0+c1+c2

# function to solve <x**2>_{thermal} problem

def avg_x2(eta, avx):
   
    iT = beta/(2.*sigmaeff2)
    
    bm = -beta*lamda*eta
    bp = beta*lamda*eta
    
    tp = np.sqrt(iT)*(-eta+lamda1)
    tm = np.sqrt(iT)*(eta+lamda1)

    nom1 = (tp/(2.*iT*np.sqrt(iT)) - (lamda1+avx)/iT)*exp(-tp**2) + (1./(2.*iT*np.sqrt(iT)) + (lamda1+avx)**2/np.sqrt(iT))*(np.sqrt(np.pi)/2.)*erfc(tp)
    nom2 = (tm/(2.*iT*np.sqrt(iT)) - (lamda1-avx)/iT)*exp(-tm**2) + (1./(2.*iT*np.sqrt(iT)) + (lamda1-avx)**2/np.sqrt(iT))*(np.sqrt(np.pi)/2.)*erfc(tm)
    dnom1 = .5*np.sqrt(np.pi/iT)*erfc(tp)
    dnom2 = .5*np.sqrt(np.pi/iT)*erfc(tm)

    integrand = (exp(bm)*nom1+exp(bp)*nom2)/(exp(bm)*dnom1+exp(bp)*dnom2)   

    return integrand  


def fdq(x0):

    c0=0
    c1=0
    c2=0
    

    for i in range(len(roots_leg)):

        c0 += weights_leg[i]*coeff*(lstd/qstd)*exp(-(lstd*roots_leg[i]-x0)**2/qstd**2)*avx2_leg0[i]
    
    for i in range(len(roots_lag)):

        c1 += weights_lag[i]*coeff*exp(roots_lag[i] - (qstd*roots_lag[i]+ lstd+x0)**2/qstd**2)*avx2_lag1[i]
        
        c2 += weights_lag[i]*coeff*exp(roots_lag[i] - (qstd*roots_lag[i]+ lstd-x0)**2/qstd**2)*avx2_lag2[i]

    return c0+c1+c2

# sigma_{effective} computed at each run

def sigma_eff2(dq):
    return sigma2*(1.+ beta*dq/(sigma2*delta))

# sigma_{eta} computed at each run

def sigma_eta2(q):
    return q/delta


#Main
  
if __name__=='__main__':

    
  #N=int(raw_input("Enter column value: "))
  #K=int(raw_input("Enter sparsity value: "))

  N = 200
  K = 10

  mp.dps = 30   # set high precision to 30 

  avg_x = np.vectorize(avg_x)
  avg_x2 = np.vectorize(avg_x2)


  Leg = leg.leggauss(100)
  roots_leg = Leg[0]
  weights_leg = Leg[1]

  Lag = lag.laggauss(100)
  roots_lag = Lag[0]
  weights_lag = Lag[1]  


  mu = 0
  sigma2 = .001           # \sigma parameter I vary sigma2 rather than lamda as I found making smaller lambda doesnt lead to a convergent result
  sigmax02 =1.
  lamda = 1.               #\lambda parameter

  M = K
  while(M<=N):


      ##############
      #producing sample sparse array
      s = np.zeros(K) 

      for i in range(K):
          s[i] = 2*random.random()-1

      ###################
     
      #s = np.array([1.03602, -0.288999, -1.36459, -0.213365, -0.0476977, -2.5811, -0.0263163, -0.749791, 2.02066, -0.0240022])  # imported sparse coefficients
      
      delta = M/float(N)
      rho = K/float(N)

      sigmaeff2 = 1.1   # initial \sigma_eff
      sigmaeta2 = lamda*sigmaeff2   # initial \sigma_eta


      # compute the fluctuation between two adjacent steps for \sigma_eff and \sigma_eta  

      fluc_seff2 = [1.]  
      fluc_seta2 = [1.]

      m=0
      beta = .1                    # initial inverse temperature

      #q = delta*sigmaeta2
      #dq = (sigmaeff2 - 1.)*(2.*delta/beta)

      while(fluc_seff2[m]>1e-3):

          if beta <1.:
              beta +=.1
          else:
              beta+=.1


          m+=1
          # update initial \sigma_eff and \sigma_eta
          seff2_ini = sigmaeff2   
          seta2_ini = sigmaeta2

          #############################################################
          # the part of the code that compute the thermal average <x>^2 and <(x-<x>)^2>


          a = np.sqrt(4.)  # just a coefficient
          coeff = 1./np.sqrt(np.pi)
          qstd = np.sqrt(2.*sigmaeta2)
          lamda1 = lamda*sigmaeff2
          lstd =  lamda1 + a*np.sqrt(sigmax02)



          avx_leg0 = avg_x(lstd*roots_leg)
          avx_lag1 = avg_x(-qstd*roots_lag-lstd)
          avx_lag2 = avg_x(qstd*roots_lag+lstd)

          avx2_leg0 = avg_x2(lstd*roots_leg, avx_leg0)
          avx2_lag1 = avg_x2(-qstd*roots_lag-lstd, avx_lag1)
          avx2_lag2 = avg_x2(qstd*roots_lag+lstd, avx_lag2)


          ##############################################################
          
            
          q0 = [fq(0)]*(N-K)        #q is computed when sparse coefficients are zero         
          dq0 = [fdq(0)]*(N-K)
          

          for i in range(len(s)):
              x0 = s[i]
              q0.append(fq(x0))
              dq0.append(fdq(x0))
              
          q0 = np.array(q0, dtype=np.float)
          dq0 = np.array(dq0, dtype=np.float)

          #print q0


          q = np.sum(q0)/N     # compute <x>^2 (MSE) averaged over thermal and quenched disorder, this is the q term in the notes
          dq = np.sum(dq0)/N   # compute <(x - <x>)^2> averaged over thermal and quenched disorder, this is the Q-q term in the notes

          # compute new \sigma_eta^{2} and \sigma_eff^{2}
          sigmaeta2 = sigma_eta2(q)  
          sigmaeff2 = sigma_eff2(dq)

          fluc_seff2.append(np.abs(sigmaeff2-seff2_ini)/sigmaeff2)
          fluc_seta2.append(np.abs(sigmaeta2-seta2_ini)/sigmaeta2)


          print 1./beta, q, dq, sigmaeff2/beta, sigmaeff2, sigmaeta2    # keeping track of the variables on the screen
          #print ""
          #print np.abs(sigmaeff2-seff2_ini)/sigmaeff2
          #print np.abs(sigmaeta2-seta2_ini)/sigmaeta2
          


      with open("dat.txt", "a") as myfile:
          # write error data for each different \delta = M/N
          myfile.write(str(M/float(N))+' '+str(q)+'\n')
          myfile.close()
          
      M+=5    # increment of undersampling step
      
print "run is finished"
