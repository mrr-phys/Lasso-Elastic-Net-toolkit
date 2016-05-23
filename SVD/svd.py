import numpy as np
import operator
from os import listdir
import matplotlib
import matplotlib.pyplot as plt
from numpy.linalg import *
from scipy.stats.stats import pearsonr
from numpy import linalg as la
import xlrd

def svd(data, S):
     
    #calculate SVD
    U, s, Vt = linalg.svd( data )
    V = Vt.T
    #take out columns we don't need
    Sig = np.eye(S)*s[:S]
    newU = U[:,:S]
    newV = V[:,:S]

    #print np.sqrt(Sig[0])
    #print newV

   # np.savetxt('dat.txt', np.sqrt(Sig[0,0])*newU[:,0])

     
    # Retrieve dataset 
    print Sig[9,9]
   # print newU
    print (np.sqrt(n)-np.sqrt(m))*(1./np.sqrt(m))
    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1)
    colors = ['blue','red']
    
   # for i in xrange(num_Sel):
   #     ax.scatter(np.sqrt(Sig[0,0])*newU[i,0],np.sqrt(Sig[1,1])*newU[i,1], color= 'blue')
   # for i in xrange(num_UnS):
    #    ax.scatter(np.sqrt(Sig[0,0])*newU[i+num_Sel,0],np.sqrt(Sig[1,1])*newU[i+num_Sel,1], color= 'red')


  #  for i in xrange(n):
    #    ax.scatter(Sig[0,0]*newV[i,0],Sig[1,1]*newV[i,1], color= 'gray')
   

  #  plt.xlabel('SVD1')
  #  plt.ylabel('SVD2')
  #  plt.show()    
    
#Main
  
if __name__=='__main__':
 
# load data points
    n=1000
    
    m=10
    
    
    H=(1./np.sqrt(m))*np.random.randn(m,n)
    
    
    
    
    
    
    

        

  #  H = np.array(H,dtype=int)
  #  mean = H.mean(axis=0)
  #  mean.shape = (1,mean.shape[0])
  #  H = H-mean #takeout column mean
  #  mean = H.mean(axis=1)
   # mean.shape = (mean.shape[0], 1)
   # H = H-mean  #takeout raw mean  
    #fly,snp = np.shape(H)

    svd(H,m)




    
