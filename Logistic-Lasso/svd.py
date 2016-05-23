import numpy as np
import operator
from os import listdir
import matplotlib
import matplotlib.pyplot as plt
from numpy.linalg import *
from scipy.stats.stats import pearsonr
from numpy import linalg as la
import xlrd

def svd(data, S=2):
     
    #calculate SVD
    U, s, Vt = linalg.svd( data )
    V = Vt.T
    #take out columns we don't need
    Sig = np.eye(S)*s[:S]
    newU = U[:,:S]
    newV = V[:,:S]

    #print np.sqrt(Sig[0])
    #print newV

    np.savetxt('dat.txt', np.sqrt(Sig[0,0])*newU[:,0])

     
    # Retrieve dataset 
    # new = U[:,:2]*Sig*V[:2,:]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ['blue','red']
    
   # for i in xrange(num_Sel):
   #     ax.scatter(np.sqrt(Sig[0,0])*newU[i,0],np.sqrt(Sig[1,1])*newU[i,1], color= 'blue')
   # for i in xrange(num_UnS):
    #    ax.scatter(np.sqrt(Sig[0,0])*newU[i+num_Sel,0],np.sqrt(Sig[1,1])*newU[i+num_Sel,1], color= 'red')


    for i in xrange(23):
        ax.scatter(np.sqrt(Sig[0,0])*newV[i,0],np.sqrt(Sig[1,1])*newV[i,1], color= 'gray')
    for i, txt in enumerate(SNPS[0]):
        ax.annotate(txt, (np.sqrt(Sig[0,0])*newV[i,0],np.sqrt(Sig[1,1])*newV[i,1]))


    plt.xlabel('SVD1')
    plt.ylabel('SVD2')
    plt.show()    
    
#Main
  
if __name__=='__main__':
 
# load data points
    gen11 = xlrd.open_workbook("/home/mramezanali/Dropbox/Convex Optimization/BackgroundLiterature/Code/logit/uniqueGenotypes_gen25_allImputations.xls")
    MOR = gen11.sheet_by_name("MOR")
    MUL = gen11.sheet_by_name("MUL")
    LEW = gen11.sheet_by_name("LEW")
    BRI = gen11.sheet_by_name("BRI")
    DOB = gen11.sheet_by_name("DOB")
    STU = gen11.sheet_by_name("STU")
    Selected = gen11.sheet_by_name("Selected")
    Unselected = gen11.sheet_by_name("Unselected")


    SNPS= []
    SNPS.append(MOR.row_values(0,start_colx=0,end_colx=23))

    H = []   
    for r in range(MOR.nrows-1):
        H.append(MOR.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(MUL.nrows-1):
        H.append(MUL.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(LEW.nrows-1):
        H.append(LEW.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(BRI.nrows-1):
        H.append(BRI.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(DOB.nrows-1):
        H.append(DOB.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(STU.nrows-1):
        H.append(STU.row_values(r+1,start_colx=0,end_colx=23))          

    H = np.array(H,dtype=int)
    mean = H.mean(axis=0)
    mean.shape = (1,mean.shape[0])
    H = H-mean #takeout column mean
    mean = H.mean(axis=1)
    mean.shape = (mean.shape[0], 1)
    H = H-mean  #takeout raw mean  
    #fly,snp = np.shape(H)

    num_Sel=MOR.nrows+MUL.nrows+LEW.nrows-3
    num_UnS=BRI.nrows+DOB.nrows+STU.nrows-3

    svd(H,2)




    
