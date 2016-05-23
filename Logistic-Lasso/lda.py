import numpy as np
import operator
from os import listdir
import matplotlib
import matplotlib.pyplot as plt
from numpy.linalg import *
from scipy.stats.stats import pearsonr
from numpy import linalg as la
from scipy.linalg import inv
import xlrd
from sklearn.lda import LDA
from sklearn import preprocessing
from sklearn.decomposition import PCA as sklearnPCA

def lda(data, S):
     
    #calculate SVD
    U, s, Vt = linalg.svd( data )
    V = Vt.T
  
    egn_val=0
    for i in range(S):
        egn_val+=s[i]**2
    r = (egn_val/np.sum(s**2))*100
    #take out columns we don't need
    Sig = np.eye(S)*s[:S]
    U_sv = U[:,:S]
    V_sv = V[:,:S]

    #print np.sqrt(Sig[0])
    #print newV

    

     
    # Retrieve dataset 
    # new = U[:,:2]*Sig*V[:2,:]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ['blue','red']
    
   # for i in xrange(num_Sel):
     #   ax.scatter(U_sv[i,0],U_sv[i,1], color= 'blue')
    #for i in xrange(num_UnS):
      #  ax.scatter(U_sv[i+num_Sel,0],U_sv[i+num_Sel,1], color= 'red')

    #y = [Y[i,0] for i in range(Y.shape[0])]
      
    # LDA
    U_Sel = U_sv[:num_Sel,:]
    U_UnS = U_sv[num_Sel:num_Sel+num_UnS,:]
    
    inter_cov = np.cov(U_Sel.T)+np.cov(U_UnS.T)
    mean_diff = U_Sel.mean(axis=0)-U_UnS.mean(axis=0)
    lda_coef = inv(inter_cov).dot(mean_diff)

   # print lda_coef

   # sk_lda = LDA()
    #x_lda = sk_lda.fit_transform(U_sv, y)

    #print sk_lda.coef_

    
    # PCA
    
    #sk_pca = sklearnPCA(n_components=2)
    #x_ldapca = sk_pca.fit_transform(x_lda)
    
    V_lda = V_sv[1:,:S].dot(lda_coef)

    V_lda_sort = np.sort(V_lda)[::-1]
    V_lda_argsort = np.argsort(V_lda)[::-1]
    SNPS_sort = [SNPS[0][i] for i in V_lda_argsort]


    #print clf.predict(newU)
    
    #for i in xrange(num_Sel):
     #   ax.scatter(x_lda[i,0],x_lda[i,1], color= 'blue')
    #for i in xrange(num_UnS):
     #   ax.scatter(x_lda[i+num_Sel,0],x_lda[i+num_Sel,1], color= 'red')


   # for i in xrange(23):
    #    ax.scatter(V_lda[i],0, color= 'gray')
   # for i, txt in enumerate(SNPS[0]):
    #    ax.annotate(txt, (V_lda[i],0))


    plt.plot(range(H.shape[1]-1),V_lda_sort/100., 'o', color= 'gray')

    # Rotate x-labels on the x-axis
    fig.autofmt_xdate()

    # Label-values for x and y axis
    plt.xticks(range(H.shape[1]-1),SNPS_sort)

    plt.title('the %i leading SVD components captured %.2f percent' %( S, r))

    plt.ylabel('Contribution')
    plt.xlim(-.3, H.shape[1]-1)

    #plt.xlabel('SVD1')
   # plt.ylabel('SVD2')
    plt.show()    
    
#Main
  
if __name__=='__main__':
 
# load data points
    gen11 = xlrd.open_workbook("C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\logit\uniqueGenotypes_gen25_allImputations.xls")
    MOR = gen11.sheet_by_name("MOR")
    MUL = gen11.sheet_by_name("MUL")
    LEW = gen11.sheet_by_name("LEW")
    BRI = gen11.sheet_by_name("BRI")
    DOB = gen11.sheet_by_name("DOB")
    STU = gen11.sheet_by_name("STU")
    Sel = gen11.sheet_by_name("Selected")
    Uns = gen11.sheet_by_name("Unselected")


    SNPS= []
    SNPS.append(MOR.row_values(0,start_colx=0,end_colx=23))

    H = []   
  #  for r in range(MOR.nrows-1):
    #    H.append(MOR.row_values(r+1,start_colx=0,end_colx=23))
   # for r in range(MUL.nrows-1):
    #    H.append(MUL.row_values(r+1,start_colx=0,end_colx=23))
   # for r in range(LEW.nrows-1):
    #    H.append(LEW.row_values(r+1,start_colx=0,end_colx=23))
   # for r in range(BRI.nrows-1):
    #    H.append(BRI.row_values(r+1,start_colx=0,end_colx=23))
   # for r in range(DOB.nrows-1):
    #    H.append(DOB.row_values(r+1,start_colx=0,end_colx=23))
    #for r in range(STU.nrows-1):
    #    H.append(STU.row_values(r+1,start_colx=0,end_colx=23))          

   # num_Sel=MOR.nrows+MUL.nrows+LEW.nrows-3
   # num_UnS=BRI.nrows+DOB.nrows+STU.nrows-3

    for r in range(Sel.nrows-1):
        H.append(Sel.row_values(r+1,start_colx=0,end_colx=23))
    for r in range(Uns.nrows-1):
        H.append(Uns.row_values(r+1,start_colx=0,end_colx=23))

    num_Sel=Sel.nrows-1
    num_UnS=Uns.nrows-1       

        

    H = np.array(H,dtype=float)

    
    Y = np.vstack((np.ones((num_Sel, 1),dtype=H.dtype),-np.ones((num_UnS, 1),dtype=H.dtype)))
    H = np.hstack((Y,H))

    #preprocessing.scale(H, axis=0, with_mean=True, with_std=True, copy=False) #scaling to mean zero variance one

    
   # H_Sel = H[:num_Sel,1:]
   # H_UnS = H[num_Sel:num_Sel+num_UnS,1:]
    
   # inter_cov = np.cov(H_Sel.T)+np.cov(H_UnS.T)
    #mean_diff = H_Sel.mean(axis=0)-H_UnS.mean(axis=0)
    #print inv(inter_cov).dot(mean_diff)

    mean = H.mean(axis=0)
    mean.shape = (1,mean.shape[0])
    H = H-mean #takeout column mean
    mean = H.mean(axis=1)
    mean.shape = (mean.shape[0], 1)
    H = H-mean  #takeout raw mean  


    lda(H,5) #compute LDA with keeping the first 5 SVD component


    




    
