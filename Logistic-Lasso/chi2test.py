import numpy as np
from scipy.linalg import inv
from scipy.stats import chisquare
#from scipy.stats import chi2_contingency
import random
#from openpyxl import Workbook
import xlrd




if __name__=='__main__':

    #load the dataset
    #H = np.loadtxt('fly_test.txt', delimiter=',')

    gen25 = xlrd.open_workbook("C:\Users\mohammad\Dropbox\Convex Optimization\BackgroundLiterature\Code\logit\uniqueGenotypes_gen25_allImputations.xls")
    MOR = gen25.sheet_by_name("MOR")
    MUL = gen25.sheet_by_name("MUL")
    LEW = gen25.sheet_by_name("LEW")
    BRI = gen25.sheet_by_name("BRI")
    DOB = gen25.sheet_by_name("DOB")
    STU = gen25.sheet_by_name("STU")
    #Sell = gen25.sheet_by_name("Selected")
    #UnSl = gen25.sheet_by_name("Unselected")

    
    Sel=[]
    UnS=[]
    for r in range(MOR.nrows-1):
        Sel.append(MOR.row_values(r+1,8,9)[0])
    for r in range(DOB.nrows-1):
        UnS.append(DOB.row_values(r+1,8,9)[0])

    
  #  f_Obs=np.array([[Sel.count(0), Sel.count(1), Sel.count(2)],\
   #        [UnS.count(0), UnS.count(1), UnS.count(2)]])
    

    #chi2 = chisquare(f_obs=np.array([2.*Sel.count(0)+ Sel.count(1), Sel.count(1)+2.*Sel.count(2)]),\
      #               f_exp=np.array([2.*UnS.count(0)+ UnS.count(1), UnS.count(1)+2.*UnS.count(2)]))


    LOR = np.log(2.*Sel.count(0)+ Sel.count(1)) + np.log(UnS.count(1)+2.*UnS.count(2))\
          - np.log(Sel.count(1)+2.*Sel.count(2)) - np.log(2.*UnS.count(0)+ UnS.count(1))

    print np.exp(LOR)

    #SE_01 = np.sqrt(np.sum(1/f_Obs.astype(np.float64)))
    
