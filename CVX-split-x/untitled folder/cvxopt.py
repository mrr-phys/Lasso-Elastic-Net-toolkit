import numpy as np
import cvxopt
import cvxopt.solvers as solvers
import random
from pylab import plot,show,legend


def l1m(A, b, e=0):
	if e<0:
		print 'l1m: e must be nonnegative'
		return None
	m = A.size[0]
	n = A.size[1]
	if b.size[0]!=m:
		print "A must have the same number of rows as the size of b"
		return None
	if A.typecode=='z' or b.typecode=='z':
		return cl1m_socp(A, m, n, b, e)
	else:
		if e>0:
			return rl1m_socp(A, m, n, b, e)
		else: # e == 0
			return rl1m_lp(A, m, n, b)


def rl1m_lp(A, m, n, b):
    c = cvxopt.matrix(0.0, (2*n,1))
    c[0:n] = 1.0
    c[n:2*n] = -1.0
    B = cvxopt.matrix(0.0, (m,2*n))
    B[:,0:n] = A		
    I = cvxopt.matrix(np.eye(n))
    G = cvxopt.matrix(0.0, (2*n,2*n))
    G[0:n, 0:n] = -I
    G[0:n, n:2*n] = -I
    G[n:2*n, 0:n] = I
    G[n:2*n, n:2*n] = -I
    h = cvxopt.matrix(0.0, (2*n,1))
    sol = solvers.lp(c, G, h, B, b)
    s = sol['x']
    return s[0:n]
    



def norm2(a):
    # return ||a||^{2}
    return a.transpose().dot(a)

def mse(a):
    # return mean square error
    return (1./np.size(a))*norm2(a)

#Main
  
if __name__=='__main__':


    n = 200  # signal dimension
    k = 1   # number of Fourier spikes
    #m = 5   # number of measuerments

    g_iter = 1
    s_iter = 50

    er_cutoff = 1.0e-6

    

    itr = 0
    

    while(k>=1 and k<200):

            m=k+1

            mse_cut = 1.
            
            while (mse_cut > er_cutoff):

                    itr+=1

                    print mse_cut 

                    mse_slist = []

                    A = cvxopt.matrix(np.random.randn(m,n)/np.sqrt(m))

                    for i in range(s_iter):

                            f = cvxopt.matrix(0.0, (n,1))
                            a = np.arange(n)
                            np.random.shuffle(a)

                            for i in range(k):
                                    f[int(a[i])] = np.random.randn()

                            fa = []

                            for i in range(n):

                                    fa.append(f[i])


                            y = A*f                    

                            mse_glist = []

                            for i in range(g_iter):

                                    P = cvxopt.matrix(0.0, (n,n))

                                    P[range(0,n*n,n+1)] = 1.0

                                    x = l1m(A*P, y, 1.0e-8)

                                    g = P*x

                                    ga = []

                                    for i in range(n):
                                            ga.append(g[i])

                                    mse_glist.append(mse(np.array(fa)-np.array(ga)))

                            mse_slist.append(np.sum(mse_glist)/g_iter)
                

                    mse_cut = np.sum(mse_slist)/s_iter
                    
                    m+=1

            with open("Phase.txt", "a") as myfile:
                    myfile.write(str(k/float(n))+' ' +str(m/float(n))+'\n')
                    myfile.close()

            print "hiiiiiiiiiiiiiiiiiii"


            k += 5
            



print "run is finished"


