import numpy as np
import numpy.random as rnd
import scipy.io as io
import matplotlib.pyplot as plt

def binarize(X):
    return (X > 0.5) + 0

def weight(X):
    bX = binarize(X)
    return sum(bX.flat)

def hamming(A,B):
    dAB = (A!=B)+0
    return sum(dAB.flat)

def bdl1(X,K,lu=0.0):
    """
    Binary Matrix Factorization
    Aims to decompose the binary matrix X  (mxn) into binary matrices 
    U (mxk) and V (kxn)
    The idea is to apply a binarized version of the traditional alternate descent-based
    dictionary learning 

    min ||X-UV||_0 + lu||U||_0 + lv||V||_0

    the algorithm is a binarized version of traditional 2nd order (Newton) updates
    on U and V
    """
    lv = lu
    maxit = 100
    [m,n]=X.shape
    # 0) initialize U
    U = (rnd.rand(m,K) > 0.5)+0
    V = (rnd.rand(K,n) > 0.5)+0
    print '|U|=',weight(U),' |V|=',weight(V),' |X-UV|=',hamming(X,binarize(np.dot(U,V)))
    for it in range(maxit):
        Up = U
        Vp = V
        # 1) update V
        Sv = np.diag(1./np.maximum(np.sum(U,0).flat,1.)) # normalization = count ones in cols of U
        Ev = X - ((np.dot(U,V) > 0)+0)
        Gv = -np.dot(U.T,Ev) + lv*V
        dV = np.dot(Sv,-Gv)
        V  = binarize(V + dV)
        #print 'Ev=',Ev,'\nV=',V,'\nGv=',Gv,'\nHv=',Hv,'\nSv=',Sv,'\ndV=',dV
        # 2) update U
        Su = np.diag(1. / np.maximum(np.sum(V,1).flat,1.))
        Eu = X - ((np.dot(U,V) > 0)+0)
        #Gu = np.dot(X,V.T) - np.dot(U,Hu)  + lv*U
        Gu = -np.dot(Eu,V.T)  + lu*U
	dU = np.dot(-Gu,Su)
        #print 'Eu=',Eu,'\nU=',U,'\nGu=',Gu,'\nHu=',Hu,'\nSu=',Su,'\ndU=',dU
        U  = binarize(U + dU)
        #print 'Vp=',Vp,'V=',V
        #print 'dV=',dV
        #print 'UV=',np.dot(U,V)
        dU = hamming(U,Up)
        dV = hamming(V,Vp)
        print '|U|=',weight(U),' |V|=',weight(V),' |dU|=',dU,' |dV|=',dV, ' |X-UV|=',hamming(X,binarize(np.dot(U,V)))
        if (dU+dV) == 0:
            break
    plt.figure()
    plt.spy(np.dot(U,V)>0 )
    plt.show()
    return [U,V]

def exp0():
    X = np.mat('[1,1,0,0;1,1,0,0;0,0,1,1;0,0,1,1]')
    bdl1(X,2)

def exp1(K,l):
    data= io.loadmat('../data/lines2d_a.mat')
    X = data['Abis']
    P = data['data']
    idx = data['goodIdx']-1
    P = P[:,idx]
    print weight(X)
    [U,V] = bdl1(X,K,l)
    plt.figure()
    x = []
    y = []
    C = []
    for k in range(K):
        Ik = np.nonzero(U[:,k])
        Ik = Ik[0]
        print 'k=',k,'|k|=',len(Ik)
        cm = plt.cm.rainbow
        ck = [float(k)/float(K)]*len(Ik)
        x = x + P[0,Ik].tolist()
        y = y + P[1,Ik].tolist()
        C = C + ck
    plt.scatter(x,y,s=36,c=C,marker='o')
    plt.show()
