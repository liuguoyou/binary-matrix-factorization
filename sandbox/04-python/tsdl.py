import numpy as np
import scipy.io as io
import numpy.random as rnd
import matplotlib.pyplot as plt
import numpy.linalg as npla
import universal

def tsdl_fista(X,K,lu,lv):
    """
    Twice Sparse Dictionary Learning
    find local solution to
    1/2||X-UV||_F^2 + lu*||U||_1 + lv*||V||_1
    """
    [m,n]=X.shape
    U = rnd.rand(m,K)
    V = rnd.rand(K,n)
    Up = np.copy(U)
    Vp = np.copy(V)
    #
    # first iteration using Ridge
    #
    for J in range(100):
        np.copyto(Up,U)
        U = np.dot(np.dot(X,V.T), npla.inv(np.dot(V,V.T) + lu*np.eye(K,K)))
        np.copyto(Vp,V)
        V = np.dot(npla.inv(np.dot(U.T,U) + lv*np.eye(K,K)),np.dot(U.T,X))
        f = npla.norm(X - np.dot(U,V))
        f = 0.5*f*f
        ru = npla.norm(U)
        ru = 0.5*lu*ru*ru
        rv = npla.norm(V)
        rv = 0.5*lv*rv*rv
        G = f + ru + rv # surrogate function being solved
        F = f + lu*np.sum(np.abs(U.flat)) + lv*np.sum(np.abs(V.flat)) # target function 
        dU = npla.norm(U - Up)
        dV = npla.norm(V - Vp)
        if (J % 20)==0:
            print 'J=',J,' |dU|=',dU,' |dV|=',dV,' G=',G,' F=',F
        if (dU+dV)/(npla.norm(U)+npla.norm(V)) < 1e-5:
            break
    #
    # main loop: alternate minimization between U and V
    #
    Yv = np.copy(V)
    Yu = np.copy(U)
    dU = 1.1*dU
    dV = 1.1*dV
    for J in range(1000):
        #
        # 1) solve for V given fixed U
        #
        L  = 1.0
        b  = 4.0
        np.copyto(Yv,V)
        t  = 1
        Fp = 1.1*F
        F = npla.norm(X-np.dot(U,V))
        F = 0.5*F*F + lu*np.sum(np.abs(U).flat) + lv*np.sum(np.abs(V).flat)
        if (J % 20)==0:
            print 'J=',J,' |dU|=',dU,' |dV|=',dV,' F(U,V)=',F
        if (dU+dV) < 1e-7:
#            print 'converged in argument'
            break
        elif (abs(Fp-F)/Fp)< 1e-5:
#            print 'converged in cost'
            break
        for fi in range(1000):
            #F = npla.norm(X-np.dot(U,V))
            #F = 0.5*F*F + lu*np.sum(np.abs(U).flat) + lv*np.sum(np.abs(V).flat)
            #print 'fi=',fi,' F(U,V)=',F
            Ev = X - np.dot(U,Yv)
            Gv = -np.dot(U.T,Ev)
            fv = npla.norm(X - np.dot(U,Yv))
            fv = .5*fv*fv
            for bi in range(14):
                Zv = Yv - (1./L)*Gv
                Pv = np.sign(Zv)*np.maximum(0.0,np.abs(Zv) - (1./L)*lv)
                d = Pv-Yv # ok
                gv = lv*np.sum(np.abs(Pv.flat)) # ok
                nd = npla.norm(d)
                Qv = fv + np.trace(np.dot(d,Gv.T)) + (L/2.)*nd*nd# + gv
                Fv = npla.norm(X - np.dot(U,Pv))
                Fv = 0.5*Fv*Fv# + gv
                if Fv <= Qv:
                    break
                L = L * b
            if L >= 1e8:
                break
            np.copyto(Vp,V)
            np.copyto(V,Pv)
            tp = t
            t = (1. + np.sqrt(1. + 4.*t*t))/2.
            dV = V - Vp
            Yv = V + (tp - 1.)/t*dV
            dV = npla.norm(dV)
            #if (fi % 20)==0:
            #    print 'J=',J,' fi=',fi,' L=',L,' |dV|=',npla.norm(dV),' F(V)=',Fv
        #
        # 2) solve for U given V
        #
        L  = 1.0
        b  = 4.0
        np.copyto(Yu,U)
        t  = 1
        for fi in range(1000):
            Eu = X - np.dot(Yu,V)
            Gu = -np.dot(Eu,V.T)
            fu = npla.norm(X - np.dot(Yu,V)) 
            fu = .5*fu*fu
            for bi in range(14):
                Zu = Yu - (1./L)*Gu
                Pu = np.sign(Zu)*np.maximum(0.0,np.abs(Zu)-(1./L)*lu)
                d = Pu-Yu
                nd = npla.norm(d)
                gu = lu*np.sum(np.abs(Pu.flat))
                Qu = fu + np.trace(np.dot(d.T,Gu)) + (L/2.)*nd*nd# + gu
                Fu = npla.norm(X - np.dot(Pu,V))
                Fu = 0.5*Fu*Fu# + gu
                if Fu <= Qu:
                    break
                L = L * b
            if L >= 1e8:
                break # not enough change
            np.copyto(Up,U)
            np.copyto(U,Pu)
            dU = U - Up
            tp = t
            t = (1. + np.sqrt(1. + 4.*t*t))/2.
            Yu = U + (tp - 1.)/t*dU
            dU = npla.norm(dU)
            #if (J % 10)==0:
            #    print 'J=',J,' fi=',fi,' |dU|=',npla.norm(dU),' F(U)=',Fu
    return [U,V]

def mdl_tsdl(X,Kmin,Kmax):
    [m,n] = X.shape
    lmax = np.sqrt(max(m,n))
    ls = [lmax*.0003,lmax*.001,lmax*.01,lmax*.03]
    qs = [2.,1.,.5,.25,.125]
    #ls = [1.,2.]
    #qs = [1.0,0.5]
    Lmin = 1e30
    for K in range(Kmin,Kmax+1):
        for l in ls:
            [U,V] = tsdl_fista(X,K,l,l)
            for q in qs:
                L = universal.low_rank_L(X,U,V,q)
                print 'K=',K,' l=',l,' q=',q,' L=',L
                if L < Lmin:
                    Lmin = L
                    Umin = U
                    Vmin = V
                    lmin = l
                    Kmin = K
                    qmin = q
    return [Umin,Vmin,lmin,qmin]

def exp1(K,l):
    data = io.loadmat('../data/lines2d_a.mat')
    X = data['Abis']
    P = data['data']
    idx = data['goodIdx']-1
    P = P[:,idx]
    [U,V] = tsdl_fista(X,K,l,l)
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

def exp2():
    data = io.loadmat('../data/lines2d_a.mat')
    X = data['Abis']
    P = data['data']
    idx = data['goodIdx']-1
    P = P[:,idx]
    [U,V,l,q] = mdl_tsdl(X,6,14)
    [m,K] = U.shape
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
