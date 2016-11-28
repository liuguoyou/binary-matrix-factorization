import numpy
import numpy.linalg
import numpy.random as rnd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

dot=numpy.dot
norm=numpy.linalg.norm
sign=numpy.sign
maximum=numpy.maximum
sqrt=numpy.sqrt

def norm1(v):
    return numpy.sum(abs(v))

def ssvd(X,lu,lv,u0=numpy.matrix([]),v0=numpy.matrix([])):
    """
    Given a matrix X, and parameter l, find 
    (uo,vo) = \arg\min_{(u,v)} (1/2)||X-uv^t||_F^2 + lu*||u||_1 + lv*||v||_1
    where lu=l*sqrt(m), lv=l*sqrt(n), and (m,n) are the number of rows and 
    cols of X. 
    """
    [m,n] = X.shape
    if (u0.size == 0) | (norm(u0)==0.):
        u = numpy.ones([m,1])*(1./sqrt(m))
        #u = rnd.rand(m,1)*(1./sqrt(m))
    else:
        u = u0
    if (v0.size == 0) | (norm(v0)==0.):
        v = numpy.ones([n,1])*(1./sqrt(n))
        #v = rnd.rand(n,1)*(1./sqrt(n))
    else:
        v = v0
    if lu.size==1:
        lu=lu*numpy.ones(u.shape)
    if lv.size==1:
        lv=lv*numpy.ones(v.shape)

    f0 = norm(X-dot(u,v.T),'fro')
    f0 = 0.5*f0*f0+dot(lu.T,abs(u))+dot(lv.T,abs(v))
    athres = 1e-6
    fthres = 1e-6*f0
    maxit = 1000
    # zeroth iteration
    f = f0
    #print "i=",0," f=",f0
    for i in range(1,maxit+1):
        up = u
        vp = v
        fp = f
        xu = dot(X.T,u)
        #print "nxu=", norm(xu)
        v  = sign(xu) * maximum(0,abs(xu)-lv) * (1.0/(dot(u.T,u)+1e-4))
        #v  = v*(1/maximum(norm(v),1))
        xv = dot(X,v)
        u  = sign(xv) * maximum(0,abs(xv)-lu) * (1/(dot(v.T,v)+1e-4))
        #u  = u*(1/maximum(norm(u),1))
        faux  = norm(X - dot(u,v.T),'fro')
        f = 0.5*faux*faux + dot(lu.T,abs(u)) + dot(lv.T,abs(v))
        [nu,nv] = [norm(u), norm(v)]
        [ndu,ndv] = [norm(u-up),norm(v-vp)]
        df = abs(fp-f)
        #if i % 50 == 0:
            #print "i=",i," f=",f," u=",nu," v=",nv," du=",ndu," dv=",ndv," df=",df
#        if df < fthres:
#            print "converged in function value"
#            break
        if (ndu <= athres*nu) & (ndv <= athres*nv):
            #print "converged in function argument"
            break
    if i==maxit:
        print "did not converge"
    return [u,v]

def ssvdw(X,lu,lv,u0=[],v0=[]):
    [u,v]=ssvd(X,lu,lv,u0,v0)
    lu=1/(abs(u)+numpy.mean(abs(u))+1e-6)
    lv=1/(abs(v)+numpy.mean(abs(v))+1e-6)
    [u,v]=ssvd(X,lu,lv,u,v)
    return [u,v]

def quantize(A,q):
    """ 
    Quantizes a vector/matrix A to precision q.
    The resulting matrix has integer values which, when multiplied by q,
    give the quantized entries of A. 
    The function returns integer so that they can be easily encoded later.
    """
    if q > 0.:
        return numpy.round((1./q)*A)
    else:
        return A;
        
def binarize(A):
    """ 
    binarizes a vector A so that every absolute value below 0.5 is 0,
    and 1 otherwise.
    """
    B=numpy.zeros(A.shape)
    numpy.putmask(B,A>=0.5,1)
    return B

def uni_bernoulli(A):
    """
    universal codelength for an IID Bernoulli sequence .
    Uses the approximate closed form expression (via Stirling's formula)
    of the enumerative code for Bernoulli sequences (Cover'91)
    """
    n = A.size
    k = numpy.count_nonzero(A.flat)
    if n == 0:
        return 0 # null matrix: nothing to encode
    elif k == n:
        return 0.5*numpy.log2(n)
    elif k > 1:
        return 0.5*numpy.log2(n) + \
            -0.5*(numpy.log2(numpy.pi)+1) + \
            (n+0.5)*numpy.log2(n) \
            -(n-k+0.5)*numpy.log2(n-k)-(k+0.5)*numpy.log2(k) # general case
    elif k == 1:
        return 0.5*numpy.log2(n) + numpy.log2(n) # parameter cost + value of parmeter
    else:
        return 0.5*numpy.log2(n) # only parameter cost

def uni_ber_uni(A,Q):
    """
    universal codelength for encoding a matrix with sparse
    entries, where the non-zeros are uniformly distributed in
    the subset of integers {-Q,-Q+1,...,-1,1,...,Q-1,Q} 
    This is naturally described as the entries of A being Bernoulli-Uniform
    variables, which can be described with no loss in two parts: first the Bernoulli
    part using en enumerative code gives the locations of the non-zeros, then the
    non-zeros are described with ceil(log_2(Q))+1 bits each
    """
    Lber = uni_bernoulli(A)
    k = numpy.count_nonzero(A.flat)
    if (k > 0) & (Q > 0):
        Luni = k*(numpy.ceil(numpy.log2(Q))+1)
    else:
        Luni = 0
    return Lber + Luni

def bin_svd_codelength(A,u,v,q):
    """
    Compute the codelength of encoding a low rank binary matrix A
    as E,uq,vq where
    uq = [u]_q, u quantized uniformly to a precision of q
    vq = [v]_q
    E = [A - uq * vq']_1, where [1] is binary thresholding
    """
    uq = quantize(u,q) 
    vq = quantize(v,q)
    Q=max([max(abs(uq)),max(abs(vq))])
    E = binarize(A - (q*uq)*(q*vq.T))
    return [uni_bernoulli(E), uni_ber_uni(uq,Q), uni_ber_uni(vq,Q)] 
    
def bin_svd_opt_codelength(A,u,v):
    qtop=numpy.ceil(numpy.log2(max(A.flat)))
    qs=numpy.power(2.,range(int(qtop)-6,int(qtop)+1))
    Lmin=1e30
    qmin=0
    for q in qs:
        [Le,Lu,Lv]=bin_svd_codelength(A,u,v,q)
        L=Le+Lu+Lv
        #print "q=",q," Le=",Le," Lu=",Lu," Lv=",Lv," L=",L    
        if L < Lmin:
            Lmin = L
            Lumin = Lu
            Lvmin = Lv
            Lemin = Le
            qmin = q
    #print "best q=",qmin," best L=",Lmin
    return [Lmin,Lumin,Lvmin,Lemin,qmin]

def bin_svd_mdl_step(X): 
    [m,n]=X.shape
    norm0=numpy.count_nonzero
    codelength=bin_svd_opt_codelength
    ls=numpy.power(0.5,numpy.arange(1.,8.,1.))
    lu0=numpy.sqrt(m)
    lv0=numpy.sqrt(n)
    lus=lu0*ls
    lvs=lv0*ls
    u=numpy.ones([m,1])/numpy.sqrt(m)
    v=numpy.ones([n,1])/numpy.sqrt(n)
    umin=u
    vmin=v
    Lmin=X.size
    LM=numpy.zeros([ls.size,ls.size])
    i = 0
    for lu in lus:
        j = 0
        for lv in lvs:
            u0=u
            v0=v
            #[u,v] = ssvd(X,lu,lv,u0,v0)
            [u,v] = ssvd(X,lu,lv)
            [L,Lu,Lv,Le,q] = codelength(X,u,v)
            L=int(L)
            print 'lu=',lu,' lv=',lv,' q=',q,' |u|_0=',norm0(u),' |v|_0=',norm0(v), 'L=',L, 'raw=',X.size
            LM[i,j]=L
            if L < Lmin:
                Lmin = L
                Lumin = Lu
                Lvmin = Lv
                Lemin = Le
                umin = u
                vmin = v
                lumin = lu
                lvmin = lv
                qmin = q
            j = j + 1
        i = i + 1
    print 'BEST:lu=',lumin,' lvmin=',lvmin,' q=',qmin,' |u|_0=',norm0(umin),' |v|_0=',norm0(vmin), 'L=',Lmin, 'raw=',X.size
    #fig = plt.figure()
    #ax  = fig.gca(projection='3d')
    #[lumesh,lvmesh]=numpy.meshgrid(lu,lv)
    #ax.plot_surface(lumesh,lvmesh,LM)
    #plt.savefig('LM.svg')
    #plt.show()
    return [umin,vmin,Lmin,Lumin,Lvmin,Lemin]
    
def bin_svd_mdl(X,kmax=0): 
    """
    decompose X as UV', where U and V are matrices whose columns are sparse
    """
    [m,n]=X.shape
    Lu_accu=0 # cost of describing previous u's
    Lv_accu=0 # cost of describing previous v's
    k=-1
    E=X
    Lp=uni_bernoulli(E)
    L = Lp - 1e-3 # trick to enter while at least once...
    if kmax==0:
        kmax = min(m,n)
    while L < Lp:
        if k >= kmax:
            break
        k = k + 1
        print 'k=',k,'L=',L
        Lp = L
        [u,v,L,Lu,Lv,Le] = bin_svd_mdl_step(E)
        if k > 0:
            U = numpy.concatenate((U,u),1)
            V = numpy.concatenate((V,v),1)
        else:
            U = u
            V = v
        Lu_accu = Lu_accu + Lu
        Lv_accu = Lv_accu + Lv
        L = Le + Lu_accu + Lv_accu
        print 'newL=',L
        E = binarize(E-numpy.outer(u,v))
    return [U[:,:k],V[:,:k],Lp]
