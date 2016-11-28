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

def ber_L(A):
    """
    universal codelength for an IID Bernoulli sequence .
    Uses the approximate closed form expression (via Stiring's formula)
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

def ber_uni_L(A,Q):
    """
    universal codelength for encoding a matrix with sparse
    entries, where the non-zeros are uniformly distributed in
    the subset of integers {-Q,-Q+1,...,-1,1,...,Q-1,Q} 
    This is naturally described as the entries of A being Bernoulli-Uniform
    variables, which can be described with no loss in two parts: first the Bernoulli
    part using en enumerative code gives the locations of the non-zeros, then the
    non-zeros are described with ceil(log_2(Q))+1 bits each
    """
    Lber = ber_L(A)
    k = numpy.count_nonzero(A.flat)
    if (k > 0) & (Q > 0):
        Luni = k*(numpy.ceil(numpy.log2(Q))+1)
    else:
        Luni = 0
    return Lber + Luni

def low_rank_L(A,U,V,q):
    """
    Compute the codelength of encoding a low rank binary matrix A
    as E,uq,vq where
    uq = [u]_q, u quantized uniformly to a precision of q
    vq = [v]_q
    E = [A - uq * vq']_1, where [1] is binary thresholding
    """
    [m,K] = U.shape
    Uq = quantize(U,q) 
    Vq = quantize(V,q)
    Qv = numpy.max(Vq)
    Qu = numpy.max(Uq)
    Q = max(Qv,Qu)
    E = binarize(A -numpy.dot(q*Uq,q*Vq))
    L = ber_L(E)
    for k in range(K):
        L = L + ber_uni_L(Uq[:,k],Q) + ber_uni_L(Vq[k,:],Q)
    return L
    
