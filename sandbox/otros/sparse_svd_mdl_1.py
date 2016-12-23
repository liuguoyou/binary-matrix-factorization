def bin_svd_codelength(A,u,s,v,q):
    """
    Compute the codelength of encoding a low rank binary matrix A
    as E,uq,s,vq where
    uq = [u]_q, u quantized uniformly to a precision of q
    vq = [v]_q
    E = [A - uq * s * vq']_1, where [1] is binary thresholding
    """
    
