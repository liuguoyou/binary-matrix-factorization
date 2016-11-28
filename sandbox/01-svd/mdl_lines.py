import scipy.io
import mdlsvd 
import matplotlib.pyplot as plt
import numpy 
from mpl_toolkits.mplot3d import Axes3D
dot=numpy.dot

data=scipy.io.loadmat('data/lines2d_a.mat')
X=data['Abis']
[U,V,L]=mdlsvd.bin_svd_mdl(X)
