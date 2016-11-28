import scipy.io
import mdlsvd 
import matplotlib.pyplot as plt
import numpy 
from mpl_toolkits.mplot3d import Axes3D

data=scipy.io.loadmat('data/lines2d_a.mat')
X=data['Abis']
[u,v,lu,lv,LM]=mdlsvd.bin_svd_mdl(X)
numpy.savetxt('LM.ascii',LM)
numpy.savetxt('u.ascii',u)
numpy.savetxt('v.ascii',v)

fig = plt.figure()
ax  = fig.gca(projection='3d')
[lu,lv]=numpy.meshgrid(lu,lv)
ax.plot_surface(lu,lv,LM)
ply.savefig('L.svg')
plt.show()

plt.figure()
plt.spy(dot(u,v.T))
plt.savefig('spy.svg')
plt.show()
