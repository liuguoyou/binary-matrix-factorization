from plot_bands import *
import matplotlib.pyplot as plt

for du in range(1):
	for lmi in range(4):
		for mi in range(4):
			name = 'du{0:d}_lmi{1:d}_mi{2:d}'.format(du,lmi,mi)
			D = np.loadtxt('tabla_du{0:d}_lmi{1:d}_mi{2:d}.txt'.format(du,lmi,mi))
			D = D / (289*10000)
			plt.figure()
			maxK = range(1,(D.shape[0]*20)+1,20)
			plot_bands(maxK,D,axis=1)
			plt.grid(True)
			plt.xlabel('number of atoms')
			plt.ylabel('codelength')
			plt.title(('Model selection score for '+name))
			plt.axis((0,maxK[-1],0.4,0.8))
			plt.savefig('{0:s}.png'.format(name))
