#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
print fname
C=np.loadtxt(fname+'.asc')
print C.shape
hE,hD,hA = C[:,1],C[:,2],C[:,3]
LE,LD,LA,LX = C[:,4],C[:,5],C[:,6],C[:,7]
plt.figure(1)
plt.plot(C[:,1:4],'*-')
plt.grid(True)
plt.legend(('h(E)','h(A)','h(D)','h(X)'))
plt.xlabel('iteration')
plt.ylabel('Hamming weight')
plt.savefig(fname+'-h.pdf',dpi=300)
plt.figure(2)
plt.plot(C[:,4:8],'*-')
plt.grid(True)
plt.legend(('L(E)','L(A)','L(D)','L(X)'))
plt.xlabel('iteration')
plt.ylabel('codelength (bits)')
plt.savefig(fname+'-L.pdf',dpi=300)

