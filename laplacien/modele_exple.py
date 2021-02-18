# *- coding: UTF-8 -*

import numpy as np
import modele_bocf as bocf
import poisson
import edo
import matplotlib.pyplot as plt

# Coefficient de diffusion
D = 1.e-3

# discretisation en espace
Omega = [0,1,0,1]
N = 10
P = (N+1)**2
y1 = np.zeros((4,P)) # u,v,w,s
y0 = np.zeros((4,P)) # u,v,w,s

A = D*poisson.matrix_neumann2D(Omega,N,N) # Matrice de -D*Laplacien(u) avec CL Neumann

# discretisation en temps
h = 0.1 # pas de temps
M = 500 # nombre d'itérations en temps

# point de départ
t0 = 0.
y0[0,:] = 10  # u
y0[1,:] = 1.  # v
y0[2,:] = 1.  # w
y0[3,:] = 0.  # s

# Pour l'affichage on construit les X et Y
dx = (Omega[1]-Omega[0])/N
Y,X = np.mgrid[Omega[0]:Omega[1]+dx:dx, Omega[2]:Omega[3]+dx:dx]
U = np.reshape(y0[0,:],(N+1,N+1))
plt.contourf(X,Y,U, vmin=0,vmax=1)
plt.show()

# Solution approchée par la méthode d'Euler
for n in np.arange(0,M+1):
  y1 = y0 + h*bocf.G(y0)
  y1[0,:] = y1[0,:] + h* (A*y0[0,:])
  y0 = y1

U = np.reshape(y0[0,:],(N+1,N+1))
plt.contourf(X,Y,U, vmin=0,vmax=1)
plt.show()
  
