# *- coding: UTF-8 -*

import numpy as np
import edo
import modele_bocf as bocf
import modele_exple as me
#on importe u afin de pouvoir travailler le laplacien
t, y = edo.euler(t0,y0,G,h,N)
u = y[0,:]

#on crée une matrice ligne de taille -1 comparé à u
l_u = np.zeros((1,N))

#définition du laplacien
for i in np.range(1,n-1):
    l_u[i] = (u[i+1] - 2*u[i] + u[i+1])/((me.h)**2)

def laplacien(t0,y0,G,h,N):
    t, y = edo.euler(t0,y0,G,h,N)
    u = y[0,:]
    u1 = (u[1] - u[0])/h
    for n in range(np.arange(1,N-1)):
           
