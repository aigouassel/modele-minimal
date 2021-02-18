# *- coding: UTF-8 -*

import numpy as np
import matplotlib.pyplot as plt

# clf;

N = 101

# Générer un matrice 2x2 aléatoire avec des coefficients entre -3 et 3
A = np.round(6*np.random.random((2,2))-3)

# Une donnée initiale aléatoire entre -0.5 et 0.5
y0 = np.random.random((2))-0.5

# Trouver les valeurs propres (et les vecteurs propres)
D,U = eig(A)
print(D)

#TRacer plan de phase ne utilisant matplotlib.pyplot.quiver, et np.mgrid

# Définition de la fonction f(y) au second membre de l'équation différentielle y'(t) = f(t,y(t))
def equa_lin(y,t):
    y1 = y[0]
    y2 = y[1]
    return A[0,0]*y1+A[0,1]*y2, A[1,0]*y1+A[1,1]*y2

# Plan de phase sur (-1,1)-1,1)
#y,x = np.mgrid[-1:1:21j,-1:1:21j] # Fabrique les tableaux de coordonnées x,y où on veut tracer un vecteur
#U,V = equa_lin(x,y) # Calcule la valeur de la fonction aux points de la grille

# Résolution approchée de l'équation différentielle
t = np.linspace(0.,1.,101)
y = odeint(equa_lin,np.array([1,1]),t)

# Graphique
plt.plot(t,y[:,0],'r-')
plt.plot(t,y[:,1],'b:')
plt.show()
