# *- coding: UTF-8 -*

import numpy as np


# On programme une fonction qui résout une système différentiel de la
# forme dy/dt = G(y) avec y(0)=y0

def euler(t0,y0,G,h,N):
    """Résolution par la méthode d'Euler explicite de l'équation
    y'(t)=G(y(t)) et y(t0)=y0. Ici on utilise un pas de temps h et on
    calcule N itérations.
    """
    # (y_{n+1} - y_n)/h = G(y_n)
    n_y = np.size(y0)       # Nombre d'inconnues du système d'équations
    y = np.zeros((n_y,N+1)) # Tableau qui va contenir les solutions approchées aux instants t_n = n*h
    t = np.linspace(t0,t0+N*h,N+1)
    y[:,0] = y0             # La première colonne de y est la donnée initiale, y0
    for n in np.arange(0,N): # n vaut successivement 0, 1, 2... N-1
        y[:,n+1] = y[:,n] + h*G(y[:,n])
    # On souhaite renvoyer t et y aux différents instants
    return t, y
