# *- coding: UTF-8 -*


##DISCRETISATION DE LESPACE

import numpy as np
from numpy import meshgrid
import matplotlib.pyplot as plt


x0 = 0.; y0 = 1.
xmax = 0.; ymax = 1.
Ne = 151 ##nombre de points voulus
de = 1/(Ne-1.) ##pas d'espace
space_steps = Ne-1 ##nombre d'intervalles obtenus de la discrétisation



##ne pas tenir compte
x = np.linspace(x0,xmax,Ne)
y = np.linspace(y0,ymax,Ne)
vx, vy = meshgrid(x,y,indexing='ij')
##plt.plot(vx,vy,'ro')
##plt.axis([x0,xmax,y0,ymax])
##plt.show()



##DISCRETISATION DU TEMPS

t0 = 0; tmax = 0.5; ##en secondes = 500 ms
Nt = 151 ##nombre de points voulus
dt = 1/(Nt-1.) ##pas de temps
time_steps = Nt - 1##nombre d'intervalles obtenus de la discrétisation

t = np.linspace(t0,tmax,time_steps)


##CONDITIONS INITIALES
u0 = 0.; v0 = 0.01; w0 = 0.01; s0 = 0.
