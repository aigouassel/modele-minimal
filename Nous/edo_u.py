# *- coding: UTF-8 -*

import numpy as np
import matplotlib.pyplot as plt
import discretisation as dis
import CI

n = dis.space_steps - 1

#random initialization
Uo = np.random.rand(dis.space_steps,dis.space_steps)
Uo[0,0] = Uo[n,n] = dis.u0
Uo[0,n] = Uo[n,0] = dis.u0

K = Uo

#On définit le laplacien avec DF à trois points centrée
def laplacien(Uo):
    Uh = Uo[0:-2,1:-1] # U(i,j+1)'s
    Ub = Uo[2:,1:-1] # U(i,j-1)'s
    Ud = Uo[1:-1,2:] # U(i+1,j)'s
    Ug = Uo[1:-1,0:-2] # U(i-1,j)'s
    Uc = Uo[1:-1,1:-1] # U(i,j)'s
    return (Uh + Ub + Ug + Ud - 4*Uc)/(dis.de**2) #laplacien en différence finie

L = laplacien(Uo)

##Equation u(t_n+1) = u(t_n) + dt*laplacien(u(t_n))
for t in range(dis.time_steps):
  for i in range(1,dis.space_steps-1):
    for j in range(1,dis.space_steps-1):
        Uo[i,j] = dis.dt*L[i-1,j-1]

plt.figure(1)
plt.imshow(Uo, extent=[-1,1,-1,1]);
plt.show()
