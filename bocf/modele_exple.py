# *- coding: UTF-8 -*

import numpy as np
import modele_bocf as bocf
import edo
import matplotlib.pyplot as plt

# point de départ
t0 = 0.
y0 = np.array([0,1.,1.,0.]) # u,v,w,s

# Solution approchée par la méthode d'Euler
t,z = edo.euler(t0,y0, bocf.G, 0.01, 20000)

# COURBES
plt.figure(1)

plt.subplot(221)
plt.plot(t,z[0,:])
plt.title('u')
plt.grid(True)

plt.subplot(222)
plt.plot(t,z[1,:])
plt.title('v')
plt.grid(True)

plt.subplot(223)
plt.plot(t,z[2,:])
plt.title('w')
plt.grid(True)

plt.subplot(224)
plt.plot(t,z[3,:])
plt.title('s')
plt.grid(True)

plt.show()
