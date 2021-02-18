# *- coding: UTF-8 -*

import numpy as np
import discretisation as dis
import CI
import matplotlib.pyplot as plt
import equations_d0 as ed


U = np.zeros(dis.time_steps) #création d'une matrice colonne de taille "time_steps"
V = np.zeros(dis.time_steps)
S = np.zeros(dis.time_steps)
W = np.zeros(dis.time_steps)

U[0] = dis.u0
V[0] = dis.v0
W[0] = dis.w0
S[0] = dis.s0

for i in range(dis.time_steps-1):
        U[i+1] = U[i]   -   dis.dt*ed.J(U[i],V[i],W[i],S[i])
        V[i+1] = V[i]   +   dis.dt*(1.*((CI.theta_v_m >= U[i])   -   V[i])*(CI.theta_v >= U[i])  -   V[i]*(U[i] >= CI.theta_v)/CI.tv_p)
        W[i+1] = W[i]   +   dis.dt*((ed.w8(U[i])   -   W[i])*(CI.theta_w >= U[i])/CI.tw_inf - W[i]*(U[i] >= CI.theta_w)/CI.tw_inf)
        S[i+1] = S[i]   -   (dis.dt/ed.tau_s(U[i]))*(S[i]   -   (1 + np.tanh(CI.ks*(U[i]-CI.us)))/2)


plt.figure(1)
plt.xlabel('time')
plt.ylabel('U(t)')
plt.plot(dis.t, U)
plt.show()

plt.figure(2)
plt.xlabel('time')
plt.ylabel('V(t)')
plt.plot(dis.t, V)
plt.show()

plt.figure(3)
plt.xlabel('time')
plt.ylabel('W(t)')
plt.plot(dis.t, W)
plt.show()

plt.figure(4)
plt.xlabel('time')
plt.ylabel('S(t)')
plt.plot(dis.t, S)
plt.show()
