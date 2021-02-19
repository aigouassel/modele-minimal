# *- coding: UTF-8 -*

#Importation des modules
import numpy as np
import matplotlib.pyplot as plt

#Importation des fichiers qui vont nous être utiles
import CI
import discretisation as dis
import edo_u as edo
import equations_d0 as ed


a = 1


##On définit la taille des matrices U, V, W et S en fonction du pas de temps choisit dans le fichier discretisation.py
n = dis.space_steps-1


##Importation de U du fichier edo_o.py, le fichier réduit au problème du laplacien 
U = a*edo.Uo  
##On définit les matrices V, W et S aléatoirement (grâce à np.random) au temps zéro de taille n
V = a*np.random.rand(dis.space_steps,dis.space_steps) 
W = a*np.random.rand(dis.space_steps,dis.space_steps) 
S = a*np.random.rand(dis.space_steps,dis.space_steps) 


##Condition de Dirichlet: les CI sont insérées dans les matrices V,W et S sur les bords du domaine spatial
V[0,0] = V[n,n] = V[0,n] = V[n,0] = dis.v0
W[0,0] = W[n,n] = W[0,n] = W[n,0] = dis.w0
S[0,0] = S[n,n] = S[0,n] = S[n,0] = dis.s0


##Equations X(t_n+1) = X(t_n) + dt*... obtenues avec les Différences Finies
for i in range(dis.space_steps):
    for j in range(dis.space_steps):
        for t in range(dis.time_steps):
            V[i,j] = V[i,j]   +   dis.dt*(1.*((CI.theta_v_m >= U[i,j])   -   V[i,j])*(CI.theta_v >= U[i,j])  -   V[i,j]*(U[i,j] >= CI.theta_v)/CI.tv_p)
            W[i,j] = W[i,j]   +   dis.dt*((ed.w8(U[i,j])   -   W[i,j])*(CI.theta_w >= U[i,j])/CI.tw_inf - W[i,j]*(U[i,j] >= CI.theta_w)/CI.tw_inf)
            S[i,j] = S[i,j]   -   (dis.dt/ed.tau_s(U[i,j]))*(S[i,j]   -   (1 + np.tanh(CI.ks*(U[i,j]-CI.us)))/2)


##On graphe les matrices obtenues U, V, W, et S au temps t_n
plt.figure(1)
plt.title('U(t)')
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(U,extent=[dis.x0,dis.xmax,dis.y0,dis.ymax])
plt.show()
plt.figure(2)
plt.title('S(t)')
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(V,extent=[dis.x0,dis.xmax,dis.y0,dis.ymax])
plt.show()
plt.figure(3)
plt.title('W(t)')
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(W,extent=[dis.x0,dis.xmax,dis.y0,dis.ymax])
plt.show()
plt.figure(4)
plt.title('S(t)')
plt.xlabel('x')
plt.ylabel('y')
plt.imshow(S,extent=[dis.x0,dis.xmax,dis.y0,dis.ymax])
plt.show()
