# *- coding: UTF-8 -*

import numpy as np

# Le modèle "BOCF" (voir article) sans le terme Laplacien a 4
# fonctions inconnues, y = (u,v,w,s), qui sont solutions d'un système
# d'équations différentielles de la forme
# du/dt = I(u,v,w,s) := - (Jfi+Jso+Jsi)
# dv/dt = g_v(u,v,w,s) := (v_infini(u) - v) / tau_v(u) ou v_infini et tau_v sont des fonctions de u
# dw/dt = g_w(u,v,w,s) := (w_infini(u) - v) / tau_w(u) ou w_infini et tau_w sont des fonctions de u
# ds/dt = g_s(u,v,w,s)

# On va programmer les 4 fonctions I, g_v, g_w, g_s, de telle sorte
# que le système s'écrit de manière abstraite
# dy/dt = G(y) avec G(y) = (I(...), g_v, g_w, g_s)

# Paramètres du modèle BOCF (cf article) -- cas "endo"
u0 = 0
uu = 1.56
theta_v = 0.3
theta_w = 0.13
theta_v_m = 0.2
theta_0 = 0.006
tau_v1_m = 75
tau_v2_m = 10
tau_v_p = 1.4506
tau_w1_m = 6
tau_w2_m = 140
kw_m = 200
uw_m = 0.016
tau_w_p = 280
tau_fi = 0.1
t01 = 470
t02 = 6
tau_so1 = 40
tau_so2 = 1.2
kso = 2
uso = 0.65
tau_s1 = 2.7342
tau_s2 = 2
ks = 2.0994
us = 0.9087
tau_si = 2.9013
tau_w_inf = 0.0273
w_inf_et = 0.78

# Équations des tau
def tau_0(u):
    return t01*(theta_0 > u) + t02*(u >= theta_0)
def tau_so(u):
    return tau_so1 + (tau_so2 - tau_so1)*(1 + np.tanh(kso*(u - uso)))
def tau_vm(u):
    return tau_v1_m*(theta_v_m > u) + tau_v2_m*(u >= theta_v_m)
def tau_wm(u):
    return tau_w1_m + (tau_w2_m - tau_w1_m)*(1 + np.tanh(kw_m*(u -uw_m)))/2.
def tau_s(u):
    return tau_s1*(theta_w > u) + tau_s2*(u >= theta_w)
def v8(u):
    return 1.*(u< theta_v_m)
def w8(u):
    return (1. - u/tau_w_inf)*(theta_0 > u) + w_inf_et*(u >= theta_0)


# Équations de courants J
def Jfi(u,v):
    return v*(u >= theta_v)*(theta_v - u)*(uu - u)/tau_fi
def Jso(u):
    return (u - u0)*(theta_w > u)/tau_0(u) + 1.*(u >= theta_w)*(1/tau_so(u))
def Jsi(u,w,s):
    return - 1.*(u >= theta_w)*w*s/tau_si
def I(u,v,w,s):
    return - ( Jfi(u,v) + Jso(u) + Jsi(u,w,s) )

# Équations v, w, s
def g_v(u,v,w,s):
    return (v8(u) - v)*(theta_v > u)/tau_vm(u) - v*(u >= theta_v)/tau_v_p
def g_w(u,v,w,s):
    return (w8(u) - w)*(theta_w > u)/tau_wm(u) - w*(u >= theta_w)/tau_w_p
def g_s(u,v,w,s):
    return ( 0.5*(1 + np.tanh(ks*(u - us))) - s ) / tau_s(u)

# Fonction G de l'ensemble du système de 4 équations
def G(y):
    # y = (u,v,w,s) et G = (I,g_v,g_w,g_s), attention à l'ordre !  Et
    # je suppose que Y est un tableau Y = Y[0:3], et le résultat G est
    # construit de la même manière, G = G[0:3]
    u = y[0]
    v = y[1]
    w = y[2]
    s = y[3]
    return np.stack((I(u,v,w,s), g_v(u,v,w,s), g_w(u,v,w,s), g_s(u,v,w,s)))
