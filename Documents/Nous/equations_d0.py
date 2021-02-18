# *- coding: UTF-8 -*

# Dimension spatiale = 0

import numpy as np

import CI

def tau_0(u):
    return CI.t01*(CI.theta_0 >= u) + CI.t02*(u >= CI.theta_0)

def tau_so(u):
    return CI.tso1 + (CI.tso2 - CI.tso1)*(1 + np.tanh(CI.kso*(u - CI.uso)))

def tau_vm(u):
    return CI.tv1_m*(CI.theta_v_m >= u) + CI.tv2_m*(u >= CI.theta_v_m)

def tau_wm(u):
    return CI.tw1_m + (CI.tw2_m - CI.tw1_m)*(1 + np.tanh(CI.kw_m*(u -CI.uw_m)))/2.

def tau_s(u):
    return CI.tsi*(CI.theta_w >= u) + CI.ts2*(u >= CI.theta_w)

def w8(u):
    return (1 - u/CI.tw_inf)*(CI.theta_0 >= u) + CI.w_inf_et*(u >= CI.theta_0)

def Jfi(u,v):
    return v*(u >= CI.theta_v)*(CI.theta_v - u)*(CI.uu - u)/CI.tfi
    
def Jso(u):
    return (u - CI.u0)*(CI.theta_w >= u)/tau_0(u) + 1*(u >= CI.theta_w)*(1/tau_so(u))

def Jsi(u,w,s):
    return - 1.*(u >= CI.theta_w)*w*s/CI.tsi

def J(u,v,w,s):
    return Jfi(u,v) + Jso(u) + Jsi(u,w,s)


