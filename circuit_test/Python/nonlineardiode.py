# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:54:47 2022

@author: Roxas
"""

import numpy as np

Vp = 5 #V
R1 = 1e3 # ohms
R2 = 4e3 # ohms
R3 = 3e3 # ohms
R4 = 2e3 # ohms
I0 = 3e-9 # A
VT = 0.05 # V

tol = 1e-9

def f1(V):
    v1 = V[0]; v2 = V[1]
    return( (v1-Vp)/R1 + v1/R2 + I0*(np.exp((v1-v2)/VT)-1) )

def f2(V):
    v1 = V[0]; v2 = V[1]
    return( (v2-Vp)/R3 + v2/R4 - I0*(np.exp((v1-v2)/VT)-1) )

def j11(V):
    v1 = V[0]; v2 = V[1]
    return 1/R1 + 1/R2 + I0/VT * np.exp((v1-v2)/VT)

def j12(V):
    v1 = V[0]; v2 = V[1]
    return -I0/VT * np.exp((v1-v2)/VT)

def j21(V):
    v1 = V[0]; v2 = V[1]
    return -I0/VT * np.exp((v1-v2)/VT)

def j22(V):
    v1 = V[0]; v2 = V[1]
    return 1/R3 + 1/R4 + I0/VT * np.exp((v1-v2)/VT)

# initial guesses
v1 = 0
v2 = 0
iter = 0
V = np.array( [v1,v2] )
F = np.array( [ f1( V ) , f2( V ) ] )
J = np.array( [ [ j11( V ) , j12( V ) ] , [ j21( V ) , j22( V )] ] )
DV = np.dot( np.linalg.inv(J) , F )
estimate = V - DV
err = np.abs(estimate - V)
while ( err > tol ).any():
    F = np.array( [ f1(estimate) , f2(estimate) ] )
    J = np.array( [ [ j11(estimate) , j12(estimate) ] , [ j21(estimate) , j22(estimate) ] ] )
    DV = np.linalg.solve(J , F) # f(x)/f'(x)
    new_estimate = estimate - DV
    err = np.max(np.abs(new_estimate - estimate))
    estimate = new_estimate
    iter = iter + 1
    print(estimate, iter)
print("percentage error for ngspice V1:{:.6E}" .format((new_estimate[0] - 3.458168)/new_estimate[0]*100))
print("percentage error for ngspice V2:{:.6E}" .format((new_estimate[1] - 2.812747)/new_estimate[1]*100))
print("percentage error for LTspice V1:{:.6E}" .format((new_estimate[0] - 3.45802)/new_estimate[0]*100))
print("percentage error for LTspice V2:{:.6E}" .format((new_estimate[1] - 2.81297)/new_estimate[1]*100))
print("percentage error for ngspice D1:{:.6E}" .format(((new_estimate[0]-new_estimate[1]) - (3.458168 - 2.812747))/(new_estimate[0]-new_estimate[1])*100))
print("percentage error for LTspice D1:{:.6E}" .format(((new_estimate[0]-new_estimate[1]) - (3.45802 - 2.81297))/(new_estimate[0]-new_estimate[1])*100))
print("V1={}V\tV2={}V".format(*new_estimate))
print("Voltage across forward biased diode: {:.4E}V".format(new_estimate[0]-new_estimate[1])) #0.645421