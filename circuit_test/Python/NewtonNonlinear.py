# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:51:11 2022

@author: Roxas
"""
import numpy as np
import time as tm
import matplotlib.pyplot as plt
def F(x):
        return np.array(
            [(x[0]-x[1])/1000 + x[2],
             (x[1]-x[0])/1000 + 0.001*(x[1])**3,
             x[0] - 1])

def J(x):
        return np.array(
            [[1/1000, -1/1000, 1],
            [-1/1000, (1+(x[1])**2)/1000, 0], # Need to add 3 as coefficient of x[1] to simulate the hard-coded function
            [1, 0, 0]])

def Newton_system(x_guess, eps):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    RHS = np.array([0,0,1])
    x_old = x_guess
    F_x = np.matmul(J(x_old),x_old)-RHS # "Generic" expanded Jacobian input
    error = 9e9
    #F_value = F(x_old)  # Hard-coded input
    F_value = F_x  # "Generic" expanded Jacobian input
    iteration_counter = 0
    while error > eps and iteration_counter < 5000:
        delta = np.linalg.solve(J(x_old), F_value)
        x_new = x_old - delta
        #F_value = F(x_new) # Hard-coded input
        F_value = np.matmul(J(x_new),x_new)-RHS # "Generic" expanded Jacobian input
        error = np.max(np.abs(x_new - x_old))
        x_old = x_new
        iteration_counter += 1
        print(x_new,iteration_counter)

    # Here, either a solution is found, or too many iterations
    if error > eps:
        iteration_counter = -1
    return x_new, iteration_counter


expected = np.array([1, 0.7, -0.003])
st = tm.time()
ans, n = Newton_system(x_guess = np.array([1, 0.6823, -0.003]), eps=1e-9)
et = tm.time()
diff = et - st
print("The amount of time taken is", diff)