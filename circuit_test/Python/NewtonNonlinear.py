# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:51:11 2022

@author: Roxas
"""
import numpy as np
import matplotlib.pyplot as plt

def F(x):
        return np.array(
            [(x[0]-x[1])/1000 + x[2],
             (x[1]-x[0])/1000 + 0.001*(x[1])**3,
             x[0] - 1])

def J(x):
        return np.array(
            [[1/1000, -1/1000, 1],
            [-1/1000, (1+3*(x[1]**2))/1000 + (3*(x[1]**2))/(1000*0.1), 0],
            [1, 0, 0]])

def Newton_system(F, J, x, eps):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    F_value = F(x)
    F_norm = np.linalg.norm(F_value, ord=2)  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < 5000:
        delta = np.linalg.solve(J(x), -F_value)
        x = x + delta
        F_value = F(x)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1
        print(x,iteration_counter)

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    return x, iteration_counter


expected = np.array([1, 0.7, -0.003])
#tol = 1e-8
ans, n = Newton_system(F, J, x=np.array([1, 0, -0.001]), eps=1e-16)
#print(n, ans)
error_norm = np.linalg.norm(expected - ans, ord=2)
#assert error_norm < tol, 'norm of error =%g' % error_norm
print('norm of error =%g' % error_norm)
