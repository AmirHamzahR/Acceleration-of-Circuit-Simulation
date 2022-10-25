# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 21:53:18 2022

@author: Roxas
"""

import numpy as np

def plu(A):
    
    #number of rows
    n = A.shape[0]
    
    #Create P, L, and U matrices
    U = A.copy()
    L = np.eye(n, dtype=np.double)
    P = np.eye(n, dtype=np.double)
    
    #Loop between rows
    for i in range(n):
        
        #Permute rows if needed
        for k in range(i, n): 
            if ~np.isclose(U[i, i], 0.0):
                break
            U[[k, k+1]] = U[[k+1, k]]
            P[[k, k+1]] = P[[k+1, k]]
            
        #Eliminate entries below i with row 
        #operations on U and reverse the row 
        #operations to manipulate L
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]
        
    return P, L, U

def forward_substitution(L, b):
    
    #number of rows
    n = L.shape[0]
    
    #Creating the solution vector
    y = np.zeros_like(b, dtype=np.double);
    
    #Here we perform the forward-substitution.  
    #Initializing  with the first row.
    y[0] = b[0] / L[0, 0]
    
    #Looping over rows in reverse (from the bottom  up),
    #starting with the second to last row, because  the 
    #last row solve was completed in the last step.
    for i in range(1, n):
        y[i] = (b[i] - np.dot(L[i,:i], y[:i])) / L[i,i]
        
    return y

def back_substitution(U, y):
    
    #number of rows
    n = U.shape[0]
    
    #Creating the solution vector
    x = np.zeros_like(y, dtype=np.double);

    #Here we perform the back-substitution.  
    #Initializing with the last row.
    x[-1] = y[-1] / U[-1, -1]
    
    #Looping over rows in reverse (from the bottom up), 
    #starting with the second to last row, because the 
    #last row solve was completed in the last step.
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]
        
    return x

def plu_solve(A, b):
    
    P, L, U = plu(A)
    
    y = forward_substitution(L, np.dot(P, b))
    
    return back_substitution(U, y)

A = np.array([ [1/1000, 0, 0, 1], 
                 [0, 1/1000+1/1000, -1/1000, -1], 
                 [0, -1/1000, 1/1000+1/1000, 0], 
                 [1, -1, 0, 0] ])
b = np.array([[0], [0],  [0], [1]])

MNA=plu_solve(A,b)
print(MNA)