"""
This code section will try to generate MNA matrices for simple linear circuit and solve them
using the linear solver in Python.
Expansions:
1) Try to use LU decomposition to solve the circuit
2) Increase the size of the circuit
"""

import numpy as np

# A function that defines the different stamps for the resistor, voltage source, and current source
def cond(R):
    g = 1/R
    return g

def matrix(maxi,maxj):
    size = (maxi,maxj)
    return np.zeros(size)

def Is_assigner(x,node1,node2,maxi,maxj):
    a = matrix(maxi,maxj)
    if(node1 == 0):
        a[node2-1][0] = x
    elif(node2 == 0):
        a[node1-1][0] = -x
    else:
        a[node1-1][maxj-1] = x
        a[node2-1][maxj-1] = -x
    return a



def R_assigner(x,node1,node2,maxi,maxj):
    a = matrix(maxi,maxj)
    if(node1 == 0):
        a[node2-1,node2-1] = x
    elif(node2 == 0):
        a[node1-1,node1-1] = x
    else:
        a[node1-1,node1-1] = x
        a[node1-1,node2-1] = -x
        a[node2-1,node1-1] = -x
        a[node2-1,node2-1] = x
    return a

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

#Values of the variables
Va = 0
Vb = 0

I1 = 1
I2 = 2
I3 = 3

R1 = 1000
R2 = 2000
R3 = 3000
R4 = 4000
R5 = 5000
R6 = 6000

Ia = Is_assigner(I1,1,3,5,1)
Ib = Is_assigner(I2,2,0,5,1)
Ic = Is_assigner(I3,5,4,5,1)
#Need to make this more generic
RHS = Ia + Ib + Ic
print(RHS)

a = R_assigner(cond(R1),1,2,5,5)
b = R_assigner(cond(R2),1,0,5,5)
c = R_assigner(cond(R3),3,0,5,5)
d = R_assigner(cond(R4),2,4,5,5)
e = R_assigner(cond(R5),4,0,5,5)
f = R_assigner(cond(R6),3,5,5,5)

#Need to make this more generic
LHS = a + b + c + d + e + f
print(LHS)

# A function that could analyse and use the stamps to create an MNA circuit matrix

# Solving the MNA matrix by using linear solvers in python
ans1 = np.linalg.solve(LHS,RHS)
ans2 = plu_solve(LHS,RHS)
print(ans1)
print(ans2)