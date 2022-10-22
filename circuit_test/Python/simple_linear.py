"""
This code section will try to generate MNA matrices for simple linear circuit and solve them
using the linear solver in Python.
Expansions:
1) Try to use LU decomposition to solve the circuit
2) Increase the size of the circuit
"""

import time as tm
import numpy as np
import matplotlib.pyplot as plt

# Function for the iterations of a set of unknowns that we are trying to solve
def X(i):
    x = np.zeros((len(i),1))
    return x

# Function to add the matrices within each other after each iteration?



# Functions that defines the different stamps for the resistor, voltage source, and current source
def cond(R):
    g = 1/R
    return g

# A function that creates a matrix of zeros
def matrix(maxi,maxj):
    size = (maxi,maxj)
    return np.zeros(size)

# A function that extends the branch of the matrix
def mat_ext(M):
    a = np.shape(M)
    M1=np.zeros((a[0]+1,a[1]+1),M.dtype)
    M1[:-1,:-1]=M
    return M1

# Function that sums up the matrices for the total overall matrix
def mat_sum(list_of_mat):
    a = list_of_mat[0]
    for i in range(1, len(list_of_mat)):
        a += list_of_mat[i]
    return a

# A function that initiliaze the current sources value and input into the RHS matrices based on MNA stamps
def Is_assigner(node1,node2,x,maxi,maxj):
    maxj = 1
    a = matrix(maxi,maxj)
    if(node1 == 0):
        a[node2-1][0] = x
    elif(node2 == 0):
        a[node1-1][0] = -x
    else:
        a[node1-1][maxj-1] = x
        a[node2-1][maxj-1] = -x
    return a

# A function that initiliaze the resistor values and input into the LHS matrices based on MNA stamps
def R_assigner(node1,node2,x,maxi,maxj):
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

def Diode_assigner(node1,node2,Is,VT,maxi,maxj):
    v1 = 0
    v2 = 0
    x = Is*(np.exp((v1-v2)/VT)-1)
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
    return a,v1,v2

# A function that adds in voltage for the RHS based on MNA stamps
def Vs_assigner(M,x):
    Value = np.array([[x]])
    
    a = np.concatenate((M, Value), axis=0)
    return a

# A function that extends the matrix with the values needed for the LHS based on MNA stamps
def branch_ext(M,node1,node2):
    # Create the column for voltage stamp
    size_a = np.shape(M)
    va= np.zeros((size_a[0],1))
    if(node1 == 0):
        va[node2-1,0] = 1
    elif(node2 == 0):
        va[node1-1,0] = -1
    else:
        va[node1-1,0] = 1
        va[node2-1,0] = -1
    
    # Create the row for voltage stamp
    zero_ext = np.zeros((1,1))
    ha = va[..., None]
    haz = np.concatenate((ha,zero_ext),axis=None)
    
    #Contatenate both for the new matrix
    M1 = np.hstack((M,va))
    M2 = np.vstack((M1,haz))
    
    return M2
    

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
    y = np.zeros_like(b, dtype=np.double)
    
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
    x = np.zeros_like(y, dtype=np.double)

    #Here we perform the back-substitution.  
    #Initializing with the last row.
    x[-1] = y[-1] / U[-1, -1]
    
    #Looping over rows in reverse (from the bottom up), 
    #starting with the second to last row, because the 
    #last row solve was completed in the last step.
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]
        
    return x

def Diode_system(LHS, RHS, node1, node2, Is, VT, maxi,maxj):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    eps = 1e-9
    error = 9e9
    iteration_counter = 0
    v = np.ones((maxi+1,1))
    
    x = (Is/VT)*(np.exp((v[node1-1]-v[node2-1])/VT))
    
    x1 = x*(v[node1-1]-v[node2-1])-Is*(np.exp((v[node1-1]-v[node2-1])/VT)-1)
    
    b = matrix(maxi+1,1)
    if(node1 == 0):
        b[node2-1][0] = x1
    elif(node2 == 0):
        b[node1-1][0] = -x1
    else:
        b[node1-1][0] = x1
        b[node2-1][0] = -x1
    F_x = RHS + b

    a = matrix(maxi+1,maxj+1)
    if(node1 == 0):
        a[node2-1,node2-1] = x
    elif(node2 == 0):
        a[node1-1,node1-1] = x
    else:
        a[node1-1,node1-1] = x
        a[node1-1,node2-1] = -x
        a[node2-1,node1-1] = -x
        a[node2-1,node2-1] = x
    J_x = LHS + a
    
    while error > eps and iteration_counter < 5000:
        
        delta = np.linalg.solve(J_x, np.matmul(J_x,v) - F_x)
        error = np.max(np.abs(delta))
        v = v - delta
        
        x = (Is/VT)*(np.exp((v[node1-1]-v[node2-1])/VT))
        
        x1 = x*(v[node1-1]-v[node2-1])-Is*(np.exp((v[node1-1]-v[node2-1])/VT)-1)
        if(node1 == 0):
            b[node2-1][0] = x1
        elif(node2 == 0):
            b[node1-1][0] = -x1
        else:
            b[node1-1][0] = x1
            b[node2-1][0] = -x1
        F_x = RHS + b
        
        if(node1 == 0):
            a[node2-1,node2-1] = x
        elif(node2 == 0):
            a[node1-1,node1-1] = x
        else:
            a[node1-1,node1-1] = x
            a[node1-1,node2-1] = -x
            a[node2-1,node1-1] = -x
            a[node2-1,node2-1] = x
        J_x = LHS + a
        iteration_counter += 1
        

    # Here, either a solution is found, or too many iterations
    if error > eps:
        iteration_counter = -1
    return v


def plu_solve(A, b):
    
    P, L, U = plu(A)
    
    y = forward_substitution(L, np.dot(P, b))
    
    return back_substitution(U, y)

#Values of the variables
Is = 3e-9 # A
VT = 0.05 # V
Vs = [
    5
]

I = [
    0
]

R = [
    1e3,4e3,3e3,2e3
]

# total number of nodes and voltage sources that is contained to build the base matrix 
# that contains only zeros which will then be occupied with values from the components
T_nodes = 3

# Size of matrix
Maxi = T_nodes
Maxj = Maxi


# Resistor values in a similar format of SPICE simulators's netlist
I_stamp = [
    Is_assigner(0,0,I[0],Maxi,Maxj)
]
# Adding the current values from the stamp to create the overall RHS matrix 
T_RHS = mat_sum(I_stamp)

#Adding the voltage source values
RHS = Vs_assigner(T_RHS,Vs[0])

print(RHS)

R_stamp = [
    R_assigner(1,3,cond(R[0]),Maxi,Maxj),
    R_assigner(1,0,cond(R[1]),Maxi,Maxj),
    R_assigner(2,3,cond(R[2]),Maxi,Maxj),
    R_assigner(2,0,cond(R[3]),Maxi,Maxj)
]

# Adding the resistor values from the stamp to create the overall RHS matrix
T_R = mat_sum(R_stamp)

# Adding the branch values from different stamp sources (eg. Voltage source, VCCS)
LHS = branch_ext(T_R,0,3)
print(LHS)

# Create the function matrices of the variables that are going to be used
# new_LHS = np.dot(LHS, X(LHS))


# A function that could analyse and use the stamps to create an MNA circuit matrix
# get the start time
#st1 = time.time()
# Solving the MNA matrix by using linear solvers in python
#ans1 = np.linalg.solve(LHS,RHS)
# get the end time
#et1 = time.time()

st2 = tm.time()
ans = Diode_system(LHS, RHS, 1, 2, Is, VT, Maxi, Maxj)
# get the end time
et2 = tm.time()

#print(ans1)
print(ans)
print(ans[1]-ans[0])

# get the execution time
# elapsed_time1 = et1 - st1
elapsed_time2 = et2 - st2
print('Execution time:', elapsed_time2, 'seconds')

# plt.plot(t, x[0], color = 'red', label = '$x_1$')
# plt.plot(t, x[1], color = 'blue', label = '$x_2$')
# plt.xlabel('Time ($s$)')
# plt.ylabel('Variables ($x_1$, $x_2$)')
# plt.title('Variation of the Variables ($x_1$, $x_2$) with Time\n')
# plt.legend()
# plt.grid(which='major')
# plt.minorticks_on()
# plt.grid(which='minor', alpha=0.2)

# plt.show()