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

# Doing the sine wave analysis
def V_sine(V_o,V_a,t_d,theta,f,t1):
    V_t1 = 0
    if(t1 < t_d):
        V_t1 = V_o
    else:
        V_t1 = V_o + V_a*np.exp(-theta*(t1-t_d))*np.sin(2*np.pi*f*(t1-t_d))
    return V_t1

# Doing the pulse wave analysis
def V_pulse(V1,V2,t1,td,tr,tf,tpw,tper):
    V_t1 = 0
    if(t1 >= 0 and t1 < td):
        V_t1 = V1
    elif(t1 >= td and t1 < (td + tr)):
        V_t1 = V1 + (V2-V1)*(t1-td)/tr
    elif(t1 >= (td + tr) and t1 < (td + tr + tpw)):
        V_t1 = V2
    elif(t1 >= (td + tr + tpw) and t1 < (td + tr + tpw + tf)):
        V_t1 = V2 + (V1 - V2)*(t1-(td+tr+tpw))/tf
    elif(t1 >= (td + tr + tpw + tf) and t1 <= (td + tper)):
        V_t1 = V1
    else:
        t1 = td
    return V_t1, t1

# The amount of iterations for the timestep, the higher the more accurate but uses more computing resources
n = 5001

# Defining the Time and Timestep for Transient Simulation
X1 = np.ones((n,1))
t_start = 0
t_end = 12e-3
t1 = 0
h = (t_end - t_start) / (n-1)
t = np.arange(t_start, t_end + h, h)
simul = 4

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
def Is_assigner(node_x,node_y,x,maxi,maxj):
    maxj = 1
    a = matrix(maxi,maxj)
    if(node_x == 0):
        a[node_y-1][0] = x
    elif(node_y == 0):
        a[node_x-1][0] = -x
    else:
        a[node_x-1][maxj-1] = x
        a[node_y-1][maxj-1] = -x
    return a

# A function that initiliaze the resistor values and input into the LHS matrices based on MNA stamps
def R_assigner(node_x,node_y,x,maxi,maxj):
    a = matrix(maxi,maxj)
    if(node_x == 0):
        a[node_y-1,node_y-1] = x
    elif(node_y == 0):
        a[node_x-1,node_x-1] = x
    else:
        a[node_x-1,node_x-1] = x
        a[node_x-1,node_y-1] = -x
        a[node_y-1,node_x-1] = -x
        a[node_y-1,node_y-1] = x
    return a

def Diode_assigner(node_x,node_y,Is,VT,maxi,maxj):
    v1 = 0
    v2 = 0
    x = Is*(np.exp((v1-v2)/VT)-1)
    a = matrix(maxi,maxj)
    if(node_x == 0):
        a[node_y-1,node_y-1] = x
    elif(node_y == 0):
        a[node_x-1,node_x-1] = x
    else:
        a[node_x-1,node_x-1] = x
        a[node_x-1,node_y-1] = -x
        a[node_y-1,node_x-1] = -x
        a[node_y-1,node_y-1] = x
    return a,v1,v2

# A function that adds in voltage for the both LHS and RHS based on MNA stamps
def Vs_assigner(RHS,V_value,G_mat, node_y, node_x):
    Value = np.array([[V_value]])
    # Extending the branch at the LHS matrix
    New_LHS = branch_ext(G_mat,node_y,node_x)
    # Assigning the value at RHS
    New_RHS = np.concatenate((RHS, Value), axis=0)
    size_RHS = np.shape(New_RHS)
    return New_RHS, New_LHS, size_RHS[0]

# A function that extends the matrix with the values needed for the LHS based on MNA stamps
def branch_ext(M,node_x,node_y):
    # Create the column for voltage stamp
    size_a = np.shape(M)
    va= np.zeros((size_a[0],1))
    if(node_x == 0):
        va[node_y-1,0] = 1
    elif(node_y == 0):
        va[node_x-1,0] = -1
    else:
        va[node_x-1,0] = 1
        va[node_y-1,0] = -1
    
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

def Diode_system(LHS, RHS, node_x, node_y, Is, VT, maxi,maxj):
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
    
    x = (Is/VT)*(np.exp((v[node_x-1]-v[node_y-1])/VT))
    
    x1 = x*(v[node_x-1]-v[node_y-1])-Is*(np.exp((v[node_x-1]-v[node_y-1])/VT)-1)
    
    b = matrix(maxi+1,1)
    if(node_x == 0):
        b[node_y-1][0] = x1
    elif(node_y == 0):
        b[node_x-1][0] = -x1
    else:
        b[node_x-1][0] = x1
        b[node_y-1][0] = -x1
    G_x = RHS + b

    a = matrix(maxi+1,maxj+1)
    if(node_x == 0):
        a[node_y-1,node_y-1] = x
    elif(node_y == 0):
        a[node_x-1,node_x-1] = x
    else:
        a[node_x-1,node_x-1] = x
        a[node_x-1,node_y-1] = -x
        a[node_y-1,node_x-1] = -x
        a[node_y-1,node_y-1] = x
    J_x = LHS + a
    
    while error > eps and iteration_counter < 5000:
        
        delta = np.linalg.solve(J_x, np.matmul(J_x,v) - G_x)
        error = np.max(np.abs(delta))
        v = v - delta
        
        x = (Is/VT)*(np.exp((v[node_x-1]-v[node_y-1])/VT))
        
        x1 = x*(v[node_x-1]-v[node_y-1])-Is*(np.exp((v[node_x-1]-v[node_y-1])/VT)-1)
        if(node_x == 0):
            b[node_y-1][0] = x1
        elif(node_y == 0):
            b[node_x-1][0] = -x1
        else:
            b[node_x-1][0] = x1
            b[node_y-1][0] = -x1
        G_x = RHS + b
        
        if(node_x == 0):
            a[node_y-1,node_y-1] = x
        elif(node_y == 0):
            a[node_x-1,node_x-1] = x
        else:
            a[node_x-1,node_x-1] = x
            a[node_x-1,node_y-1] = -x
            a[node_y-1,node_x-1] = -x
            a[node_y-1,node_y-1] = x
        J_x = LHS + a
        iteration_counter += 1
        

    # Here, either a solution is found, or too many iterations
    if error > eps:
        iteration_counter = -1
    return v, J_x, G_x


def plu_solve(A, b):
    
    P, L, U = plu(A)
    
    y = forward_substitution(L, np.dot(P, b))
    
    return back_substitution(U, y)

# The variables for the sine wave voltage source
t_d = 1e-3
f = 1e3
theta = 400
V_a = 0.5
V_o = 1

# The variables for the pulse wave voltage source
V1 = 0
V2 = 1
t1 = 0
td = 1e-3
tr = 0.5e-3
tf = 0.2e-3
tpw = 2e-3
tper = 4e-3

#Values of the variables
Is = 3e-9 # A
VT = 0.05 # V
V1, temp = V_pulse(V1,V2,t1,td,tr,tf,tpw,tper)
Vs = [
    V1
]

I = [
    0
]

R = [
    2e3,3e3,4e3,2e3,1e3
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

R_stamp = [
    R_assigner(1,2,cond(R[0]),Maxi,Maxj),
    R_assigner(2,0,cond(R[1]),Maxi,Maxj),
    R_assigner(2,3,cond(R[2]),Maxi,Maxj),
    R_assigner(3,0,cond(R[3]),Maxi,Maxj),
    R_assigner(3,0,cond(R[4]),Maxi,Maxj)
]

# Adding the resistor values from the stamp to create the overall RHS matrix
T_LHS = mat_sum(R_stamp)

# Adding the branch values from different stamp sources (eg. Voltage source, VCCS) which affects both LHS and RHS


# Adding the voltage source
RHS, G_mat, column_val = Vs_assigner(T_RHS, Vs[0], T_LHS, 0, 1)

print(G_mat)
print(RHS)
print(column_val)
# Create the function matrices of the variables that are going to be used
# new_LHS = np.dot(LHS, X(LHS))


# A function that could analyse and use the stamps to create an MNA circuit matrix
# get the start time
#st1 = time.time()
# Solving the MNA matrix by using linear solvers in python
#ans1 = np.linalg.solve(LHS,RHS)
# get the end time
#et1 = time.time()



# st2 = tm.time()
# The function that would calculate the G_matrix and RHS of the diode using Newton-Raphson 
# temp, J_x, G_x = Diode_system(G_mat, RHS, 1, 2, Is, VT, Maxi, Maxj)
# ans, J_x, G_x = Diode_system(J_x, G_x, 1, 3, Is, VT, Maxi, Maxj)
# get the end time
# et2 = tm.time()

simul = 1

volt1 = np.ones((n,1))
volt2 = np.ones((n,1))
volt3 = np.ones((n,1))
current = np.ones((n,1))
# get the start time
st = tm.time()
match simul:
    # 1 does the pulse wave voltage source
    case 1:
        for i in range (0, len(t)):
            t1 = t1 + h
            V_t, t1 = V_pulse(V1,V2,t1,td,tr,tf,tpw,tper)
            v_value = np.array([[V_t]])
            RHS[column_val-1,0] = v_value
            x_solution = np.linalg.solve(G_mat,RHS)
            # RHS = RHS + Vs_assigner(+)
            
            # Storing the final results at the given timestep
            volt1[i] = x_solution[0]
            volt2[i] = x_solution[1]
            volt3[i] = x_solution[2]
            current[i] = x_solution[3]
        
            # get the end time     
        et = tm.time()   
    
        # Plotting the Graph
        
        plt.plot(t, volt1, color = 'red', label = '$nodal voltage 1$')
        plt.plot(t, volt2, color = 'blue', label = '$nodal voltage 2$')
        plt.plot(t, volt3, color = 'green', label = '$nodal voltage 3$')
        plt.plot(t, current, color = 'yellow', label = '$current$')
        plt.xlabel('Time ($s$)')
        plt.ylabel('Variables ($v1$, $v2$, $v3$, current)')
        plt.title('Variation of the Variables ($v1$, $v2$, $v3$, current) with Time\n')
        plt.legend()
        plt.grid(which='major')
        plt.minorticks_on()
        plt.grid(which='minor', alpha=0.2)
    
    # 2 does the sine wave voltage source
    case 2:
        for i in range (0, len(t)):
            t1 = t1 + h
            V_t = V_sine(V_o,V_a,t_d,theta,f,t1)
            # Storing the final results at the given timestep
            X1[i] = V_t
        
            # get the end time     
        et = tm.time()   
        # Plotting the Graph
        plt.plot(t, X1, color = 'red', label = '$V_s$')
        plt.xlabel('Time ($s$)')
        plt.ylabel('Variables ($V_s$)')
        plt.title('Variation of the Variables ($V_s$) with Time\n' 
                    + 'Timestep, $h$ = {}'.format(h))
        plt.legend()
        plt.grid(which='major')
        plt.minorticks_on()
        plt.grid(which='minor', alpha=0.2)
    case _:
        print("Doing the OP analysis")
        print(np.linalg.solve(G_mat,RHS))
        et = tm.time() 


#print(ans1)
# print(ans)
# print(ans[1]-ans[0])

# get the execution time
# elapsed_time1 = et1 - st1
elapsed_time2 = et - st
print('Execution time:', elapsed_time2, 'seconds')

plt.show()