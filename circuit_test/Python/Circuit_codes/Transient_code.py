"""
This code section can now generate a transient simulation by performing OP analysis which is used for the initial conditions.
The initial conditions is then transferred into the Newton-Raphson solver which solves for the nodal voltages, currents, and
variables that are needed for the transient simulation.
Expansions:
1) Try adding more components such as inductors, VCCS, and increase circuit size
2) Incorporate with C++ to be used for acceleration
"""

import time as tm
import numpy as np
import matplotlib.pyplot as plt

# The amount of iterations for the timestep, the higher the more accurate but uses more computing resources
n = 5001 #5001 seems to be the sweet spot

# Choosing the pulse wave for the transient simulation
simul = 1

# Defining the Time and Timestep for Transient Simulation
X1 = np.ones((n,1))
t_start = 0
t_end = 12e-2
t1 = 0
h = (t_end - t_start) / (n-1)
t = np.arange(t_start, t_end + h, h)

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
    V_t1 = V1
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
def Is_assigner(node_x,node_y,I,maxi,maxj):
    maxj = 1
    a = matrix(maxi,maxj)
    if(node_x == 0):
        a[node_y-1][0] = I
    elif(node_y == 0):
        a[node_x-1][0] = I
    else:
        a[node_x-1][maxj-1] = -I
        a[node_y-1][maxj-1] = I
    return a

# A function that initiliaze the resistor values and input into the LHS matrices based on MNA stamps
def R_assigner(node_x,node_y,R,maxi,maxj):
    a = matrix(maxi,maxj)
    if(node_x == 0):
        a[node_y-1,node_y-1] = R
    elif(node_y == 0):
        a[node_x-1,node_x-1] = R
    else:
        a[node_x-1,node_x-1] = R
        a[node_x-1,node_y-1] = -R
        a[node_y-1,node_x-1] = -R
        a[node_y-1,node_y-1] = R
    return a

def C_assigner(node_x,node_y,C,LHS,RHS,init,h):
    # this if else statment uses trapezoidal formula
    size_LHS = np.shape(LHS)
    x = C/h
    if(node_x == 0):
        x1 = C*init[node_y-1]/h
    elif(node_y == 0):
        x1 = C*init[node_x-1]/h
    else:
        x1 = C*(init[node_x-1]-init[node_y-1])/h
    
    # Matrix stamp for a capacitor on RHS
    a = Is_assigner(node_x,node_y,x1,size_LHS[0],size_LHS[1])
    # Matrix stamp for a capacitor on LHS
    b = R_assigner(node_x,node_y,x,size_LHS[0],size_LHS[1]) 
    New_LHS = LHS + b
    New_RHS = RHS + a
    
    return New_LHS, New_RHS

def Capacitor_system(C,RHS,LHS,h,init):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    eps = 1e-8
    error = 9e9
    iteration_counter = 0
    solution = init
        
    LHS, RHS = C_assigner(2,0,C[0],LHS,RHS,init,h)
    LHS, RHS = C_assigner(3,0,C[1],LHS,RHS,init,h)
    
    while error > eps and iteration_counter < 5000:
        
        delta = np.linalg.solve(LHS, np.matmul(LHS,solution) - RHS)
        error = np.max(np.abs(delta))
        solution -= delta
        
        iteration_counter += 1

    return solution
    
def Diode_assigner(node_x,node_y,Is,VT,LHS,RHS,solution):
    size_LHS = np.shape(LHS)
    x = (Is/VT)*(np.exp((solution[node_x-1]-solution[node_y-1])/VT))
    
    x1 = x*(solution[node_x-1]-solution[node_y-1])-Is*(np.exp((solution[node_x-1]-solution[node_y-1])/VT)-1)
    
    b = Is_assigner(node_x,node_y,-x1,size_LHS[0],size_LHS[1])
    
    if(node_x == 0):
        x = (Is/VT)*(np.exp((solution[node_y-1])/VT))
        x1 = x*(solution[node_y-1])-Is*(np.exp((solution[node_y-1])/VT)-1)
    elif(node_y == 0):
        x = (Is/VT)*(np.exp((solution[node_x-1])/VT))
        x1 = x*(solution[node_x-1])-Is*(np.exp((solution[node_x-1])/VT)-1)
    else:
        x = (Is/VT)*(np.exp((solution[node_x-1]-solution[node_y-1])/VT))
        x1 = x*(solution[node_x-1]-solution[node_y-1])-Is*(np.exp((solution[node_x-1]-solution[node_y-1])/VT)-1)
        
    # Matrix stamp for a capacitor on RHS
    a = Is_assigner(node_x,node_y,-x1,size_LHS[0],size_LHS[1])
    # Matrix stamp for a capacitor on LHS
    b = R_assigner(node_x,node_y,x,size_LHS[0],size_LHS[1]) 
    
    New_LHS = LHS + b
    New_RHS = RHS + a
    
    return New_LHS, New_RHS

def Diode_system(LHS, RHS):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    eps = 1e-9
    error = 9e9
    iteration_counter = 0
    size_LHS = np.shape(LHS)
    solution = np.ones((size_LHS[0],1))
    
    LHS, RHS = Diode_assigner(2,0,Is,VT,LHS,RHS,solution)
    
    while error > eps and iteration_counter < 5000:
        
        delta = np.linalg.solve(LHS, np.matmul(LHS,solution) - RHS)
        error = np.max(np.abs(delta))
        solution -= delta
        
        iteration_counter += 1
        
    return solution

# A function that adds in voltage for the both LHS and RHS based on MNA stamps
def Vs_assigner(V_value, node_x, node_y, LHS, RHS):
    Value = np.array([[V_value]])
    # Extending the branch at the LHS matrix
    New_LHS = branch_ext(LHS,node_y,node_x)
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

def plu_solve(A, b):
    
    P, L, U = plu(A)
    
    y = forward_substitution(L, np.dot(P, b))
    
    return back_substitution(U, y)

# Values of the variables
# Initial voltage before turning on:
V1 = 2
Vs = [
    V1
]

I = [
    0
]

C = [
    0.4e-6,0.8e-6
]

R = [
    3e3,1e3
]

# total number of nodes and voltage sources that is contained to build the base matrix 
# that contains only zeros which will then be occupied with values from the components
T_nodes = 3

# Size of matrix
Maxi = T_nodes
Maxj = Maxi

# Resistor values in a similar format of SPICE simulators's netlist
I_stamp = [
    # default state (for RHS) - Is_assigner(0,0,I[0],Maxi,Maxj)
    Is_assigner(0,0,0,Maxi,Maxj)
]

# Adding the current values from the stamp to create the overall RHS matrix 
RHS = mat_sum(I_stamp)
R_stamp = [
    # default state (for LHS) - R_assigner(0,0,0,Maxi,Maxj)
    R_assigner(2,1,cond(R[0]),Maxi,Maxj),
    R_assigner(3,2,cond(R[1]),Maxi,Maxj)
]

# Adding the resistor values from the stamp to create the overall RHS matrix
LHS = mat_sum(R_stamp)

# Adding the branch values from different stamp sources (eg. Voltage source, VCCS) which affects both LHS and RHS
RHS, LHS, Vs_locate = Vs_assigner(Vs[0], 1, 0, LHS, RHS)

# OP analysis for the initial conditions for the transient simulation
init = np.linalg.solve(LHS,RHS)

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
            
            # The functions and variables needed for the pulse wave voltage source
            t1 = t1 + h
            V_t, t1 = V_pulse(V1 = V1,V2 = 6,t1 = t1,td = 1e-3,tr = 0.5e-3,tf = 0.2e-3,tpw = 2e-3,tper = 4e-3)
            v_value = np.array([[V_t]])
            RHS[Vs_locate-1,0] = v_value
            
            # Using Newton-Raphson and backward euler stamp to simulate the capacitors in the circuit
            ans = Capacitor_system(C,RHS,LHS,h,init)
            x_solution = ans
            
            # Storing the final results at the given timestep
            volt1[i] = x_solution[0]
            volt2[i] = x_solution[1]
            volt3[i] = x_solution[2]
            current[i] = x_solution[3]
            
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
        # get the end time     
           
    
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
        print(np.linalg.solve(LHS,RHS))
        et = tm.time() 

elapsed_time2 = et - st
print('Execution time:', elapsed_time2, 'seconds')

plt.show()