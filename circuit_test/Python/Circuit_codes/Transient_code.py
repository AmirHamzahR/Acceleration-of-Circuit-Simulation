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
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt
from numpy import inf

from numpy import count_nonzero
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

print(h)
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
def Is_assigner(node_x,node_y,I,RHS):
    size_RHS = np.shape(RHS)
    a = matrix(size_RHS[0],size_RHS[1])
    if(node_x == 0):
        a[node_y-1][0] = I
    elif(node_y == 0):
        a[node_x-1][0] = -I
    else:
        a[node_x-1][0] = -I
        a[node_y-1][0] = I
        
    New_RHS = RHS + a
    return New_RHS

# A function that initiliaze the resistor values and input into the LHS matrices based on MNA stamps
def R_assigner(node_x,node_y,R,LHS):
    size_LHS = np.shape(LHS)
    a = np.zeros((size_LHS[0],size_LHS[1]))
    
    if(node_x == 0):
        a[node_y-1][node_y-1] = R
        
    elif(node_y == 0):
        a[node_x-1][node_x-1] = R
    else:
        a[node_x-1][node_x-1] = R
        a[node_x-1][node_y-1] = -R
        a[node_y-1][node_x-1] = -R
        a[node_y-1][node_y-1] = R
        
    New_LHS = LHS + a
    return New_LHS

# A function that initiliaze the resistor values and input into the LHS matrices based on MNA stamps
def VCCS_assigner(node_x,node_y,node_cx,node_cy,R,LHS):
    size_LHS = np.shape(LHS)
    a = matrix(size_LHS[0],size_LHS[1])
    size_LHS = np.shape(LHS)
    a = matrix(size_LHS[0],size_LHS[1])
    if(node_x == 0):
        if(node_cx == 0):
            a[node_y-1,node_cy-1] = R
        elif(node_cy == 0):
            a[node_y-1,node_cx-1] = -R
        else:
            a[node_y-1,node_cx-1] = -R
            a[node_y-1,node_cy-1]= R
    elif(node_y == 0):
        if(node_cx == 0):
            a[node_x-1,node_cy-1] = -R
        elif(node_cy == 0):
            a[node_x-1,node_cx-1] = R
        else:
            a[node_x-1,node_cx-1] = R
            a[node_x-1,node_cy-1] = -R
    else:
        a[node_x-1,node_cx-1] = R
        a[node_x-1,node_cy-1] = -R
        a[node_y-1,node_cx-1] = -R
        a[node_y-1,node_cy-1] = R
    LHS = LHS + a
    
    LHS = LHS + a
    
    return LHS

def N_JFet_assigner(node_d,node_g,node_s,beta,Vto,LAMBDA,LHS,RHS,solution):
    VCCS_assigner(node_d,node_s,node_d,node_s,R,LHS)
    Diode_assigner(node_g,node_d,2.7e-9,0.05,4e-12,h,LHS,RHS,solution)
    Diode_assigner(node_g,node_s,2.7e-9,0.05,4e-12,h,LHS,RHS,solution)
    
    return LHS, RHS

def C_assigner(node_x,node_y,C,LHS,RHS,solution,h,mode):
    x = C/h
    if(mode == 0):
        x1 = 0
        x = 0
    else:
        if(node_x == 0):
            x1 = C*solution[node_y-1]/h
        elif(node_y == 0):
            x1 = C*solution[node_x-1]/h
        else:
            x1 = C*(solution[node_x-1]-solution[node_y-1])/h
    
    # Matrix stamp for a capacitor on RHS
    New_RHS = Is_assigner(node_x,node_y,x1,RHS)
    # Matrix stamp for a capacitor on LHS
    New_LHS = R_assigner(node_x,node_y,x,LHS)
    
    return New_LHS, New_RHS

# A function that extends the matrix with the values needed for the LHS based on MNA stamps
def branch_ext(M,node_x,node_y):
    # Create the column for voltage stamp
    size_a = np.shape(M)
    va= np.zeros((size_a[0],1))
    if(node_x == 0):
        va[node_y-1,0] = 1
    elif(node_y == 0):
        va[node_x-1,0] = 1
    else:
        va[node_x-1,0] = -1
        va[node_y-1,0] = 1
    
    # Create the row for voltage stamp
    zero_ext = np.zeros((1,1))
    ha = va[..., None]
    haz = np.concatenate((ha,zero_ext),axis=None)
    
    #Contatenate both for the new matrix
    M1 = np.hstack((M,va))
    M2 = np.vstack((M1,haz))
    
    return M2

# A function that adds in voltage for the both LHS and RHS based on MNA stamps
def Vs_assigner(V_value, node_x, node_y, LHS, RHS):
    Value = np.array([[V_value]])
    # Extending the branch at the LHS matrix
    New_LHS = branch_ext(LHS,node_x,node_y)
    # Assigning the value at RHS
    if(node_x == 0):
        New_RHS = np.concatenate((RHS, -Value), axis=0)
    else:
        New_RHS = np.concatenate((RHS, Value), axis=0)
    size_RHS = np.shape(New_RHS)
    return New_LHS, New_RHS, size_RHS[0]

def Ind_assigner(L,node_x,node_y,LHS,RHS,h,init):
    if(L == 0):
        # Matrix stamp for an inductor on LHS
        New_LHS = branch_ext(LHS,node_y,node_x)
        size_LHS = np.shape(New_LHS)
        New_LHS[size_LHS[0]-1,size_LHS[1]-1] = -L/h
        Ind_val = np.array([[0]])
        New_RHS = np.concatenate((RHS, Ind_val), axis=0)
    else:
        New_RHS = RHS
        New_LHS = LHS
        size_LHS = np.shape(LHS)
        New_LHS[size_LHS[0]-1,size_LHS[1]-1] = -L/h
        # Matrix stamp for an inductor on RHS
        New_RHS[size_LHS[0]-1,0] = -L*init[size_LHS[0]-1,0]/h
    
    # return New_LHS, New_RHS

def fet_assigner(number, node_vd, node_vg, node_vs, node_vb, solution, LHS, RHS,  mode):
    # the position of the drain, gate, source, and base voltages are hard-coded for the transistor model
    # we are using a discrete MOSFET, so the vb is connected to the source terminal
    add = 0
    for i in range(number-1):
        add +=14
        i += 1
        
    # The settings for the large signal analysis model
    # uses node_vd as starting reference node for the simulation
    if(node_vd == 0):
        vd = 0
    else:
        # LHS = R_assigner(2,5,cond(1e-3),LHS) # RD
        vd = solution[node_vd-1,0]
        
    if(node_vg == 0):
        vg = 0
    else:
        vg = solution[node_vg-1,0]
    
    if(node_vs == 0):
        vs = 0
    else:
        vs = solution[node_vs-1,0]
        
    if(node_vb == 0):
        vb = 0
    else:
        vb = solution[node_vb-1,0]
    # the other settings for fet model based on the large signal analysis
    LHS = R_assigner(1,4,cond(1.1e-3),LHS) # RG
    LHS = R_assigner(7,0,cond(1.2e-3),LHS) # RS
    LHS = R_assigner(6,0,cond(1.3e-3),LHS) # RB
    LHS = R_assigner(5,7,cond(1e9),LHS) # RDS
    LHS, RHS = Diode_assigner(6,5,10e-14,0.05,4e-12,h,LHS,RHS,solution,mode) # Diode BD
    LHS, RHS = Diode_assigner(6,7,10e-14,0.05,4e-12,h,LHS,RHS,solution,mode) # Diode BS
    LHS, RHS = C_assigner(6,5,1e-12,LHS,RHS,solution,h,mode) # Capacitor BD
    LHS, RHS = C_assigner(4,5,1e-12,LHS,RHS,solution,h,mode) # Capacitor GD
    LHS, RHS = C_assigner(4,7,1e-12,LHS,RHS,solution,h,mode) # Capacitor GBO_1
    LHS, RHS = C_assigner(4,6,1e-12,LHS,RHS,solution,h,mode) # Capacitor GBO_2
    LHS, RHS = C_assigner(7,6,1e-12,LHS,RHS,solution,h,mode) # Capacitor BSO
    
    W = 400e-6
    L = 10e-6
    kp = 200e-6
    mCox = kp
    phi = 0.6
    vt0 = 3
    LAMBDA = 0.1
    gamma = 0
    Beta = (mCox)*(W/L)
    vsb = vs - vb
    vt = vt0 + gamma*((np.sqrt(2*phi+vsb)-np.sqrt(2*phi))) # already taking into account the body effect of MOSFETs
    # vt = 2
    vgs = vg - vs
    vds = vd - vs
    # gm = K*vds
    if (vgs>vt) and (vds <= (vgs-vt)): # the transistor is in linear
        # id = gm*vgs + gds*vds
        id = Beta*(vgs-vt-1/2*vds)*vds*(1+LAMBDA*vds)
        gds = Beta*(1+LAMBDA*vds)*(vgs-vt-vds)+Beta*LAMBDA*vds*(vgs-vt-1/2*vds)
        gm = Beta*(1+LAMBDA*vds)*vds
    elif (vgs>vt) and (vds > (vgs-vt)): # the transistor is in saturation
        id = (Beta/2)*np.power((vgs - vt),2) * (1+LAMBDA*vds)
        gds = (Beta/2)*LAMBDA*np.power((vgs-vt),2)
        gm = Beta*(1+LAMBDA*vds)*(vgs-vt)
    else: # the transistor is in cutoff
        id = 0
        gds = 0
        gm = 0
    
    I_DSeq = id - gm*vgs - gds*vds # 10.190 equation 
    # id = gm*vgs + gds*vds
    # print(id)
    # print(vgs,vds,id)
    # print(vgs,vt,vds,id)
    # the position of the drain current is hard-coded for the transistor model
    '''Change the current source, add diodes and capacitors, make it as a switch, IV curve, voltage source as gate, voltage source as source,
    drain voltage, repeating for loop, directly go to analog model on the transistor, look the VT for the diode.'''
    New_RHS = Is_assigner(node_vd,node_vs,I_DSeq,RHS)
    
    LHS = R_assigner(node_vd,node_vs,gds,LHS) # assigning gds
    LHS = VCCS_assigner(node_vd,node_vs,node_vg,node_vs,gm,LHS) # assigning gm
    New_LHS = LHS
    # print(New_RHS)
    # print(New_RHS,"\n",New_LHS)
    return New_LHS,New_RHS

def Diode_assigner(node_x,node_y,Is,VT,cd,h,LHS,RHS,solution,mode):
    if(mode == 0 ): # OP analysis
        a = 0
    else: # Transient analysis
        a = cd/h
    
    if(node_x == 0):
        x = (Is/VT)*(np.exp((-solution[node_y-1])/VT)) + a
        x1 = (x*(-solution[node_y-1])-Is*(np.exp((-solution[node_y-1])/VT)-1))
    elif(node_y == 0):
        x = (Is/VT)*(np.exp((solution[node_x-1]-0)/VT)) + a
        x1 = (x*(solution[node_x-1]-0)-Is*(np.exp((solution[node_x-1]-0)/VT)-1))
    else:
        x = (Is/VT)*(np.exp((solution[node_x-1]-solution[node_y-1])/VT)) + a
        x1 = (x*(solution[node_x-1]-solution[node_y-1])-Is*(np.exp((solution[node_x-1]-solution[node_y-1])/VT)-1))
    
    # x = np.float16(x)
    # x1 = np.float16(x1)
    # Matrix stamp for a diode on RHS
    New_RHS = Is_assigner(node_x,node_y,x1,RHS)
    # Matrix stamp for a diode on LHS
    New_LHS = R_assigner(node_x,node_y,x,LHS)
    
    return New_LHS, New_RHS

# The Newton Raphson system that solves non-linear and dynamic elements in the circuit
def NewtonRaphson_system(RHS,LHS,h,init):
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
    
    # Label the non-linear and dynamic components here
    # LHS, RHS = fet_assigner(1, 1, 0, 3, 3, LHS, RHS)
    # LHS, RHS = Ind_assigner(0.5e-6,2,0,LHS,RHS,h,solution)
    # print(RHS)
      
    while error > eps and iteration_counter < 5:
        J_x, F_x = fet_assigner(1,2, 1, 0, 0, solution, LHS, RHS, 1)
        delta = np.linalg.solve(J_x, np.matmul(J_x,solution) - F_x)
        error = np.max(np.abs(delta))
        solution -= delta
        print(J_x, "\n", F_x)
        iteration_counter += 1
        # print(iteration_counter)
        
    return solution

# Function that linearize the DC components using NR algorithm during operating point, t = 0
def OPanalysis_system(LHS,RHS):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    """
    eps = 1e-8
    error = 9e9
    iteration_counter = 0
    size = np.shape(LHS)
    solution = np.zeros((size[0],1))
    
    # iteration counter set to 50 as followed by the SPICE OPUS textbook
    while error > eps and iteration_counter < 5:
        J_x, F_x = fet_assigner(1, 2, 1, 0, 0, solution, LHS, RHS, 0)
        # print(J_x,"\n",F_x)
        np.savetxt('mat_LHS.csv',LHS, fmt = '%f', delimiter=",")
        # if(sp.issparse(LHS)):
        # delta = spl.spsolve(LHS, (np.matmul(LHS,solution) - RHS))
        # else:
        delta = np.linalg.solve(J_x, np.matmul(J_x,solution) - F_x)
        error = np.max(np.abs(delta))
        solution -= delta
        
        iteration_counter += 1
        
    # saves the LHS, RHS, and solution matrices in CSV files
    np.savetxt('mat_Solution.csv',solution, fmt = '%f', delimiter=",")
    np.savetxt('mat_LHS.csv',J_x, fmt = '%f', delimiter=",")
    np.savetxt('mat_RHS.csv',F_x, fmt = '%f', delimiter=",")
    # print(1.0 - count_nonzero(LHS) / LHS.size) # just to count the non-zeros and see how sparse it is
    
    return solution

''' 
The total number of nodes and voltage sources that is contained to build the base matrix
that contains only zeros which will then be occupied with values from the components.
    - a normal 2 port component will add 1/2 nodes depending on how it is connected
    - each FET transistor adds 9 nodes
'''
T_nodes = 7

# Size of matrix
Maxi = T_nodes
Maxj = Maxi


LHS = np.zeros((Maxi,Maxj)) # LHS matrix
RHS = np.zeros((Maxi,1)) # RHS matrix

# Values of the variables

LHS = R_assigner(3, 2, cond(2e3), LHS)

# Initial voltage before turning on:
V1 = 0

# Adding the branch values from different stamp sources (eg. Voltage source, VCCS) which affects both LHS and RHS
LHS, RHS, temp = Vs_assigner(5, 3, 0, LHS, RHS)
LHS, RHS, Vs_locate = Vs_assigner(V1, 1, 0, LHS, RHS)
# LHS, RHS, Vs_locate = Vs_assigner(Vs[0], 1, 0, LHS, RHS)

# OP analysis for the initial conditions for the transient simulation
solution = OPanalysis_system(LHS,RHS)
# print("LHS:",LHS)
# print("RHS:",RHS)
print("OP analysis results:\n",solution)

volt1 = np.ones((n,1))
volt2 = np.ones((n,1))
volt3 = np.ones((n,1))
volt4 = np.ones((n,1))
current = np.ones((n,1))
RHS1 = np.ones((n,1))
RHS2 = np.ones((n,1))

# get the start time
st = tm.time()
match simul:
    # 1 does the pulse wave voltage source
    case 1:
        for i in range (0, len(t)):
            
            # The functions and variables needed for the pulse wave voltage source
            V_t, t1 = V_pulse(V1 = V1,V2 = 5,t1 = t1,td = 1e-3,tr = 1e-3,tf = 1e-3,tpw = 2e-3,tper = 4e-3)
            v_value = np.array([[V_t]])
            RHS[Vs_locate-1,0] = v_value
            t1 = t1 + h
            # Using Newton-Raphson and backward euler stamp to simulate the capacitors in the circuit
            solution = NewtonRaphson_system(RHS,LHS,h,solution)
            x_solution = solution
            
            # Storing the final results at the given timestep
            volt1[i] = x_solution[0]
            volt2[i] = x_solution[1]
            volt3[i] = x_solution[2]
            # volt4[i] = x_solution[3]
            # current[i] = x_solution[4]
            
        et = tm.time()
        # Plotting the Graph
        
        plt.plot(t, volt1, color = 'red', label = '$nodal voltage 1$')
        plt.plot(t, volt2, color = 'blue', label = '$nodal voltage 2$')
        plt.plot(t, volt3, color = 'green', label = '$nodal voltage 3$')
        # plt.plot(t, volt4, color = 'black', label = '$nodal voltage 4$')
        # plt.plot(t, current, color = 'yellow', label = '$current$')
        # plt.plot(t, RHS2-RHS1, color = 'pink', label = '$resistor current$')
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
