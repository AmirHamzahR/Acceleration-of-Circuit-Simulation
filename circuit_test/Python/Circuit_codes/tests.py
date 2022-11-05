'''Just an area to test out different functions if it works or not'''
'''import numpy as np

def mat_ext(M):
    a = np.shape(M)
    M1=np.zeros((a[0]+1,a[1]+1),M.dtype)
    M1[:-1,:-1]=M
    return M1

def Vs_assigner(M,x):
    Value = np.array([[x]])
    a = np.concatenate((M, Value), axis=0)
    return a

def add_values(M,node1,node2):
    # Create the column for voltage stamp
    size_a = np.shape(M)
    va= np.zeros((size_a[0],1))
    va[node1-1,0] = 1
    va[node2-1,0] = -1
    print(va)
    
    # Create the row for voltage stamp
    zero_ext = np.zeros((1,1))
    ha = va[..., None]
    haz = np.concatenate((ha,zero_ext),axis=None)
    
    #Contatenate both for the new matrix
    M1 = np.hstack((M,va))
    M2 = np.vstack((M1,haz))
    
    return M2

M=np.arange(16).reshape((4,4))    

M1=np.arange(4).reshape((4,1))    

print(Vs_assigner(M1,20))

print(M)
#a = mat_ext(M)
a = add_values(M,2,3)
b = add_values(a,1,2)

print(a)
print(b)'''

import time as tm
import numpy as np
import scipy as sci
from numpy.linalg import inv
import matplotlib.pyplot as plt



L = 1.5 	# Inductance value in Henry
C = 0.0001 	# capacitance in F
R = 100 		# resistance in Ohm, R is considered to be zero
h = 0.001 	# time interval

# Defining the Non-Linear Equations Function
# x1 is the current variable, i, while x2 is the voltage variable, of the capacitor.

n = 5001


# The variables for the sine wave voltage source
t_d = 1e-3
f = 1e3
theta = 400
V_a = 0.5
V_o = 1

# The variables for the pulse wave voltage source
V1 = 0
V2 = 1
td = 1e-3
tr = 0.5e-3
tf = 0.2e-3
tpw = 2e-3
tper = 4e-3

LHS = np.array([[0.000001, 0, 1],
        [0, 1, 0],
        [ 1, 0, 0]])

RHS = np.array([[0],
               [0],
               [1]])
print(np.linalg.solve(LHS,RHS))

# Outer Time Loop

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


# Defining the Time and Timestep for Transient Simulation
X1 = np.ones((n,1))
t_start = 0
t_end = 6e-3
t1 = 0
h = (t_end - t_start) / (n-1)
t = np.arange(t_start, t_end + h, h)
simul = 1

# get the start time
st = tm.time()
match simul:
    # 1 does the pulse wave voltage source
    case 1:
        for i in range (0, len(t)):
            t1 = t1 + h
            V_t, t1 = V_pulse(V1,V2,t1,td,tr,tf,tpw,tper)
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
        et = tm.time() 

	


# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')

plt.show()