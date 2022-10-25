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


# Defining the Time and Timestep for Transient Simulation
t_start = 0
t_end = 6e-3
n = 5001
h = (t_end - t_start) / (n-1)
t = np.arange(t_start, t_end + h, h)

# The variables for the sine wave voltage source
X1 = np.ones((n,1))
t_d = 1e-3
f = 1e3
theta = 400
V_a = 0.5
V_o = 1
t1 = 0

# get the start time
st = tm.time()
# Outer Time Loop

for i in range (0, len(t)):
    t1 = t1 + h
    if(t1 < t_d):
        V_t = V_o
    else:
        V_t = V_o + V_a*np.exp(-theta*(t1-t_d))*np.sin(2*np.pi*f*(t1-t_d))
        
    # Storing the final results at the given timestep
    X1[i] = V_t
	
# Plotting the Graph

# get the end time
et = tm.time()

plt.plot(t, X1, color = 'red', label = '$x_1$')
plt.xlabel('Time ($s$)')
plt.ylabel('Variables ($x_1$, $x_2$)')
plt.title('Variation of the Variables ($x_1$, $x_2$) with Time\n' 
			+ 'Timestep, $h$ = {}'.format(h))
plt.legend()
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor', alpha=0.2)




# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')

plt.show()