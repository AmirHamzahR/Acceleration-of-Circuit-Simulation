import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Nodal voltage function
def nodal_voltage(node_read):
    sol = pd.read_csv('solution.csv')
    sol = sol.apply(pd.to_numeric,downcast='float')
    solution = sol.to_numpy()
    time1 = pd.read_csv('time.csv')
    time1 = time1.apply(pd.to_numeric,downcast='float')
    Maxi = pd.read_csv('MaxI.csv')
    Maxi = Maxi.apply(pd.to_numeric,downcast='float')
    Maxi = int(Maxi.iloc[0][0])
    new_solution = np.zeros(len(time1))
    n_iter = 0
    for i in range(0,(len(time1))):
        
        new_solution[i] = solution[n_iter+(node_read)-1]
        n_iter += Maxi
        
    return new_solution

# Current of resistor function
def R_current(node_x,node_y,R_value):
    R_curr = (nodal_voltage(node_y)-nodal_voltage(node_x))/R_value
    return R_curr

# Current of voltage function
def Vs_current(Vs_num,Total_Vs):
    sol = pd.read_csv('solution.csv')
    sol = sol.apply(pd.to_numeric,downcast='float')
    solution = sol.to_numpy()
    time1 = pd.read_csv('time.csv')
    time1 = time1.apply(pd.to_numeric,downcast='float')
    Maxi = pd.read_csv('MaxI.csv')
    Maxi = Maxi.apply(pd.to_numeric,downcast='float')
    Maxi = int(Maxi.iloc[0][0])
    V_curr = np.zeros(len(time1))
    n_iter = Maxi
    for i in range(0,(len(time1))):
        
        V_curr[i] = solution[n_iter-(Total_Vs-Vs_num)-1]
        n_iter += Maxi
        
    return V_curr

# NODAL VOLTAGES FOR Y-AXIS
v1 = nodal_voltage(1)
v2 = nodal_voltage(2)
v3 = nodal_voltage(3)
v4 = nodal_voltage(4)
v5 = nodal_voltage(5)

# TIME ARRAY FOR X-AXIS
time = pd.read_csv('time.csv')
time = time.apply(pd.to_numeric,downcast='float')
time = time.to_numpy()

# PLOT SETTINGS FOR VOLTAGE GRAPHS
a = plt.figure(1)
plt.plot(time, v1, color = 'red', label = '$nodal voltage 1$')
plt.plot(time, v2, color = 'blue', label = '$nodal voltage 2$')
plt.plot(time, v3, color = 'green', label = '$nodal voltage 3$')
plt.plot(time, v4, color = 'yellow', label = '$nodal voltage 4$')
plt.plot(time, v5, color = 'black', label = '$nodal voltage 5$')

plt.xlabel('Time ($s$)')
plt.ylabel('Variables ($v1$, $v2$, $v3$, $v4$, current)')
plt.title('C++ results ($v1$, $v2$, $v3$, $v4$, current) with Time\n')
plt.legend()
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor', alpha=0.2)

a.show()

# PLOT SETTINGS FOR CURRENT GRAPHS

# # VOLTAGE SOURCE CURRENTS FOR Y-AXIS
# Note that the numbering of the voltage source currents should be the same position
# as how you have assigned the voltage sources in C++ code. V1 = 1, V_pulse = 2, etc.
# Total_Vs is the total voltage sources assigned (This does not use the node numbering system.)
V1_curr = Vs_current(Vs_num = 1,Total_Vs = 2) # V1
V2_curr = Vs_current(Vs_num = 2,Total_Vs = 2) # V_pulse

# RESISTOR CURRENTS FOR Y-AXIS
# R1_curr = R_current(1,2,2000)

b = plt.figure(2)

plt.plot(time, V1_curr, color = 'black', label = '$current voltage 1$')
plt.plot(time, V2_curr, color = 'red', label = '$ voltage pulse$')

plt.xlabel('Time ($s$)')
plt.ylabel('Variables ($V1curr$)')
plt.title('C++ results ($V1curr$) with Time\n')
plt.legend()
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor', alpha=0.2)

b.show()
plt.show()