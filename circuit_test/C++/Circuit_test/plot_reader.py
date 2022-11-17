import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# reading the csv files using pandas
v1 = pd.read_csv('V1.csv', header=1)
v2 = pd.read_csv('V2.csv', header=1)
v3 = pd.read_csv('V3.csv', header=1)
v4 = pd.read_csv('V4.csv', header=1)
current = pd.read_csv('current.csv', header=1)
time = pd.read_csv('time.csv', header=1)

v1 = v1.apply(pd.to_numeric,downcast='float')
v2 = v2.apply(pd.to_numeric,downcast='float')
v3 = v3.apply(pd.to_numeric,downcast='float')
v4 = v4.apply(pd.to_numeric,downcast='float')
current = current.apply(pd.to_numeric,downcast='float')
time = time.apply(pd.to_numeric,downcast='float')

# plot color
plt.plot(time, v1, color = 'red', label = '$nodal voltage 1$')
plt.plot(time, v2, color = 'blue', label = '$nodal voltage 2$')
plt.plot(time, v3, color = 'green', label = '$nodal voltage 3$')
plt.plot(time, v4, color = 'black', label = '$nodal voltage 4$')
plt.plot(time, current, color = 'yellow', label = '$current$')

# plot settings
plt.xlabel('Time ($s$)')
plt.ylabel('Variables ($v1$, $v2$, $v3$, $v4$, current)')
plt.title('C++ results ($v1$, $v2$, $v3$, $v4$, current) with Time\n')
plt.legend()
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor', alpha=0.2)

plt.show()