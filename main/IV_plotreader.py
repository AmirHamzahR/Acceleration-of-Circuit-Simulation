import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Nodal voltage function

# VDS STEP VOLTAGE ARRAY FOR X-AXIS
volt = pd.read_csv('VDS_step.csv')
volt = volt.apply(pd.to_numeric,downcast='float')
volt = volt.to_numpy()

# ID ARRAY FOR Y-AXIS
id_MOS = pd.read_csv('id.csv')
id_MOS = id_MOS.apply(pd.to_numeric,downcast='float')
id_MOS = id_MOS.to_numpy()
# PLOT SETTINGS FOR VOLTAGE GRAPHS
plt.plot(volt, id_MOS, color = 'blue', label = '$id (MOSFET)$')

plt.xlabel('VDS ($s$)')
plt.ylabel('id (MOSFET)')
plt.title('C++ results MOSFET characteristics\n')
plt.legend()
plt.grid(which='major')
plt.minorticks_on()
plt.grid(which='minor', alpha=0.2)

plt.show()