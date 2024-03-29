![AoCs Banner](circuit_test/pics/Banner.png)

# Basic Overview
An MEng project which creates an optimizable circuit simulator mainly in C++ ([Armadillo library](https://arma.sourceforge.net/docs.html)) but also tested in Python ([NumPy](https://numpy.org/doc/), [Matplotlib](https://matplotlib.org/stable/index.html), [Pandas](https://pandas.pydata.org/docs/) library)  and MATLAB. The purpose is to test how far softwares and hardware could go along together in achieving the best performance benchmarks. The main method in simulating these circuits is by using LU decomposition for the sparse matrices from the circuit variables. This is by converting the circuit nodes into matrices using the Modified Nodal Analysis methodology,

![MNApng](circuit_test/pics/MNA.png)

and solving them using Newton-Raphson algorithm and LU decomposition. 

![NRalgo](circuit_test/pics/NRalgo.png)

As the circuit gets bigger, the sparsity of the matrices increases which needs more speed and technique for efficiency. This project sets up a circuit simulator fine-tuned for performance through the utilization of optimization flags on the CPU. Most of the research work has been written down in the [Daybook Loggings](Daybook-loggings.md) markdown file.

# Code Flowchart

Here is the high-level diagram of how the overall circuit code flowchart:

![HighLevel](circuit_test/pics/HighLevel.png)

## Pre-requisite

Before you begin, ensure you have the following libraries and tools installed:

- **C++**: The core implementation of the circuit simulator is done in C++. You'll need a C++ compiler and build tools to compile and run the code.
- **[Armadillo Library](https://arma.sourceforge.net/docs.html)**: Used for linear algebra operations and matrix calculations in C++.
- **Python**: The simulation results are analyzed and visualized using Python scripts.
- **[NumPy](https://numpy.org/doc/)**: A fundamental package for scientific computing with Python.
- **[Matplotlib](https://matplotlib.org/stable/index.html)**: A plotting library for creating static, animated, and interactive visualizations in Python.
- **[Pandas](https://pandas.pydata.org/docs/)**: A data manipulation and analysis library for Python.
- **MATLAB**: Some testing and analysis might be done in MATLAB.

## Installation and Setup

To begin working with this project, follow these steps:

1. Clone the project repository using Git:
   ```bash
   git clone https://github.com/AmirHamzahR/Acceleration-of-Circuit-Simulation.git
   ```

2. Install the required libraries mentioned above as per your selected programming languages (C++, Python, MATLAB).

# Circuit Components and Analysis Instructions

The main cpp file to run the circuit simulator is named as [`Transient_code.cpp`](main/Transient_code.cpp). To assign the circuit components, it follows a similar assignment for any SPICE simulators.  The circuit component functions are called from the header file [`Transient_code.h`](main/Transient_code.h) and assigned inside the `main` and `DynamicNonlinear` function in `Transient_code.cpp`.

### _Resistor Assignment_
Resistor R, node x, node_y, value
```cpp
R_assigner(node_x, node_y, resistance_value, LHS, RHS);
```

### _Capacitor Assignment_
Capacitor C, node x, node y, capacitance value, timestep
```cpp
C_assigner(node_x, node_y, capacitance_value, timestep, LHS, RHS, solution, mode);
```

### _Diode Assignment_
Diode D, node x, node y, Is, VT, cd, timestep
```cpp
Diode_assigner(node_x, node_y, Is, VT, cd, timestep, LHS, RHS, solution, mode);
```

### _Current Source Assignment_
Current Source Is, node x, node y, current value
```cpp
Is_assigner(node_x, node_y, current_value, LHS, RHS);
```

### _Voltage Source Assignment_
Voltage Source Vs, node x, node y, voltage value
```cpp
Vs_assigner(node_x, node_y, voltage_value, LHS, RHS);
```

### _Voltage-Controlled Current Source (VCCS) Assignment_
VCCS, node x, node y, node cx, node cy, transconductance value
```cpp
VCCS_assigner(node_x, node_y, node_cx, node_cy, transconductance, LHS);
```

### _Pulsed Voltage Source Assignment_
Pulsed Voltage Source V, V1, V2, t1, td, tr, tf, tpw, tper, timestep
```cpp
V_pulse(V1, V2, t1, td, tr, tf, tpw, tper, timestep);
```

### _N-MOS Transistor Assignment_
N-MOS Transistor, number, node vd, node vg, node vs, node vb, timestep
```cpp
NMOS_assigner(number, node_vd, node_vg, node_vs, node_vb, timestep, solution, LHS, RHS, mode);
```

### _P-MOS Transistor Assignment_
P-MOS Transistor, number, node vs, node vg, node vd, node vb, timestep
```cpp
PMOS_assigner(number, node_vs, node_vg, node_vd, node_vb, timestep, solution, LHS, RHS, mode);
```

### _Ring Oscillator Assignment_
Ring Oscillator, W, L, R, C, timestep
```cpp
RingOscillatorStages(W, L, R, C, LHS, RHS, solution, timestep, mode);
```

## Transient Simulation Instructions

1. Define `t_start`, `t_end`, and `h` for transient simulation.
2. Configure the number of external nodes, MOSFETs, and cascaded levels.
3. Include "Transient_code.h" to define the main simulation function and matrices.

### Dynamic and Non-Linear Component Assignment
Define stamp matrices for dynamic and non-linear components using `DynamicNonLinear` function.

### Main Circuit Simulation Loop
1. Define matrix size using `Maxi` and `Maxj`.
2. Assign linear components' stamps to `LHS` and `RHS`.
3. Define initial conditions.
4. Perform OP analysis to initialize circuit state.
5. Include transient simulation loop to:
   - Update `LHS` and `RHS`.
   - Call Newton-Raphson system.
   - Store simulation results.

### Save Simulation Results
- `solution.csv`: Nodal voltages and currents.
- `MaxI.csv`: Maximum current values.
- `time.csv`: Time values.

### Plotting Simulation Results with [Transient_plotreader.py](main/Transient_plotreader.py)

1. Save the provided code as `Transientplotreader.py`.
2. Make sure you have `pandas`, `matplotlib`, and `numpy` installed (you can install them using `pip` if not).
3. Run the `Transientplotreader.py` code in your preferred Python environment.

This code reads and processes the CSV data from the circuit simulation and generates plots for nodal voltages and currents. Adjust the plot settings and the voltage/current sources in the code as needed to visualize the simulation results effectively.

### Execution and Benchmarking
- Code supports DC OP Analysis and Transient Analysis.
- Displays execution times for both analyses.

Adjust component values, nodes, time settings, and assignments as needed.

# License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT) and is open for any contributions.

> All documents are referenced from the University of Edinburgh

