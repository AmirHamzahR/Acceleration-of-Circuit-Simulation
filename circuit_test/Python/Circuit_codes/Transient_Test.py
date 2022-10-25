"""
	Program to test transient simulation
	
	Vc(t) + Ri(t) + Li'(t) = 0
	i'(t) = -Ri(t)/L - 1/L * Vc(t)
	Vc'(t) = i(t)/c
"""
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
def f1(x_old_1, x1, x2, h):
	return x_old_1 + h * (-R/L * x1 - 1/L * x2) - x1

def f2(x_old_2, x1, x2, h):
	return x_old_2 + h/C * x1 - x2


# Defining the Jacobian Function

def jacobian(x_old_1, x_old_2, x1, x2, h):

	J = np.ones((2,2))
	dx = 1e-6
	x1_prime = x1 + dx
	x2_prime = x2 + dx

	# Row 1
	J[0,0] = (f1(x_old_1, x1_prime, x2, h) - f1(x_old_1, x1, x2, h)) / dx
	J[0,1] = (f1(x_old_1, x1, x2_prime, h) - f1(x_old_1, x1, x2, h)) / dx

	# Row 2
	J[1,0] = (f2(x_old_2, x1_prime, x2, h) - f2(x_old_2, x1, x2, h)) / dx
	J[1,1] = (f2(x_old_2, x1, x2_prime, h) - f2(x_old_2, x1, x2, h)) / dx

	return J

# Defining the initial guess values
x1 = 1
x2 = 0
x3 = 0

# Defining the Time and Timestep for Transient Simulation
t_start = 0
t_end = 10 * 0.05
n = 5000
h = (t_end - t_start) / (n-1)
t = np.arange(t_start, t_end + h, h)

# Initializing an array to store the final results at each timestep
X1 = np.ones((n,1))
X2 = np.ones((n,1))

# get the start time
st = tm.time()
# Outer Time Loop

for i in range (0, len(t)):
	
	X_old = np.ones((2,1))
	X_old[0] = x1
	X_old[1] = x2

	X_guess = np.copy(X_old)
	
	F = np.ones((2,1))

	# Parameters for Newton Raphson Solver
	alpha = 1
	tol = 1e-16
	iter = 1

	# Initial error
	error = 9e9

	# Inner Loop - Starting the Newton Raphson Method
	
	while error > tol:
		
		# Computing the Jacobian
		J = jacobian(X_old[0], X_old[1], X_guess[0], X_guess[1],  h)

		# Computing the F vector
		F[0] = f1(X_old[0], X_guess[0], X_guess[1], h)
		F[1] = f2(X_old[1], X_guess[0], X_guess[1], h)

		# Computing new values
		X_new = X_guess - alpha*np.linalg.solve((J), F)

		# Computing the maximum absolute error
		error = np.max(np.abs(X_new - X_guess))

		# Updating old values
		X_guess = X_new

		# Final Results at the given timestep
		x1 = X_new[0]
		x2 = X_new[1]

		# Updating the Iteration
		iter = iter + 1

	# Storing the final results at the given timestep
	X1[i] = x1
	X2[i] = x2
	
# Plotting the Graph

# get the end time
et = tm.time()

plt.plot(t, X1, color = 'red', label = '$x_1$')
plt.plot(t, X2, color = 'blue', label = '$x_2$')
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