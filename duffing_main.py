import numpy as np
from duffing import *

# Initialize variables
a = 4
b = 0.154
r = 6.1
w = 1.2199778
y = np.array([0.1, 0.1])
yf = 40000
dt = 0.001
t_range = np.arange(0, yf, dt)

# Solve
sol = solve(a, b, r, w, y, yf, dt, t_range, True)
# Plot
#plot(t_range, sol)
# Plot Poincare section
attractor(sol, w, dt)
# Make bifurcation diagram
#bifurcation(a, b, w, y)