import numpy as np
from duffing import solve, plot

a = 4
b = 0.154
r = 2.21
w = 1.22
y = np.array([0.1, 0.1])
yf = 800
dt = 0.005
t_range = np.arange(0, yf + dt, dt)

sol = solve(a, b, r, w, y, yf, dt, t_range, True)
plot(t_range, sol)
