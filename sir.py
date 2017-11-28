# By Sam Ehrenstein
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def sir(y, t, b, k):
    s, i, r = y
    return [-b*s*i, b*s*i-k*i, k*i]


pop = 10000.0
y0 = [1, 1/pop, 0]
t = np.linspace(0, 500, 501)
sol = odeint(sir, y0, t, args=(0.3, 1/30))
plt.plot(t, sol[:,0]*pop, markersize=0.5)
plt.plot(t, sol[:,1]*pop, markersize=0.5)
plt.plot(t, sol[:,2]*pop, markersize=0.5)
plt.show()
