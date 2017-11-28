#Duffing equation project
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def duffing(y, t, a, b, r, w):
    x, xp = y
    dydt = np.zeros(2)
    dydt[0] = xp
    dydt[1] = x - b*xp - a*(x**3) +r*np.cos(w*t)
    return dydt

def rk4(yn, tn, h, a, b, r, w):
    k1 = duffing(yn, tn, a, b, r, w)
    k2 = duffing(yn + h/2*k1, tn + h/2, a, b, r, w)
    k3 = duffing(yn + h/2*k2, tn + h/2, a, b, r, w)
    k4 = duffing(yn + h*k3, tn + h, a, b, r, w)
    return yn + h/6*(k1 + 2*k2 + 2*k3 + k4)


a = 4
b = 0.154
r = 2.21
w = 1.22
y = np.array([2, 0])
yf = 800
ys = [y]
dt = 0.005

for t in np.arange(0, yf, dt):
    y = rk4(y, t, dt, a, b, r, w)
    ys.append(y)

t_range = np.arange(0, yf + dt, dt)

plt.subplot(2, 1, 1)
plt.plot(t_range, [item[0] for item in ys])
plt.subplot(2, 1, 2)
plt.plot([item[0] for item in ys], [item[1] for item in ys])
plt.show()