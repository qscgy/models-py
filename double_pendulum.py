import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, ode

# Get position via the Hamiltonian
def doublepend(y, t, l1, l2, g):
    th1, th2, p1, p2 = y
    c1 = (p1*p2*np.sin(th1-th2))/(l1*l2*(1+np.sin(th1-th2)**2))
    c2 = (l2**2*p1**2+l1**2*2*p2**2-l1*l2*p1*p2*np.cos(th1-th2))/(2*l1**2*l2**2*(1+np.sin(th1-th2)**2)**2)*np.sin(2*(th1-th2))
    dydt = np.zeros(4)
    dydt[0] = (l2*p1-l2*p2*np.cos(th1-th2))/(l1**2*l2*(1+np.sin(th1-th2)**2))
    dydt[1] = (l1*2*p2-l2*p1*np.cos(th1-th2))/(l1*l2**2*(1+np.sin(th1-th2**2)))
    dydt[2] = -2*g*l1*np.sin(th1)-c1+c2
    dydt[3] = -g*l2*np.sin(th2)+c1-c2
    return dydt

# m1, m2 = 1 to simplify things
l1 = 1
l2 = 1.2
g = -9.8
y = np.array([np.pi/6, np.pi/5, 0, 0])
yf = 10
dt = 0.01
t_range = np.arange(0, yf, dt)

sol = odeint(doublepend, y, t_range, args=(l1, l2, g))

th1 = [item[0] for item in sol]
th2 = [item[1] for item in sol]
x1 = l1*np.sin(th1)
y1 = -l1*np.cos(th1)
plt.plot(x1, y1)
plt.plot(x1+l2*np.sin(th2), y1+l2*np.cos(th2))
plt.show()
