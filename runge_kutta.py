import numpy as np
import matplotlib.pyplot as plt


def pop(t, y):
    return np.log(t+1)

def hw3(t, y):
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = -y[0]
    return dydt


def rk4(f, tn, yn, h):
    k1 = f(tn, yn)
    k2 = f(tn + h/2, yn + h/2*k1)
    k3 = f(tn + h/2, yn + h/2*k2)
    k4 = f(tn + h, yn + h*k3)
    return yn + h/6*(k1 + 2*k2 + 2*k3 + k4)


y = np.array([0, 1])
yf = 15
dt = 0.05
ys = [y]
for t in np.arange(0, yf, dt):
    y = rk4(hw3, t, y, dt)
    ys.append(y)

for y in ys:
    print(y[1])

plt.plot(np.linspace(0, yf, yf/dt+1), [item[1] for item in ys])
plt.show()
