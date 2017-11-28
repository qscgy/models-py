from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def spring(y, t, b, c):
    x, v = y
    return [v, -b * v - c * x]


def spring_euler(y0, dt, b, c):
    ys = [y0]
    for t in range(int(10 / dt)):
        y_new = ys[-1] + spring(ys[-1], t, b, c) * dt
        ys.append(y_new)
    return ys

y0 = [-0.5, 0]
t = np.arange(0, 10, 0.1)
sol = odeint(spring, y0, t, args=(0.9, 12))
plt.ylim(-0.6, 0.6)
plt.plot(t, sol[:, 0])
plt.show()
