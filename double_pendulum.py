import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import numpy as np
from scipy.integrate import odeint, ode


# Get position via the Hamiltonian
def doublepend(y, t, l1, l2, g):
    th1, th2, p1, p2 = y
    # print('{} {}'.format(th1, th2))
    c1 = (p1*p2*np.sin(th1-th2))/(l1*l2*(1+(np.sin(th1-th2)**2)))
    c2 = (l2**2*p1**2+l1**2*2.*p2**2-l1*l2*p1*p2*np.cos(th1-th2))/(2*l1**2*l2**2*(1+(np.sin(th1-th2)**2))**2)*np.sin(2*(th1-th2))
    dydt = np.zeros(4)
    dydt[0] = (l2*p1-l1*p2*np.cos(th1-th2))/(l1**2*l2*(1+(np.sin(th1-th2)**2)))
    dydt[1] = (l1*2.*p2-l2*p1*np.cos(th1-th2))/(l1*l2**2*(1+(np.sin(th1-th2)**2)))
    dydt[2] = -2.*g*l1*np.sin(th1)-c1+c2
    dydt[3] = -g*l2*np.sin(th2)+c1-c2
    return dydt


def angle_to_xy(y, l1, l2):
    return np.array([l1*np.sin(y[0]), -l1*np.cos(y[0]), l1*np.sin(y[0])+l2*np.sin(y[1]), -l1*np.cos(y[0])-l2*np.cos(y[1])])


# m1, m2 = 1 to simplify things
l1 = 2
l2 = 2
g = 9.8
y = np.array([np.pi/2.5, np.pi/2, 0, 0])
yf = 100
dt = 0.002
t_range = np.arange(0, yf, dt)
show_anim = True

if show_anim:
    sol = odeint(doublepend, y, t_range, args=(l1, l2, g))
    th1 = [item[0] for item in sol]
    th2 = [item[1] for item in sol]
    x1 = l1*np.sin(th1)
    y1 = -l1*np.cos(th1)
    x2 = x1+l2*np.sin(th2)
    y2 = y1-l2*np.cos(th2)
    pts = np.array([x1, y1, x2, y2]).transpose()
    last = y    # last is the previous value of y
    print(pts)

fig = plt.figure()
ax1 = fig.add_subplot(111)

counter = 0
step = 10
xs = []
ys = []


def rk4(yn, tn, h):
    k1 = doublepend(yn, tn, l1, l2, g)
    k2 = doublepend(yn + h/2*k1, tn + h/2, l1, l2, g)
    k3 = doublepend(yn + h/2*k2, tn + h/2, l1, l2, g)
    k4 = doublepend(yn + h*k3, tn + h, l1, l2, g)
    return yn + h/6*(k1 + 2*k2 + 2*k3 + k4)


def animate(i):
    global counter, step, last, dt, y
    # ends = pts[counter]
    ends = angle_to_xy(y, l1, l2)
    plt.subplot(111)
    p1 = ends[0:2]
    p2 = ends[2:]
    xs.append(p2[0])
    ys.append(p2[1])

    ax1.clear()
    plt.xlim(-4, 4)
    plt.ylim(-7, 1)
    ax1.add_patch(Circle((p1[0], p1[1]), radius=0.15))
    ax1.add_patch(Circle((p2[0], p2[1]), radius=0.15))
    plt.plot([0, p1[0], p2[0]], [0, p1[1], p2[1]], color='black')
    plt.plot(xs, ys, color='red')

    counter += step
    for j in range(step):   # compute where the pendulum will be in the next frame
        y = rk4(last, (counter+j) * dt, dt)
        last = y


if show_anim:
    ani = animation.FuncAnimation(fig, animate, interval=20)
else:
    plt.plot(x1, y1, x2, y2)
plt.show()
