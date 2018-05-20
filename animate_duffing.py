import duffing
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle

a = 4
b = 0.154
r = 5
w = 1.2199778
y = np.array([0.1, 0.1])
yf = 100
dt = 0.005
t_range = np.arange(0, yf + dt, dt)

sol = duffing.solve(a, b, r, w, y, yf, dt, t_range, True)

fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2, 1, 2)


def fun(x):
    # print(sol[x][0])
    return sol[x][0]


counter = 0
n = 15  # number of "periods" of the spring; purely for aesthetics
side = 1
cursor = ax2.axvline(counter*dt, color="red")


# Get the y coordinates of the points in the spring
def y_pts(n, h):
    ys = []
    cur = h
    for i in range(n):
        ys.append(h)
        h = -h  # it's a zigzag, so the next point has the opposite sign
    return ys

def animate(i):
    global x, y, counter
    y1 = fun(counter)+5  # resting position is at x = 5
    plt.subplot(2, 1, 1)
    ax1.clear()
    plt.xlim(0, 10)
    plt.ylim(-3, 3)
    # plt.plot([x for x in range(500)],[-2.15 for x in range(500)],lw=2,color="black")
    # plt.plot([0, 0], [-3, 3])
    plt.plot([0, 10], [-side/2, -side/2], color="black")    # plot the ground
    plt.plot(np.linspace(0, y1-side/2, n), y_pts(n, 0.4), color="gray")  # plot the spring
    # plt.plot(y1,0,marker="s",ms=40) # plot the mass
    ax1.add_patch(Rectangle((y1-side/2, -side/2), side, side))  # plot the mass
    counter += 30

    # ax2.clear()
    plt.subplot(2, 1, 2)
    # plt.plot(t_range, [item[0] for item in sol])
    # plt.plot([counter*dt, counter*dt],[-5, 5])
    # fig.canvas.draw()
    cursor.set_xdata(counter*dt)    # update the cursor


plt.subplot(2, 1, 2)
plt.plot(t_range, [item[0] for item in sol])
ani = animation.FuncAnimation(fig,animate,interval=30)
plt.show()
