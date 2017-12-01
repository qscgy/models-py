#Duffing equation project
import numpy as np
import matplotlib
matplotlib.use(("TkAgg"))
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#Takes the current values of the coupled system and returns their derivatives
def duffing(y, t, a, b, r, w):
    x, xp = y
    dydt = np.zeros(2)
    dydt[0] = xp
    dydt[1] = x - b*xp - a*(x**3) +r*np.cos(w*t)
    return dydt

# Custom method to run RK4; much slower than odeint
def rk4(yn, tn, h, a, b, r, w):
    k1 = duffing(yn, tn, a, b, r, w)
    k2 = duffing(yn + h/2*k1, tn + h/2, a, b, r, w)
    k3 = duffing(yn + h/2*k2, tn + h/2, a, b, r, w)
    k4 = duffing(yn + h*k3, tn + h, a, b, r, w)
    return yn + h/6*(k1 + 2*k2 + 2*k3 + k4)


# Solves the Duffing equation for given parameters
def solve(a, b, r, w, y, yf, dt, t_range, use_odeint):
    sol = [y]
    if use_odeint:
        sol = odeint(duffing, y, t_range, args=(a, b, r, w))
    else:
        for t in np.arange(0, yf, dt):
            y = rk4(y, t, dt, a, b, r, w)
            sol.append(y)
    return sol


# Generate bifurcation diagram by iterating through r and plotting the attractors
def bifurcation(a, b, w, y):
    yf = 600
    dt = 2*np.pi/w*0.001
    t_range = np.arange(0, yf + dt, dt)
    r_vals = []
    attractors = []

    for r in np.arange(5.2, 5.9, 0.001):
        print(r)
        sol = solve(a, b, r, w, y, yf, dt, t_range, True)
        sol_x = np.array([item[0] for item in sol]) # Extract x values into a numpy array (for speed)
        ys = np.array([item[1] for item in sol])
        '''
        for i in range(1, sol.shape[0]-1):
            if sol[i-1] < sol[i] > sol[i+1]:
                r_vals.append(r)
                attractors.append(sol[i])
        
        '''
        # Sample every 2pi/w
        for i in np.arange(0, sol_x.shape[0], 1000):
            # print(i)
            if i > sol_x.shape[0] / 2:  # Transient should have dissipated by then
                r_vals.append(r)
                attractors.append(sol_x[int(i)])    # i should be an integer value, but has data type numpy64_float

    plt.plot(r_vals, attractors, 'go', ms=0.5)
    plt.show()


# Plots the solution
def plot(t_range, sol):
    plt.subplot(2, 1, 1)
    plt.plot(t_range, [item[0] for item in sol])
    plt.subplot(2, 1, 2)
    plt.plot([item[0] for item in sol[int(0.75*len(sol)):]], [item[1] for item in sol[int(0.75*len(sol)):]])
    plt.show()


# Plot Poincare section
def attractor(sol, w, dt):
    xs = [item[0] for item in sol]
    ys = [item[1] for item in sol]
    xt = []
    yt = []
    for i in np.arange(0, sol.shape[0], 2*np.pi/w/dt):
        print(i)
        if i > sol.shape[0]/2:  # Steady state behavior should be dominant by then
          xt.append(xs[int(i)])
          yt.append(ys[int(i)])

    print(xt)
    plt.plot(xt, yt, 'ro', ms=1)
    plt.show()
