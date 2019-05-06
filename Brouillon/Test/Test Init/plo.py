import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 3, 100)
k = 2*np.pi
w = 2*np.pi
dt = 0.01  

t = 0
for i in range(50):
    y = np.cos(k*x - w*t)
    plt.plot(x, y)
    plt.pause(0.01) # pause avec duree en secondes
    t = t + dt

plt.show()

