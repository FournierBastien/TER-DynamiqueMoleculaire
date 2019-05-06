from multiprocessing import Process
import numpy as np
import matplotlib.pyplot as plt

def display(x,y):
    plt.plot(x,y)
    plt.show()

if __name__ == '__main__':
    x = np.arange(0, 2*np.pi, 0.1)
    y = np.sin(x)
    p = Process(target=display, args=(x,y))
    p.start()

    test=input('test = ')
    print(test)

    p.join()
