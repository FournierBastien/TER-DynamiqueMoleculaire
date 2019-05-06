# -*- coding: cp1252 -*-
import numpy as np
from time import *

# creation de 5 arrays Numpy 
dimension = 40
f = np.ndarray(shape =(dimension),dtype=np.float)
h = np.ndarray(shape =(dimension),dtype=np.float)
w = np.ndarray(shape =(dimension),dtype=np.float)
e = np.ndarray(shape =(dimension),dtype=np.float)
t = np.ndarray(shape =(dimension),dtype=np.float)
# remplissage des arrays Numpy
h.fill (10.)
w.fill (1000.)
e.fill (290.)
t.fill (300.)

t1 = time()  # recuperation du temps système
i = 0
while (i < 1000000):  # execution de un million de calcul
    #f = h * w * ( e - t ) # calcul
    i = i +1
    print ('temps pour')

print ('temps pour un million de passages:', (time() - t1))
print(f)

