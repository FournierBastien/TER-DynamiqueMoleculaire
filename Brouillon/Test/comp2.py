# -*- coding: cp1252 -*-
import numpy as np
from time import *

# creation de 5 arrays Numpy 
dimension = 40
f = [0] * dimension#np.ndarray(shape =(dimension),dtype=np.float)
h = [10.] * dimension
w = [1000.] * dimension
e = [290.] * dimension
t = [300.] * dimension
# remplissage des arrays Numpy
"""h.fill (10.)
w.fill (1000.)
e.fill (290.)
t.fill (300.)"""


t1 = time()  # recuperation du temps système
for i in range(1000000):  # execution de un million de calcul
	print ('temps pour')
   # f[i] =h[i] * w[i] * ( e[i] - t ) # calcul
 
print ('temps pour un million de passages:', (time() - t1))
print (f)
