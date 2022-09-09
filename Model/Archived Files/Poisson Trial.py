# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 15:20:20 2022

@author: gareyes3
"""

import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()
np.random.binomial(1000,10**-2)
print("--- %s seconds ---" % (time.time() - start_time))


start_time = time.time()
np.random.poisson(10**2, 1000).sum()
print("--- %s seconds ---" % (time.time() - start_time))




1000*10**-2

list_1= []
list_2=[]
for i in range(100000):

    result = np.random.poisson(10**-2, 1000).sum()
    list_1.append(result)

list_2=np.array(list_1)
list_2.mean()

plt.hist(list_2)



list_1= []
for i in range(100000):

    result = np.random.binomial(1000,10**-2)
    list_1.append(result)

list_2=np.array(list_1)
list_2.mean()

plt.hist(list_2)




np.random.poisson(10**reduction, CFU).sum()