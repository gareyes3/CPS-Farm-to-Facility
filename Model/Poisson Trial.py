# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 15:20:20 2022

@author: gareyes3
"""

import numpy as np

1000*10**-2

list_1= []
for i in range(10000):

    result = np.random.poisson(10**-2, 1000).sum()
    list_1.append(result)

list_2=np.array(list_1)
list_2.mean()



np.random.poisson(10**reduction, CFU).sum()