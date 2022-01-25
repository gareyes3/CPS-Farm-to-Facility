# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 13:23:55 2021

@author: gareyes3
"""

import numpy as np
import math
import statistics
from matplotlib import pyplot as plt
import Funz



totals = []

for i in list(range(1,1000)):        
    time = 6
    
    seg1 = Funz.F_DieOff1()
    seg2 = Funz.F_DieOff2()
    break_p =Funz.Func_NormalTrunc(0.11,3.71, 0.68,0.98)
    
    if time>break_p: 
        Timeseg1 = time-break_p
        Timeseg2 = time-Timeseg1
        Dieoff1 = seg1*Timeseg1
        Dieoff2 = seg2*Timeseg2
        Total = Dieoff1 +Dieoff2
    elif time<break_p:
        Total = time*Dieoff1
    
    totals.append(Total) 


Reduction = -((2/(2.45/24))**0.3)


0.5*6




totals= []
for i in list(range(1,1000)):
    time = 2
    die=Funz.Func_NormalTrunc(-1.04,-0.33, -0.77,0.21)
    total  = time*die
    totals.append(total)
    
plt.hist(totals)
statistics.median(totals)
