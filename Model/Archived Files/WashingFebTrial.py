# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 09:07:30 2022

@author: gareyes3
"""
import numpy as np
import pandas as pd
import math


XL = 10000
part_weight = 5


V = 2_500
Blw = np.random.triangular(0.38,0.75,2.2)
alpha = np.random.uniform(0.5,0.75)
Xs = np.random.triangular(0.003,0.055,0.149)
Rpw = 1/25
WT = 1
FCmin = 5


sigma = XL/(5*454)

Xw = 0





#def F_Washing(df,Blw,alpha, Xs, Rpw, ) :
    
df= pd.DataFrame({
    "Cont":np.zeros(100)
        })
    
df.Cont=10000
PartWeight = 5
    
#Contamination Info:
    #Initial water levels
Xw = 0

XL = df.at[1,"Cont"]
sigma = XL/(PartWeight*454)
#Model Parameters
V = 2_500
Blw = np.random.triangular(0.38,0.75,2.2)
alpha = np.random.uniform(0.5,0.75)
Xs = np.random.triangular(0.003,0.055,0.149)
Rpw = 1/55 #3500 L * and 45 kg rate per min
WT = 1
FCmin = 5

#intermediate Outs
Bws = sigma * (1-Xs)* (WT/V)
    #Change in water Conamination
Change_XW = Bws - (Blw * Xw *  Rpw) - (alpha * Xw * FCmin)
if Change_XW>0:
    MaxCFUs = Xw*10**math.ceil(Change_XW)
    Difference =  math.ceil(Change_XW) - Change_XW
    Xw =np.random.binomial(MaxCFUs, 10**-Difference)
else:
    Xw= np.random.binomial(Xw,10**Change_XW)
    
    #Change in XL
Change_XL = (Blw * Xw ) - (alpha * sigma * FCmin)  - (sigma/WT)

if Change_XL>0:
    MaxCFUs = XL*10**math.ceil(Change_XL)
    Difference =  math.ceil(Change_XW) - Change_XW
    XL =np.random.binomial(MaxCFUs, 10**-Difference)
else:
    XL= np.random.binomial(XL,10**Change_XL)


df.at[1,"Cont"] = XL


Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
(Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
((AvCont- AvContAfter)*Rate)/V

