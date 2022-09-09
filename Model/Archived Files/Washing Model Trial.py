# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:38:01 2021

@author: gareyes3


"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import math
import random
import itertools
from matplotlib import pyplot as plt


'''
#Setup
Time_Wash=Agg_WashTime
Times = np.arange(0, Time_Wash+0.1, 0.1).tolist()
Times = [round(num, 1) for num in Times]

C= 0 # ini #(mg/L) #Concentration of Free Chrloride available
O = 0  # Initial Oxygen deman
Pre_runningT = 0
RK = 1 #Dosing Period. 
Ro = 12 #Chlorine dosing period, ever 12 minutes
Ro0 = 2 #Minutes duration of dose
SigC = 1.70*(10**-3)
BC = 5.38*(10**-4)
K0 = 32.3 # free chroliride demand per min

r1= 12.75 #(mg/(ml/min**2))
r2 = 7.47 #(mg/(ml/min**2))
r3 = 5.56 #(mg/(ml/min**2))


#Inputs

V = 3200 *1000 #L #From Luo et al 2012. 
Rate = 45 #kg/min #From Luo et al 2012. 
Wash_Time = 2.3 #min 
c1 = 1/Wash_Time #Reciprocal of average time. 
L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
'''

def F_Chloride_lvl (Time_Wash):
    #Function Inputs. 
    #Changing times to 0.1 increments.
    Times = np.arange(0, Time_Wash+0.1, 0.1).tolist()
    Times = [round(num, 1) for num in Times]
    #Addition Rates
    r1= 12.75 #(mg/(ml/min**2))
    r2 = 7.47 #(mg/(ml/min**2))
    r3 = 5.56 #(mg/(ml/min**2))
    #Dose
    Ro = 12 #Chlorine dosing period, ever7 12 minutes
    Ro0 = 2 #Minutes duration of dose
    #Time
    Pre_runningT = 0 #Runing time variable
    K0 = 32.3 # free chrolirine demand per min 
    C= 0 # initial #(mg/L) #Concentration of Free Chrloride available
    O = 0  # Initial Oxygen demand
    #Other parameters
    SigC = 1.70*(10**-3) #Natural decay of FC
    BC = 5.38*(10**-4) #Depletion rate of FC in water. 
    A_Per =0
    List_Time_Ints = list(range(Ro,500,Ro))
    List_C=[]
    for i in Times: 
        Running_Time = i
        if(Running_Time in List_Time_Ints):
            A_Per=A_Per+1
        Time_Interval = Running_Time-(Pre_runningT)
        if 0<= Running_Time <= (0+Ro0) :
            Rate = r1
            X = 1
        elif Ro <= Running_Time <= (Ro+Ro0) :
            Rate = r2
            X = 1
        elif 2*Ro <= Running_Time <= (2*Ro+Ro0) : 
            Rate = r3
            X = 1
        elif (A_Per*Ro) <= Running_Time <= (A_Per*Ro+Ro0) : 
            Rate = r3
            X = 1
        else: 
            X = 0
        dO = K0*Time_Interval #Demand per time interval
        O = O+dO # Current oxygen demand
        decay = ((-SigC*Time_Interval)*C) - ((BC*Time_Interval)*O*C)  #Decay due to demand of chlorine
        Increase = (Rate*X*Time_Interval) #increase due to dosing period. 
        dC = decay + Increase #Total chanfe in Free Chlorine
        C = C+dC #FRee Chlorine after set time.
        if C < 0:
            C = 0 
        Pre_runningT = i #Running Time.
        List_C.append(C)
    return C

gb3 = gb2

def F_Partitioning_ProcLines(gb3 , NPartitions):
    for j in gb3:
        j = F_Partitioning_W(DF= j,NPartitions= NPartitions) 
        List_GB3.append(j)
    return List_GB3





def Washing_ProcLines (List_GB):
    List_GB3=[]
    for j in gb3:
        j = F_Partitioning_W(DF= j,NPartitions= 40) 
        List_GB3.append(j)
        
    
    for j in List_GB3:
        WashT = len(j.index)
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        
        Blw = 0.47 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
        alpha = 0.52#Inactivation rate of pathogen via FC L/mgmin
        V = (3200 *1000) #L #From Luo et al 2012. 
        Rate = 45.45 #kg/min #From Luo et al 2012. 
        Wash_Time = 2.3 #min 
        c1 = 1/Wash_Time #Reciprocal of average time. 
        L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
        Xl = 0
        Xw =0  #pathogen in process water MPN/ml
        
        L_Xw = []
        L_Xl = []
        for i in Times_W:
            #Defining Initial Contamination
            Time = i
            AvCont = j.at[i,"CFU"] /(j.at[i,"Weight"]*454)
            print(AvCont)
            AvContAfter = AvCont*10**-0.8
            print(AvContAfter)
            C= F_Chloride_lvl(Time_Wash= Time)
            print(C)
            Bws = ((AvCont- AvContAfter)*Rate)/V
            print(Bws)
            CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
            print(CXw)
            Xw = Xw+CXw
            if Xw<0:
                Xw = 0
            print(CXw)
            L_Xw.append(Xw)
            
            Xl = AvCont
            print(Xl)
            CXl = (Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
            print(CXl)
            Xl =Xl +CXl
            if Xl < 0:
                Xl = 0
            print(Xl)
            L_Xl.append(Xl)
            AvCont = Xl
            CFU_2 = AvCont*((j.at[i,"Weight"]*454))
            j.at[i,"CFU"] =  CFU_2
            
    return (List_GB3)
    
    

    









    WashT_Per_Pallet = 41
    Times_W = np.arange(0, (WashT_Per_Pallet+1)*No_Pallets, 1).tolist()
    Times_W = [round(num, 1) for num in Times_W]
    
    Blw = 0.47 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
    alpha = 0.52#Inactivation rate of pathogen via FC L/mgmin
    V = (3200 *1000) #L #From Luo et al 2012. 
    Rate = 45.45 #kg/min #From Luo et al 2012. 
    Wash_Time = 2.3 #min 
    c1 = 1/Wash_Time #Reciprocal of average time. 
    L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
    Xl = 0
    Xw =0  #pathogen in process water MPN/ml
    
    L_Xw = []
    L_Xl = []
    for i in Times_W:
        #Defining Initial Contamination
        Time = i
        AvCont = Trials_2.at[i,"CFU"] /(Trials_2.at[i,"Weight"]*454)
        print(AvCont)
        AvContAfter = AvCont*10**-0.8
        print(AvContAfter)
        C= F_Chloride_lvl(Time_Wash= Time)
        print(C)
        Bws = ((AvCont- AvContAfter)*Rate)/V
        print(Bws)
        CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
        print(CXw)
        Xw = Xw+CXw
        if Xw<0:
            Xw = 0
        print(CXw)
        L_Xw.append(Xw)
        
        Xl = AvCont
        print(Xl)
        CXl = (Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
        print(CXl)
        Xl =Xl +CXl
        if Xl < 0:
            Xl = 0
        print(Xl)
        L_Xl.append(Xl)
        AvCont = Xl
        CFU_2 = AvCont*((Trials_2.at[i,"Weight"]*454))
        Trials_2.at[i,"CFU"] =  CFU_2 
    
    
sum(TrialDF.CFU)  
sum(Trials_2.CFU)
plt.plot(L_Xl)


def F_Partitioning_W(DF,NPartitions):
    AllParts_Cont = []
    for i, row in DF.iterrows():
        Cont = DF.at[i,'CFU']
        PartCont=np.random.multinomial(Cont,[1/NPartitions]*NPartitions,size=1)
        PartCont = PartCont[0]
        AllParts_Cont.append(PartCont)
    b_flat = [j for i in AllParts_Cont for j in i]
    newdf = pd.concat([DF]*NPartitions,axis=0)
    newdf=newdf.sort_values(by=['PalletNo'])
    #Pallet_List=(list(range(1,NPartitions+1)))
    newdf["Part"] =list(range(1,len(newdf.index)+1))#np.tile(Pallet_List, len(newdf)//NPartitions)
    newdf = newdf.reset_index(drop=True)  
    newdf.Weight=newdf.Weight/NPartitions
    newdf.CFU = b_flat
    return newdf





    CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
    CXl = (Blw*Xw) -( alpha*Xl*C) - (C1*Xl)
















C = F_Chloride_lvl(Time_Wash= 120)
List_3 = C[1]
plt.plot(List_3)

print(C)


TrialDF["CFU"] = 0
TrialDF.at[0,"CFU"] = 200000




WaterCont_List = []


TrialDF = gb2[0]
Agg_WashTime =0
WashT_Per_Pallet = 40 #min
Agg_WashTime = Agg_WashTime+WashT_Per_Pallet


V = (3200 *1000) #L #From Luo et al 2012. 
Rate = 45 #kg/min #From Luo et al 2012. 
Wash_Time = 2.3 #min 
c1 = 1/Wash_Time #Reciprocal of average time. 
L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
Xl = 0
Xw =0  #pathogen in process water MPN/ml

Weight_Pallet = TrialDF["Weight"]*2.2

for i, row in TrialDF.iterrows():
    AvCont = TrialDF.at[i,"CFU"] /(TrialDF.at[i,"Weight"]*454)
    AvContAfter = AvCont*10**-0.8
    WashT_Per_Pallet = (TrialDF.at[i,"Weight"]*2.2)/45
    round(WashT_Per_Pallet,1)
    print(WashT_Per_Pallet)
    Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
    C= F_Chloride_lvl(Time_Wash= 120)
    alpha = 0.75
    Bws = ((AvCont- AvContAfter)*Rate)/V
    C1 = 2.3/ WashT_Per_Pallet
    
    CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
    Xw = Xw+CXw
    
    Xl = AvCont
    CXl = (Blw*Xw) -( alpha*Xl*C) - (C1*Xl)
    Xl =Xl +CXl
    print(Xl)
    AvCont = AvCont+CXl
    CFU_2 = AvCont*((TrialDF.at[i,"Weight"]*454))
    TrialDF.at[i,"CFU"] =  CFU_2

print(WaterCont_List)
WaterCont_List




#%%

for i in Times: 
    Running_Time = i
    Time_Interval = Running_Time-(Pre_runningT)
    if 0<= Running_Time <= (0+Ro0) :
        Rate = r1
        X = 1
        print("Yes")
    elif Ro <= Running_Time <= (Ro+Ro0) :
        Rate = r2
        X = 1
    elif 2*Ro <= Running_Time <= (2*Ro+Ro0) : 
        Rate = r3
        X = 1
    else: 
        X = 0
    dO = K0*Time_Interval
    O = O+dO

    decay = ((-SigC*Time_Interval)*C) - ((BC*Time_Interval)*O*C) 
    Increase = (Rate*X*Time_Interval)
    dC = decay + Increase
    L_dC.append(dC)
    C = C+dC
    L_c.append(C)
    Pre_runningT = i




































