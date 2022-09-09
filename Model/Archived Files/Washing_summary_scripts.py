# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 08:31:43 2022

@author: gareyes3
"""

#Determining the effect of chlorine wash: 
    
#Random contmaination. No_Intervention.     

difference_wash = (Baseline_NI_Wash_1[1]["Bef Washing"] - Baseline_NI_Wash_1[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()


difference_wash = (Baseline_AI_1[1]["Bef Washing"] - Baseline_AI_1[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()


difference_wash = (Baseline_AI_1_5ppm[1]["Bef Washing"] - Baseline_AI_1_5ppm[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()

difference_wash = (Baseline_AI_1_10ppm[1]["Bef Washing"] - Baseline_AI_1_10ppm[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()


difference_wash = (Baseline_AI_1_20ppm[1]["Bef Washing"] - Baseline_AI_1_20ppm[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()

difference_wash = (Baseline_AI_1_0ppm[1]["Bef Washing"] - Baseline_AI_1_0ppm[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()

difference_wash = (Baseline_AI_1_False[1]["Bef Washing"] - Baseline_AI_1_False[1]["Bef Shaker Table"] )
difference_wash=difference_wash[difference_wash!=0] #removing 0 values
np.log10(difference_wash).describe()


Changes_1=np.log10(Baseline_NI_Wash_1[1]["Bef Washing"]/Baseline_NI_Wash_1[1]["Bef Shaker Table"]).describe()

Baseline_NI_Wash_1[1]=Baseline_NI_Wash_1[1][Baseline_NI_Wash_1[1]["Bef Washing"]!=0]
Baseline_NI_Wash_1[1]["Bef Shaker Table"][Baseline_NI_Wash_1[1]["Bef Shaker Table"] == 0] = 0.001

(10**4.9) *(1-0.055)*(2.3/(3.2*10**6))


gb4= gb2

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


def F_Chloride_lvl (Time_Wash):
    #Function Inputs. 
    #Changing times to 0.1 increments.
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
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
    O = 301 # Initial Oxygen demand, as per luos initial oxygen demand
    #Other parameters
    SigC = 1.70*(10**-3) #Natural decay of FC
    BC =5.38*(10**-4) #Depletion rate of FC in water. 
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
        if(i==10):
            print(O)
        List_C.append(C)
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': List_C,
    })
    return Cdf

plt.plot(F_Chloride_lvl(36)["C"])
plt.xlabel('Time in 1/100 of a minute')
plt.ylabel('FC Levels')
    
    
import numpy as np
import pandas as pd


    
def F_Washing_ProcLines_2 (List_GB3, Wash_Rate, Cdf):
    for j in List_GB3:
        WashT = len(j.index)
        #DF_Clvl = F_DF_Clvl(WashT)
        
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        V = (3200*1000) #L #From Luo et al 2012. 
        
        Xl = 0
        Xw =0  #pathogen in process water MPN/ml
 
        L_Xw = []
        L_Xl = []
        for i in Times_W:
            Xs = np.random.triangular(0.003,0.055,0.149)
            alpha = np.random.uniform(0.5,0.75)#Inactivation rate of pathogen via FC L/mgmin
            Blw =  np.random.triangular(0.38,0.75,2.2) #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
            Wash_Time = np.random.choice([0.5,0.75,1.0]) 
            c1 = 1/Wash_Time #Reciprocal of average time.
            Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
            L = (Rate*1000)/(c1) #g of lettuce in the tank at the same time
            
            #Defining Initial Contamination
            Time = i #indicating the minute
            AvCont = j.at[i,"CFU"] /(j.at[i,"Weight"]*454) #Contmaination in partition CFU/g
            Cont_CFU= j.at[i,"CFU"]
            #Xl = j.at[i,"CFU"]
            C =   float(Cdf.loc[Cdf['Time'] == Time, 'C']) #Chlorine in ppm
            Bws =Cont_CFU*(1-Xs)*(Wash_Time/V)
            print(Bws, "bws")
            CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C) #Change in Xw: Concentration of ecoli in water
            print(CXw)
            Xw = Xw+CXw #Updated Concentration of ecoli in water
            print(Xw, "Xw")
            if Xw<0:
                Xw = 0
            L_Xw.append(Xw)
            Xl = j.at[i,"CFU"] /(j.at[i,"Weight"]*454) 
            CXl = (Blw*Xw) - (alpha*Xl*C) - (Xl/Wash_Time) #change in lettuce pieces
            print(CXl, "CXL")
            Xl =Xl +CXl #updated pathogen level in lettuce pieces
            print(Xl, "Xl")
            if Xl < 0:
                Xl = 0
            L_Xl.append(Xl)
            #AvCont = Xl
            #CFU_2 = AvCont*((j.at[i,"Weight"]*454))
            j.at[i,"CFU"] =  Xl*((j.at[i,"Weight"]*454))
    return (List_GB3) 

def F_Washing_ProcLines (List_GB3, Wash_Rate, Cdf):
    for j in List_GB3:
        WashT = len(j.index)
        #DF_Clvl = F_DF_Clvl(WashT)
        
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        
        Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
        alpha = 0.75#Inactivation rate of pathogen via FC L/mgmin
        V = (3200 *1000) #L #From Luo et al 2012. 
        Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
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
            AvContAfter = AvCont*10**-0.8
            C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
            #C= F_Chloride_lvl(Time_Wash= Time)
            Bws = ((AvCont- AvContAfter)*Rate)/V
            print(Bws)
            CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
            Xw = Xw+CXw
            if Xw<0:
                Xw = 0
            L_Xw.append(Xw)
            Xl = AvCont
            CXl = (Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
            Xl =Xl +CXl
            if Xl < 0:
                Xl = 0
            L_Xl.append(Xl)
            AvCont = Xl
            CFU_2 = AvCont*((j.at[i,"Weight"]*454))
            j.at[i,"CFU"] =  CFU_2
    return (List_GB3) 


for i in gb4:
    i["CFU"] = 100
    
for i in gb4:
  print( i["CFU"].sum())

for i in gb4:
    i["CFU"][1] = 10000    

chlorine_levs =  F_Chloride_lvl_Constant(Time_Wash = 300, C_level = 0)
chlorine_levs =  F_Chloride_lvl (Time_Wash=300)

gb5 = F_Washing_ProcLines(List_GB3 = gb4, Wash_Rate = 100, Cdf = chlorine_levs)

for i in gb5:
  print( i["CFU"].sum())

np.log10(9918/28000)

2528/(100*454)



chlorine_levs =  F_Chloride_lvl (Time_Wash=37)

10**4.9
79432*.908

#####Validation Section####
#This functions validate the model in https://doi.org/10.1016/j.fm.2015.05.010

def F_Chloride_lvl_Constant(Time_Wash, C_level):
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': C_level,
    })
    return Cdf

df_conts = pd.DataFrame(
    {'CFU': [0]*3600,
     'Weight': 1,
    })

def F_Washing_ProcLines (df , Wash_Rate, Cdf):
    WashT = 36#len(df.index)
    #DF_Clvl = F_DF_Clvl(WashT)
    
    Times_W = np.arange(0, WashT, 0.01).tolist()
    Times_W = [round(num,2) for num in Times_W]
    
    Timeint = 0.01
    
    Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
    alphablw = 0.75*Timeint#Inactivation rate of pathogen via FC L/mgmin
    alpha = 0.75
    V = (3200 *1000) #L #From Luo et al 2012. 
    Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
    Wash_Time = 0.46 #min 
    c1 = 1/Wash_Time #Reciprocal of average time. 
    L = (Rate*1000)/(c1) #g of lettuce in the tank at the same time
    Xl = 0
    Xw =0  #pathogen in process water MPN/ml
    
    L_Xw = []
    L_Xl = []
    for i in Times_W:
        index =int(i*100)
        #Defining Initial Contamination
        Time = i
        AvCont = df.at[index,"CFU"] /(df.at[index,"Weight"]*454)
        print(AvCont)
        C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
        #C= F_Chloride_lvl(Time_Wash= Time)
        Bws = 1.95*Timeint#((AvCont- AvContAfter)*Rate)/V
        CXWfirst = Bws - (Blw*Xw*(L/V))
        CXw =  CXWfirst - (alphablw*Xw*C)
        Xw = Xw+CXw
        if Xw<0:
            Xw = 0
        L_Xw.append(Xw)
        Xl = AvCont 
        CXl = (Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
        Xl =Xl +CXl
        if Xl < 0:
            Xl = 0
        L_Xl.append(Xl)
        #AvCont = Xl
        #CFU_2 = AvCont*((df.at[index,"Weight"]*454))
        #df.at[index,"CFU"] =  CFU_2
        outs = [df, L_Xl, L_Xw]
    return (outs) 


chlorine_levs =  F_Chloride_lvl (Time_Wash=36)

chlorine_levs =F_Chloride_lvl_Constant(36,10)


outs_val = F_Washing_ProcLines (df =df_conts , Wash_Rate = 100, Cdf =chlorine_levs )

#Plot of Free chlorine levels
plt.plot(chlorine_levs['C'])

#Plot of the Xl
plt.plot(outs_val[1])
plt.xlabel('Time in 1/100 of a minute')
plt.ylabel('XL (MPN/g)')

#Plot of the  Xw
plt.plot(outs_val[2])
plt.xlabel('Time in 1/100 of a minute')
plt.ylabel('XW (MPN/g)')

#Updting Washing -----------------------------------------------------------------

def F_Chloride_lvl_Constant(Time_Wash, C_level):
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': C_level,
    })
    return Cdf



df_conts = pd.DataFrame(
    {'CFU': [1]*300,
     'Weight': 100,
    })

before_CFU = df_conts["CFU"].sum()

def F_Washing_ProcLines (df , Wash_Rate, Cdf):
    WashT = len(df.index) #1 minute intervals
    #DF_Clvl = F_DF_Clvl(WashT)
    
    Times_W = np.arange(0, WashT,1).tolist()
    Times_W = [round(num,1) for num in Times_W]
    
    
    Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
    alphablw = 0.75#Inactivation rate of pathogen via FC L/mgmin
    alpha = 0.75
    V = (3500 *1000) #L #From Luo et al 2012. 
    Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
    Wash_Time = 0.46 #min 
    c1 = 1/Wash_Time #Reciprocal of average time. 
    L = (Rate*1000)/(c1) #g of lettuce in the tank at the same time
    Xl = 0
    Xw =0  #pathogen in process water MPN/ml
    
    L_Xw = []
    L_Xl = []
    for i in Times_W:
        Xs = np.random.triangular(0.003,0.055,0.149)
        index =i
        #Defining Initial Contamination
        Time = i
        AvCont = df.at[i,"CFU"] /(df.at[i,"Weight"]*454)
        #AvCont_CFU = df.at[i,"CFU"]
        #AvContAfter = AvCont*10**-0.8
        C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
        Bws = (((AvCont)-(AvCont*Xs))*(Rate*1000))/V
        
        #Bws = ((AvCont- AvContAfter)*Rate)/V
        #print(Bws)
        CXWfirst = Bws - (Blw*Xw*(L/V))
        CXw =  CXWfirst - (alphablw*Xw*C)
        print(CXw, "CXw")
        Xw = Xw+CXw
        if Xw<0:
            Xw = 0
        L_Xw.append(Xw)
        Xl = (AvCont*Xs)
        #print(Xl)
        #CXL23t = (alpha*Xl*C) - (c1*Xl)
        #print(CXL23t, "CXL23t")
        if C>0.5:
            CXL23t = 0.214*np.log(C)+0.220
        if C<0.5:
            CXL23t = 0
        print(C, "C")
        print(CXL23t, "CXL23")
        Xl = Xl*(10**-CXL23t)
        CXl = (Blw*Xw)  #- (alpha*Xl*C) - (c1*Xl)
        #print(Blw*Xw, "fist section")
        #print(CXl, "CXL")
        #print(Xl, "XL")
        Xl =Xl +CXl
        if Xl < 0:
            Xl = 0
        L_Xl.append(Xl)
        AvCont = Xl
        CFU_2 = AvCont*((df.at[index,"Weight"]*454))
        print(CFU_2," CFU_2")
        if CFU_2<1:
            CFU_2= np.random.binomial(1,CFU_2)
        elif CFU_2>1:
            partial =math.modf(CFU_2)
            part1= np.random.binomial(1,partial[0])
            part2= partial[1]
            CFU_2 = part1+part2
        df.at[index,"CFU"] = CFU_2
        print(df.at[index,"CFU"])
        outs = [df, L_Xl, L_Xw]
    return (outs) 


#chlorine_levs =  F_Chloride_lvl (Time_Wash=300)

chlorine_levs =F_Chloride_lvl_Constant(Time_Wash = 300, C_level = 10)

outs_val = F_Washing_ProcLines (df =df_conts , Wash_Rate = 100, Cdf =chlorine_levs )

after_CFU = outs_val[0]["CFU"].sum()

np.log10(before_CFU/after_CFU)

#Plot of Free chlorine levels
plt.plot(chlorine_levs['C'])

#Plot of the Xl
plt.plot(outs_val[1])

#Plot of the  Xw
plt.plot(outs_val[2])

19/3500000
90*128
90*0.128


0.46/3500000
