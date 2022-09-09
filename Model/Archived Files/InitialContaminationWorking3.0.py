# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 08:47:26 2022

@author: gareyes3
"""

# Libraries
import numpy as np
import SCInputz
import random
import Funz
from scipy.stats import beta
import InFunz
import ContScen
import Inputz
import math
import Inputz
from importlib import reload  
import matplotlib.pyplot as plt
import pandas as pd
from numpy.random import Generator, PCG64
rng = Generator(PCG64())

#Assisting Functions
def pert(a, b, c, *, size=1, lamb=4):
    r = c - a
    alpha = 1 + lamb * (b - a) / r
    beta = 1 + lamb * (c - b) / r
    return float(a + np.random.beta(alpha, beta, size=size) * r)

def betagen(a,b,mini,maxi):    
    return float(beta.rvs(a , b)* (maxi - mini) + mini)


#Contamination from irrigation water pang et al
def F_InitialCont():
    #Using base pert distribution
    #Calculation of total CFUs
    Cont_CFU_g = float(pert(0,0,634))
    g_field = SCInputz.Field_Weight*454 #454 g in 1 lb. 
    Final_Cont = int(Cont_CFU_g*g_field)
    return Final_Cont

def transfer_1(val):
    if val<=0.619:
        trans = 1.3*(val-0)/(0.619-0)
    elif (val>=0.619) and (val<=0.946):
        trans = 1.3+((340-1.3)*(val-0.619)/(0.946-0.619))
    else:
       trans = 340+((230000-340)*(val-0.946)/(1-0.946)) 
    return trans


def lognormal_max(mean,sigma):
    while True:
        age = np.random.lognormal(mean=mean,sigma=sigma,size=1)
        if age>=0.84 and age <= 8.4:
            break
    return age

#%%

#Scenario 1: Uniform Contamination
def Irrigation_Water_Cont():
    #Ecoli_soil = np.random.normal(0.549,0.816) #log10 CFU/g
    #P_Ecoli_Water = 0.35
    #PosYN_Ecoli_Water = np.random.binomial(1,P_Ecoli_Water)
    C_Ecoli_Water = (10**np.random.normal(0.604,0.357))/10 #CFU/ml water
    
    Tr_Irr_Water = np.random.uniform(1.8,21.6) #ml/g produce
    #Pr_Irr_Splashing = pert(0.02,0.04,0.05)
    #Pr_Rain_Splashing = 1
    #Soil_Transfer = betagen(0.4,0.8,0.05,16.4)
    #P_Soil_Plant = np.random.uniform(0.35,0.9)
    
    Increase_Irrigation =C_Ecoli_Water*Tr_Irr_Water #CFU/g
    Increase_IrrigationCFU = Increase_Irrigation*(454*100_000)
    return Increase_IrrigationCFU

#Scenario 2: Uniform Contamination
def Irrigation_Soil_Splash():
    Ecoli_soil = np.random.normal(0.549,0.816) #log10 CFU/g
    #Pr_Irr_Splashing = pert(0.02,0.04,0.05)
    Soil_Transfer = betagen(0.4,0.8,0.05,16.4)
    P_Soil_Plant = np.random.uniform(0.35,0.9)
    Cont = (10**Ecoli_soil)*Soil_Transfer*P_Soil_Plant
    Increase_IrrigationCFU = Cont*(454*100_000)
    return int(Increase_IrrigationCFU) 

#Scenario 3
def Feces_Addition(Month): #Wrose case scenrio, High level of contamination.
    Month_Index = Month-1
    Cont_ferral_means =[-29.13,-29.13,-29.13,-29.13,-29.13,-29.13,-1.26,-2.1,-29.13,-2.75,-2.35,-2.13] 
    Cont_ferral_sds = [9.72,9.72,9.72,9.72,9.72,9.72,2.2,2.44,9.72,2.62,2.51,2.46]
    Ferral_Day = 4260 #g of feces per day
    #tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
    Cont_Ferral=lognormal_max(Cont_ferral_means[Month_Index],Cont_ferral_sds[Month_Index])
    Cont_Ferral_CFU=10**Cont_Ferral
    CFU_ferralloc_Soil_CFU = Ferral_Day*Cont_Ferral_CFU #CFU/g
    return float(CFU_ferralloc_Soil_CFU)
    #CFU_ferralloc_CFU*tr #CFU/g trafered to Sublot
    
def Feces_Splash(df,Soil_Slots):
    tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
    Cont_Trans = Soil_Slots.apply(lambda x: rng.binomial(x,tr))
    df["CFU"] = df["CFU"]+Cont_Trans
    return df

#Scenario 4: 
def Feces_Runoff(Month): 
    amount_feralprec = 1_000 #g
    amount_cattleprec = 10_000 #g
    Month_Index = Month-1
    Irrigation_m = [0.861194671,0.921357972,0.972496777,0.484588038,0.135367426,0.03824667,0,0.012462398,0.03008165,0.280189085,0.496347228,0]
    Irrigation_m_choice = Irrigation_m[Month_Index]
    Cont_ferral_means =[-29.13,-29.13,-29.13,-29.13,-29.13,-29.13,-1.26,-2.1,-29.13,-2.75,-2.35,-2.13] 
    Cont_ferral_sds = [9.72,9.72,9.72,9.72,9.72,9.72,2.2,2.44,9.72,2.62,2.51,2.46]
    Cont_cattle_means = [-10.62,-7.07,-25.99,3.42,-6.15,-36.01,-1.15,-3.19,-3.62,-1.68,-6.75,-8.68]
    Cont_cattle_sds= [5.12,4.16,9.25,3.18,3.91,11.95,2.57,3.11,3.23,2.71,4.07,4.56]
    
    Cont_ferral = lognormal_max(Cont_ferral_means[Month_Index],Cont_ferral_sds[Month_Index])
    Cont_ferral_CFU = 10**Cont_ferral
    Cont_cattle =lognormal_max(Cont_cattle_means[Month_Index],Cont_cattle_sds[Month_Index])
    Cont_cattle_CFU = 10**Cont_cattle
    
    amount_cattleprec_Final=amount_cattleprec*Irrigation_m_choice
    amount_ferralprec_Final = amount_feralprec*Irrigation_m_choice
    
    CFU_ferralprec = amount_ferralprec_Final*Cont_ferral_CFU
    CFU_cattleprec = amount_cattleprec_Final*Cont_cattle_CFU
    
    Total_Soil_Add  =  CFU_ferralprec + CFU_cattleprec
    return float(Total_Soil_Add)
    


#%%
df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                  Field_Weight = SCInputz.Field_Weight, 
                  slot_number = SCInputz.slot_number)


Soil_Slots = df["CFU"]
#Contaminating 10 random partitions
Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Addition(4),[1/2000]*2000,1)
Soil_Slots=Updates_values[0]
Soil_Slots=pd.Series(Soil_Slots)
Soil_Slots = Soil_Slots.apply(lambda x: np.random.binomial(x,10**-0.1744)) #applying daily die-off.            
if 3 in Inputz.Final_Irrigation_Days:
    df = Feces_Splash(df,Soil_Slots)




#%% Contamination Function Pre-Harvest. 
def Cont_Ini_PHS(df):
    TotalTime = 0 #Time Tracker
    Total_CFU_v = []
    #Contamination part 1: From Beginning-to-PHS
    for i in range(Inputz.Time_I_PHS_FullD): #Loop to repeat over the full days before PHS
        #Determining Contamination Scenario
        
        #Scenario #1
        if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
            #Calculating total CFU in Field for irrigation fields
            if TotalTime in Inputz.Final_Irrigation_Days:
                Cont = Irrigation_Water_Cont()
                #Adding Uniform Contmaination to the field
                df = ContScen.F_systematic_C(df =df, 
                              Hazard_lvl = Cont,
                              No_Cont_Clusters=1,
                              Cluster_Size=SCInputz.Field_Weight, 
                              Partition_Weight= SCInputz.Partition_Weight)
                

        #Scenario #2
        if Inputz.Scenario_no == 2:
            if TotalTime in Inputz.Final_Irrigation_Days:
                Cont = Irrigation_Soil_Splash()
                #Adding Uniform Contmaination to the field
                df = ContScen.F_systematic_C(df =df, 
                              Hazard_lvl = Cont,
                              No_Cont_Clusters=1,
                              Cluster_Size=SCInputz.Field_Weight, 
                              Partition_Weight= SCInputz.Partition_Weight)
        
        #Scenario 3
        if Inputz.Scenario_no == 3:
            Soil_Slots = df["CFU"]
            #Contaminating 10 random partitions
            Index_Cont=Soil_Slots.sample(n=10).index
            Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
            Soil_Slots[Index_Cont]=Updates_values[0]
            Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
            if TotalTime in Inputz.Final_Irrigation_Days:
                df = Feces_Splash(df,Soil_Slots)
        #Scenario 4
        if Inputz.Scenario_no == 4:
            Soil_Slots = df["CFU"]
            #Contaminating 10 random partitions
            Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
            Soil_Slots=Updates_values[0]
            Soil_Slots=pd.Series(Soil_Slots)
            Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
            if TotalTime in Inputz.Final_Irrigation_Days:
                df = Feces_Splash(df,Soil_Slots)        
            
        TotalCFU = df["CFU"].sum()
        Total_CFU_v.append(TotalCFU)
        
        #Applying dieoff to the pants every day.
        df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(1))
        #Adding time to the time counter.
        TotalTime = TotalTime+1
        print(TotalTime)
    
    
        
    if Inputz.Time_I_PHS_PartD >0: #Only if the partial day is greater than 0
        #Part 2: Partial Day, deciding if contamination occurs in the first or second part of the day
        if rng.uniform(0,1)<Inputz.Time_I_PHS_PartD: 
            Happens_YN_1part=1 
        else: 
            Happens_YN_1part=0
        
        if Happens_YN_1part==1:
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)

    
            if Inputz.Scenario_no == 2: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
            
            if Inputz.Scenario_no == 3:  
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_I_PHS_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
            
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_I_PHS_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)  
                        
    
            
        TotalCFU = df["CFU"].sum()  
        Total_CFU_v.append(TotalCFU)
        #Dieoff for first part of the day:
        df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(Inputz.Time_I_PHS_PartD))
        TotalTime = TotalTime+Inputz.Time_I_PHS_PartD
        print(TotalTime)
        
    if Inputz.Time_I_PHS_PartD ==0:
        Happens_YN_1part = 0
        
        
    return [df, Total_CFU_v, TotalTime, Happens_YN_1part]

#%%

def Cont_PHS_End(df, Happens_YN_1part, TotalTime, Total_CFU_v ):
    #Part 2 After Sampling
    if Inputz.Time_PHS_H_PartD > 0: #If Time between PHS-H is less than 1
        if Happens_YN_1part==0:
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
                    
            if Inputz.Scenario_no == 2: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
                    
            if Inputz.Scenario_no == 3:  
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_PHS_H_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
                
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_PHS_H_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots) 
            
        TotalCFU = df["CFU"].sum()  
        Total_CFU_v.append(TotalCFU)
        df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(Inputz.Time_PHS_H_PartD))
        TotalTime = TotalTime+Inputz.Time_PHS_H_PartD
    
    
    for i in range(Inputz.Time_PHS_H_FullD):
            #Determining Contamination Scenario
            
            #Scenario #1
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                #Calculating total CFU in Field for irrigation fields
                if TotalTime in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
    
            #Scenario #2
            if Inputz.Scenario_no == 2:
                if TotalTime in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
            
            #Scenario 3
            if Inputz.Scenario_no == 3:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
                if TotalTime in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
            #Scenario 4
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
                if TotalTime in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)        
                
            TotalCFU = df["CFU"].sum()
            Total_CFU_v.append(TotalCFU)
            #Applying dieoff to the pants every day.
            df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(1))
            #Adding time to the time counter.
            TotalTime = TotalTime+1
    return [df, Total_CFU_v, TotalTime, Happens_YN_1part]

#%%
reload(Inputz)
Inputz.Scenario_no =3
df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                  Field_Weight = SCInputz.Field_Weight, 
                  slot_number = SCInputz.slot_number)
df["CFU"] ==1000

Outs = Cont_Ini_PHS(df)
Outs[1]
df=Outs[0]
plt.plot(Outs[1] ,color ="green")
Inputz.Final_Irrigation_Days
