# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 13:59:11 2022

@author: Gustavo Reyes
"""

#Initial Contamination figuring it out
import numpy as np
import SCInputz
import random
import Funz
from scipy.stats import beta
import InFunz
import ContScen
import Inputz
import math


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

#Daily Event

Holding_Time = int(np.random.triangular(2,4,8))
probability_irrigation_day = (1/7)
days=list(range(1,(14-Holding_Time)+1))
random.choices(days,[0.8]*len(days))

days_picked = [i for i in days if probability_irrigation_day>np.random.uniform(0,1)]

F_InitialCont()





#Scenario 1, irrigation water contaminated.
def Irrigation_Water_Cont():
    Ecoli_soil = np.random.normal(0.549,0.816) #log10 CFU/g
    P_Ecoli_Water = 0.35
    PosYN_Ecoli_Water = np.random.binomial(1,P_Ecoli_Water)
    C_Ecoli_Water = (10**np.random.normal(0.604,0.357))/10 #CFU/ml water
    
    Tr_Irr_Water = np.random.uniform(1.8,21.6) #ml/g produce
    Pr_Irr_Splashing = pert(0.02,0.04,0.05)
    Pr_Rain_Splashing = 1
    Soil_Transfer = betagen(0.4,0.8,0.05,16.4)
    P_Soil_Plant = np.random.uniform(0.35,0.9)
    
    Increase_Irrigation =C_Ecoli_Water*Tr_Irr_Water*PosYN_Ecoli_Water #CFU/g
    Increase_IrrigationCFU = Increase_Irrigation*(454*100_000)
    return Increase_IrrigationCFU
        
#Scenario 2 Rain splashing
Increase_RainSplash = (10**Ecoli_soil)*Soil_Transfer*P_Soil_Plant*Pr_Rain_Splashing
Increase_RainSplashCFU = Increase_RainSplash*(454*100_000)

#Scenario 3 Irrigation water splash
Increase_IrrigationSplash = (10**Ecoli_soil)*Soil_Transfer*P_Soil_Plant*Pr_Irr_Splashing
Increase_IrrigationSplashCFU = Increase_IrrigationSplash *(454*100_000)

#Scenario 4 Ferral Swine Defecation foliar irrigation: (Clustered Contamination)
amount_ferralloc = 4260 #mass of feces a day #456 g
amount_feralprec = 1_000 #g
amount_cattleprec = 10_000 #g

Cont_Ferral=np.random.lognormal(-1.26,2.2)
Cont_Ferral_CFU=10**Cont_Ferral
CFU_ferralloc_CFU = amount_ferralloc*Cont_Ferral_CFU
print(CFU_ferralloc_CFU)

def transfer_1(val):
    if val<=0.619:
        trans = 1.3*(val-0)/(0.619-0)
    elif (val>=0.619) and (val<=0.946):
        trans = 1.3+((340-1.3)*(val-0.619)/(0.946-0.619))
    else:
       trans = 340+((230000-340)*(val-0.946)/(1-0.946)) 
    return trans

#Transfer Ratio
tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
print(tr)

CFU_ferralloc_CFU*tr

#Scenario 5 Ferral swine defectation runoff
Month = np.random.choice([1,2,3,4,5,6,7,8,9,10,11,12 ])
Month_Index = Month-1
Irrigation_m = [0.861194671,0.921357972,0.972496777,0.484588038,0.135367426,0.03824667,0,0.012462398,0.03008165,0.280189085,0.496347228,0]
Irrigation_m_choice = Irrigation_m[Month_Index]
Cont_ferral_means =[-29.13,-29.13,-29.13,-29.13,-29.13,-29.13,-1.26,-2.1,-29.13,-2.75,-2.35,-2.13] 
Cont_ferral_sds = [9.72,9.72,9.72,9.72,9.72,9.72,2.2,2.44,9.72,2.62,2.51,2.46]
Cont_cattle_means = [-10.62,-7.07,-25.99,3.42,-6.15,-36.01,-1.15,-3.19,-3.62,-1.68,-6.75,-8.68]
Cont_cattle_sds= [5.12,4.16,9.25,3.18,3.91,11.95,2.57,3.11,3.23,2.71,4.07,4.56]

Cont_ferral = np.random.lognormal(Cont_ferral_means[Month_Index],Cont_ferral_sds[Month_Index])
Cont_ferral_CFU = 10**Cont_ferral
Cont_cattle = np.random.lognormal(Cont_cattle_means[Month_Index],Cont_cattle_sds[Month_Index])
Cont_cattle_CFU = 10**Cont_cattle

amount_cattleprec_Final=amount_cattleprec*Irrigation_m_choice
amount_ferralprec_Final = amount_feralprec*Irrigation_m_choice

CFU_ferralprec = amount_ferralprec_Final*Cont_ferral_CFU
CFU_cattleprec = amount_cattleprec_Final*Cont_cattle_CFU


tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
Total_CFU=CFU_ferralprec*tr
print(Total_CFU)



#Creating th DF: 
    
df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                  Field_Weight = SCInputz.Field_Weight, 
                  slot_number = SCInputz.slot_number)






#Beginning. 


Holding_Time = int(pert(0,2,8)) #Irrigation

Total_Days = 14
PHS_Time = 0.1666
#Times from Initial to PHS:
Time_I_PHS = Total_Days-PHS_Time
Time_I_PHS_FullD=int(math.modf(Time_I_PHS)[1])
Time_I_PHS_PartD=math.modf(Time_I_PHS)[0]

#Times from PHS to H
Time_PHS_H =Total_Days-Time_I_PHS
Time_PHS_H_FullD =int(math.modf(Time_PHS_H)[1])
Time_PHS_H_PartD =math.modf(Time_PHS_H)[0]

#Holding Time Days
Time_Cutoff_Irrigation = Total_Days-Holding_Time


#function inputs
scenario_no = 1

TotalTime= 0

#PreHarvest Section 1: From Beginning - PHS

#Contamination part 1: From Beginning-to-PHS
for i in range(Time_I_PHS_FullD): #Loop to repeat over the dats
    #Determining Contamination Scenario
    if scenario_no ==1:
        #Calculating total CFU in Field
        Cont = Irrigation_Water_Cont()
        if TotalTime<Time_Cutoff_Irrigation:
            #Adding Uniform Contmaination to the field
            df = ContScen.F_systematic_C(df =df, 
                          Hazard_lvl = Cont,
                          No_Cont_Clusters=1,
                          Cluster_Size=100_000, 
                          Partition_Weight= 50)
            
        df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(1))
    TotalTime = TotalTime+1
    
    
if Time_I_PHS_PartD >0: #Only if the partial day is greater than o
    #Part 2: Partial Day, deciding if contamination occurs in the first or second part of the day
    if np.random.uniform(0,1)<Time_I_PHS_PartD: 
        Happens_YN_1part=1 
    else: 
        Happens_YN_1part=0
    
    if Happens_YN_1part==1:
            Cont = Irrigation_Water_Cont()   
            if TotalTime<Time_Cutoff_Irrigation:
                #Adding Uniform Contmaination to the field
                df = ContScen.F_systematic_C(df =df, 
                              Hazard_lvl = Cont,
                              No_Cont_Clusters=1,
                              Cluster_Size=100_000, 
                              Partition_Weight= 50)   
    
    #Dieoff for first part of the day:
    df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(Time_I_PHS_PartD))
    TotalTime = TotalTime+Time_I_PHS_PartD
       
#Sampling Here

#Sampling Function in Here

#Part 3 Other days from PHS to 

if Time_PHS_H_PartD > 0: #If Time between PHS-H is less than 1
    if Happens_YN_1part==0:
        Cont = Irrigation_Water_Cont()
            
        if TotalTime<Time_Cutoff_Irrigation:
            #Adding Uniform Contmaination to the field
            df = ContScen.F_systematic_C(df =df, 
                          Hazard_lvl = Cont,
                          No_Cont_Clusters=1,
                          Cluster_Size=100_000, 
                          Partition_Weight= 50) 
        
        df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(Time_PHS_H_PartD))
        TotalTime = TotalTime+Time_PHS_H_PartD


for i in range(Time_PHS_H_FullD):
    if scenario_no ==1:
        #Calculating total CFU in Field
        Cont = Irrigation_Water_Cont()
        
        #Adding Uniform Contmaination to the field
        if TotalTime<Time_Cutoff_Irrigation:
            #Adding Uniform Contmaination to the field
            df = ContScen.F_systematic_C(df =df, 
                          Hazard_lvl = Cont,
                          No_Cont_Clusters=1,
                          Cluster_Size=100_000, 
                          Partition_Weight= 50) 
        
    df = Funz.Applying_dieoff(df=df, Dieoff=Funz.F_Simple_DieOff(1))
    TotalTime = TotalTime+1
    

    



    
     


















#We will do the final 14 DAYS
Contamination_Event= np.random.choice([1,2,3,4])
#1 = Contaminated Irrigation water
#2=  Contaminated feces from feral swine
#3 = Contamination from runoff of contamination 

Month = np.random.choice([1,2,3,4,5,6,7,8,9,10,11,12 ])
Irrigation_m = [0.861194671,0.921357972,0.972496777,0.484588038,0.135367426,0.03824667,0,0.012462398,0.03008165,0.280189085,0.496347228,0]


#cattle feces
np.random.normal(-1.15,2.57) # log CFU _gram
    # 2 ways to soil
    #1. Direct defectation
    #2. Precipitation runoff
    
#Precipitation: 
PrecNorm = (48-0.1)/(59.1-0.1) #normalization of precipitation for a month

amount_ferralloc = 4260 #mass of feces a day #456 g
amount_feralprec = 1_000 #g
amount_cattleprec = 10_000 #g


np.random.lognormal(-29.13,9.72,1)*454*100_000
