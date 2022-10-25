# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 09:31:30 2021

@author: Gustavo Reyes
"""
                                                       #Inputs for the Model, Initial Inputs
import ScenCondz
import ContCondz
import numpy as np  
import SCInputz 
import math
from numpy.random import Generator, PCG64
import pandas as pd
import random
rng = Generator(PCG64())
       
#Randomized Initial Contamination

#Contamination Scenario: 
    
#1 Pre-Harvest Irrigation
#2. Point Source Contamination 0.3% of the field
#3. Sporadic from soil manure
#4. Harvest Contmaination #Pending
#5. Custom Plan #Pending


#Contamination Scenario #1
if ScenCondz.Contamination_Scenario == 1: 
    Hazard_Lvl = 100_000
    Cont_Cluster = 1
    Cluster_Size = SCInputz.Field_Weight
    if ScenCondz.Holding_Time == True: #Should always be true unless scenario analysis. 
        Time_CE_H = np.random.triangular(2,4,8) #days
    elif ScenCondz.Holding_Time == False:
        Time_CE_H = np.random.triangular(0,4,8) #days
      
#Contamination Scenario #1
if ScenCondz.Contamination_Scenario == 2:
    Hazard_Lvl = 100_000
    Cont_Cluster = 1
    Cluster_Size  =1_000 
    if ScenCondz.Holding_Time == True: #Should always be true unless scenario analysis. 
        Time_CE_H = np.random.triangular(2,4,8) #days
    elif ScenCondz.Holding_Time == False:
        Time_CE_H = np.random.triangular(0,4,8) #days 
#Contamination Scenario #1
if ScenCondz.Contamination_Scenario == 3:
    Hazard_Lvl = 100_000
    Cont_Cluster = 1
    Cluster_Size = 10_000
    if ScenCondz.Holding_Time == True: #Should always be true unless scenario analysis. 
        Time_CE_H = np.random.triangular(2,4,8) #days
    elif ScenCondz.Holding_Time == False:
        Time_CE_H = np.random.triangular(0,4,8) #days 

#%% Functions 
def pert(a, b, c, *, size=1, lamb=4):
    r = c - a
    alpha = 1 + lamb * (b - a) / r
    beta = 1 + lamb * (c - b) / r
    return a + rng.beta(alpha, beta, size=size) * r

def F_Chloride_lvl (Time_Wash, Treatment):
    #Function Inputs. 
    #Changing times to 0.1 increments.
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
    #Addition Rates
    if Treatment ==1:
        r1= 12.75 #(mg/(ml/min**2))
        r2 = 7.47 #(mg/(ml/min**2))
        r3 = 5.56 #(mg/(ml/min**2))
    elif Treatment == 0:
        r1= 0 #(mg/(ml/min**2))
        r2 =0 #(mg/(ml/min**2))
        r3 =0 #(mg/(ml/min**2))
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
        #if(i==10):
        #   print(O)
        List_C.append(C)
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': List_C,
    })
    return Cdf

def F_Chloride_lvl_Constant(Time_Wash, C_level):
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': C_level,
    })
    return Cdf

#%% Contamination Scenario Inputs
'''
Scenario_no = np.random.choice([1,2,3,4])
    #1. Refers to irrigation water contaminating field
    #2. Refers to splash from irrigation water whenimproper manure is applied
    #3. Refers to splash due to ferral swine feces, clustered
    #4. Referes to splash due to cattle and ferral swine runoff due to irrigation.
Month_choice = np.random.choice([1,2,3,4,5,6,7,8,9,10,11,12 ])
'''      

#%% Inputs for Initial Contamination

#Pre-Contamination Information
'''
Total_Days = 14 #Refers to the total days that the model simulates field contamination and irrigation


#Time Options for the Sampling Scenarios 1,2,3
if ScenCondz.PHS_Int ==True :
    Time_PHS_H = 0 #Days #time from pre harvest sampling to harvest Sampling
elif ScenCondz.PHS_4h ==True:
    Time_PHS_H = 0.166 # 4Hours before harvest
elif ScenCondz.PHS_4d ==True:
    Time_PHS_H = 4 #4 days before harvest
else:
    Time_PHS_H = 0 #if not, baseline is always 4 hours

PHS_Day = math.ceil(Time_PHS_H)

#Times from Initial to PHS:
Time_I_PHS = Total_Days-Time_PHS_H
Time_I_PHS_FullD=int(math.modf(Time_I_PHS)[1])
Time_I_PHS_PartD=math.modf(Time_I_PHS)[0]

#Times from PHS to H
Time_PHS_H_FullD =int(math.modf(Time_PHS_H)[1])
Time_PHS_H_PartD =math.modf(Time_PHS_H)[0]

#Holding Time Days:

if ScenCondz.Holding_Time == True: #Should always be true unless scenario analysis. 
    Holding_Time = int(pert(2,4,8)) #Irrigation
elif ScenCondz.Holding_Time == False:
    Holding_Time = 0 #Irrigation
'''
#%% Irrigation information
'''
irrigation_frequency = int(np.random.uniform(5,7))
probability_irrigation_day = (1/irrigation_frequency) #irrigation every 7 days. 
irrigation_days=list(range(1,(14-Holding_Time+1))) #1 is 1st day, 14 is harvest.
Normalized_Irrigation=len(irrigation_days)/7
Final_Irrigation_Days = [i for i in irrigation_days if (probability_irrigation_day*Normalized_Irrigation)>np.random.uniform(0,1)]
'''

                                                                                             
#%% Die off

#this no currently used. 
# 0 Die-off in Field for first 3 contamination Events
#Break_Point=Funz.Func_NormalTrunc(0.11,3.71,0.68,0.98) #Breaking point for Die Off from Belias et al. 
#Dieoff1 = Funz.F_DieOff1() #from Belias et al. 
#Dieoff2 = Funz.F_DieOff2() #from Belias et al. 

#Lag Consumed for Growth models, always starts at 0, updated in model.
Lag_Consumed_Prev = 0

#%% Time Pre-Harvest HArvest

#Time between Contamination Event and Harvest Sampling *** Scenario Control

#Time Options for the Sampling Scenarios 1,2,3
if ScenCondz.PHS_Int ==True :
    Time_PHS_H = 0 #Days #time from pre harvest sampling to harvest Sampling
elif ScenCondz.PHS_4h ==True:
    Time_PHS_H = 0.166 # 4Hours before harvest
elif ScenCondz.PHS_4d ==True:
    Time_PHS_H = 4 #4 days before harvest
else:
    Time_PHS_H = 0 #if not, baseline is always 4 hours

    
#Here because of math
if Time_PHS_H>Time_CE_H:
    Time_CE_PHS = 0 #No time in contamination Event. 
elif Time_PHS_H<=Time_CE_H:
    Time_CE_PHS= Time_CE_H-Time_PHS_H #Days Time from Contamination Event (Irrigation) to Pre-Harvest Sampling

#print(Time_CE_H)
#%% Pre-Cooling- HArvest

#Chrlorine Pray
if SCInputz.Spray_WashYN == True:
    Harvest_Cspray_red = np.random.uniform(1.1,1.46) # log reduction achieved by chlorination spray. Prachi Pahariya 2022. #revisisit. 
if SCInputz.Spray_WashYN == False:
    Harvest_Cspray_red = 0 # log reduction achieved by chlorination spray. Prachi Pahariya 2022. #revisisit. 

#PreCooling. 

Time_H_PreCooling = np.random.uniform(2,4) #time between harevest and precooling
Temperature_H_PreCooling = np.random.uniform(15,17) #koselki et al paper

if SCInputz.Pre_CoolingYN==True:
    Time_PreCooling = np.random.uniform(0.4,0.6)  #Hours 25-35 minutes
    Temperature_PreCooling = np.random.uniform(3,4) #Precooling 3 C https://www.buschvacuum.com/us/en/news-media/pre-cooling-lettuces-reliably-thanks-to-cutting-edge-vacuum-technology.html
if SCInputz.Pre_CoolingYN==False:
    Time_PreCooling = 0  #Hours 25-35 minutes
    Temperature_PreCooling = np.random.uniform(3,4) #Precooling 3 C https://www.buschvacuum.com/us/en/news-media/pre-cooling-lettuces-reliably-thanks-to-cutting-edge-vacuum-technology.html

#%% Field Pack
#Field Packed Lettuce:
Case_Weight_FieldPack = 25 #lb


#%% Receiving

#Storage at Receiving
Time_Storage_R = np.random.uniform(0,1)#np.random.uniform(5,7)#need reference assumed
Temperature_Storage_R = np.random.triangular(0,8,14) #refrigeration temperature at storage
#PAllet size
Pallet_Weight = 4000  #weight of pallet in lb. 
 

#%% Processing
#Time_R_PC = 2  #Days around 2hr, Time between 

Processing_Lines =4 #np.random.choice([1,2,3,4]) #Number of processing lines of product

# B Cold Storage
#Infor for cold storage die-off. 
Time_ColdStorage = 24 #Houts
Temperature_ColdStorage = np.random.normal(5.1,0.27) #C 

#Food Processing Events, Transfer Rates from 
#TR1 = np.random.triangular(0,0.01,0.02) #Transfer (%) from contaminated lettuce to flume
#TR2 = np.random.triangular(0,0.02,0.02) #Transfer (%) from contaminated lettuce to shredder
#TR3 = np.random.triangular(0,0.01,0.02) #Transfer (%) from contaminated lettuce to shaker
#TR4 = np.random.triangular(0.01,0.04,0.08) #Transfer (%) from contaminated lettuce to centrifuge
#TR5 = np.random.triangular(0,0.01,0.024) #Transfer (%) from contaminated lettuce to conveyor
#Tr_P_S = TR1+TR2+TR3+TR4+TR5 #Total transfer from Contminated Lettuce to Surfaces. 
#Tr_S_P = np.random.triangular(9.9,15.33,18.83) #Total Transfer from Surfaces to new product. 
#ContS=0 #Initial Contamination in Surfaces. 

#Processing Transfer Rates in %/100:
    #Shredder
Tr_Sh_P =np.random.triangular(0.16,0.20,0.28)
Tr_P_Sh =np.random.triangular(0,0.0025,0.0053)

if SCInputz.Sanitation_YN == True:
    #Sanitation Shredder.
    Sh_Compliance = np.random.choice([0,0.25,0.50,0.75,1])
    Sh_San_freq = np.random.choice([2500,5000,7500])
    Sh_San_Eff = int(np.random.choice([-1,-2,-3,-4]))

if SCInputz.Sanitation_YN == False:
    #Sanitation Shredder.
    Sh_Compliance = 0
    Sh_San_freq = 0
    Sh_San_Eff = 0
    

    #Conveyor Belt
Tr_Cv_P =np.random.triangular(0.15,0.18,0.22)
Tr_P_Cv =np.random.triangular(0,0.0062,0.0139)

if SCInputz.Sanitation_YN == True:
    #Sanitation CV
    Cv_Compliance = np.random.choice([0,0.25,0.50,0.75,1])
    Cv_San_freq = np.random.choice([2500,5000,7500])
    Cv_San_Eff = int(np.random.choice([-1,-2,-3,-4]))

if SCInputz.Sanitation_YN == False:
    Cv_Compliance = 0
    Cv_San_freq = 0
    Cv_San_Eff = 0
   


#Flume tank washing step
Wash_Rate = 100 #lb/min

Optimized_washing_clevel = 10 #ppm

if SCInputz.Washing_YN == True:
    if SCInputz.Always_Washing_Optimized == False:
        SCInputz.Washing_Optimized=np.random.choice([True,False])
    elif SCInputz.Always_Washing_Optimized == True:
        SCInputz.Washing_Optimized=True
        
    if SCInputz.Washing_Optimized== True:
        DF_Chlevels = F_Chloride_lvl_Constant(Time_Wash = 300, C_level = Optimized_washing_clevel)
    else:
        DF_Chlevels = F_Chloride_lvl(300, Treatment =1) #Simlating Chlorine levels after time.
    
    
if SCInputz.Washing_YN == False:
    DF_Chlevels = F_Chloride_lvl(300, Treatment =0) #Simlating Chlorine levels after time. without any treament. 

    #Shaker Table
Tr_St_P =np.random.triangular(0.06,0.28,0.30)
Tr_P_St =np.random.triangular(0,0.0006,0.0038)

if SCInputz.Sanitation_YN == True:
        #Sanitation Shaker Table
    St_Compliance = np.random.choice([0,0.25,0.50,0.75,1])
    St_San_freq = np.random.choice([2500,5000,7500])
    St_San_Eff = int(np.random.choice([-1,-2,-3,-4]))
    
if SCInputz.Sanitation_YN == False:
    St_Compliance = 0
    St_San_freq = 0
    St_San_Eff = 0


    #Centrifuge
Tr_C_P =np.random.triangular(0.23,0.27,0.31)
Tr_P_C =np.random.triangular(0,0.0035,0.0159)

if SCInputz.Sanitation_YN == True:
        #Centrifuge compliance.
    C_Compliance = np.random.choice([0,0.25,0.50,0.75,1])
    C_San_freq = np.random.choice([2500,5000,7500])
    C_San_Eff = int(np.random.choice([-1,-2,-3,-4]))
    
if SCInputz.Sanitation_YN == False:
    #Centrifuge compliance.
    C_Compliance =0
    C_San_freq = 0
    C_San_Eff = 0


  
#%% Final Product

#Lb of product per hour
FP_Hour = 7_200
Pack_Weight_FP = 5 #Weight of each pack. 
N_Lots_FP = 1 #Lost of final product
#Packaging into Cases
Case_Weight = 20 #lb #Only if not field packed


#%% Post Processing Storage
Time_PostPStorage = np.random.uniform(0,1) #hr #no post processing storage https://ohioline.osu.edu/factsheet/aex-261
Temperature_ColdStorage = np.random.uniform(4,6) #C 

#%%Transportartion post processing
Trasnportation_Time = np.random.uniform(31.17,83.92) #h from Ryser et al brown all ranges. 
Transportation_Temp = np.random.triangular(0,4,6) #C

#%%Customer Storage
Time_PostPCS = np.random.uniform(24,72) #Zheng et al temperature   #hr
Temperature_PostPCS = np.random.triangular(0,8,14) #C #WZheng et al temperature , also from  mandamba model



    


