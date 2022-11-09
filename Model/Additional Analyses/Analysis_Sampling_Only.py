# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 13:53:20 2021

@author: gareyes3
"""

#Analysis of a single sampling plan
#%%
import sys, os
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

# %%
from importlib import reload
import numpy as np
import Listz
import pandas as pd
import MainModel3z
import SCInputz
import Inputz
import ContCondz
import ScenCondz
import InFunz
import OutFunz
import ContScen
import Funz
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import sys
import Trial_MainLoop_PH
import statistics

#%%

Total_grams= SCInputz.Field_Weight*454


CFU_10000g = Total_grams/10000 #-4 log
CFU_1000g = Total_grams/1000  #-3 log
CFU_100g = Total_grams/100 #-2 log
CFU_10g = Total_grams/10  #-1 log
CFU_g =  Total_grams #0 log
CFU_0_1g = Total_grams*10 #1 log

CFU_0_01g = Total_grams*100 #2log

# %% Changing to Point Source Contamination
# Contamination Challenges
ContCondz.Background_C = False
ContCondz.Point_Source_C = True
ContCondz.Systematic_C = False

# Harvester Contamination
ContCondz.Crew_C = False
ContCondz.Harvester_C = False

# Processing equipment
ContCondz.PE_C = False
ContCondz.PE_Cont_Loc = False,  # 1,2,3,4,5
# 1 = Shredder, #2 = Belt, #3 = Washing, #4 Shaker, #5Centrifuge
ContCondz.Pack_C = False



#%%
#Conditions

Tuning_SampleSize =[60,120,240,480,960] #[60,120,300,600,1200] #list(range(100,1200,200))
Tuning_Grabs =[1,30,60,120,240,480]
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
#Tuning_SampType =["4d"]# ["4d","4h","Int"]
Tuning_Cluster_Size =[100,1000,10000,25000,50000,100000]#,10000,50000]
Tuning_Contamination_levels =[CFU_10000g,CFU_1000g,CFU_100g,CFU_10g,CFU_g,CFU_0_1g]


#%%
Iteration1 = 100
Iteration2 = 100
First_Set=[]
Second_Set=[]
SS_Set = []
NG_set = []
Cluster_Size_Set=[]
Contamination_Lvl_Ser=[]

for j in list(range(1,Iteration2)):
    print(j)
    #Selecting the variable Inputs.
    Sample_Size = np.random.choice(Tuning_SampleSize)
    Number_Grabs = np.random.choice(Tuning_Grabs)
    Cluster_Size = np.random.choice(Tuning_Cluster_Size)
    Contamination_lvl = np.random.choice(Tuning_Cluster_Size)
    for i in list(range(1,Iteration1)):

        #Creation of the Data Frame to Track: 
        df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                          Field_Weight = SCInputz.Field_Weight, 
                          slot_number = SCInputz.slot_number)
        
        df=ContScen.F_systematic_C(df=df, 
                                     Hazard_lvl=Contamination_lvl,
                                     No_Cont_Clusters =1, 
                                     Cluster_Size = Cluster_Size, 
                                     Partition_Weight = SCInputz.Partition_Weight)
        Cont_Before=sum(df.CFU)
        #Sampling        
        df = Funz.F_Sampling_2(df =df,Test_Unit ="Lot", 
                                      NSamp_Unit = 1, 
                                      Samp_Size =Sample_Size, 
                                      Partition_Weight =SCInputz.Partition_Weight, 
                                      NoGrab =Number_Grabs )
        #Rejection
        df= Funz.F_Rejection_Rule3(df =df, Test_Unit = "Lot", limit = 0)
         
        Cont_Acc = sum(df.CFU)
        Cont_Rej = Cont_Before-Cont_Acc
        if Cont_Acc == 0:
            Cont_PerR = 1
        else:
            Cont_PerR = Cont_Rej/(Cont_Acc+Cont_Rej) #Percentage Rejected by Finished product sampling
            
        First_Set.append(Cont_PerR)
    mean1 = statistics.mean(First_Set)    
    Second_Set.append(mean1)
    SS_Set.append(Sample_Size)
    NG_set.append(Number_Grabs)
    Cluster_Size_Set.append(Cluster_Size)
    Contamination_Lvl_Ser.append(Contamination_lvl)
    
Data_Frame_Final  = pd.DataFrame(
    {'Sample_Size':  SS_Set,
     'Number_of_Grabs':  NG_set,
     'Cluster_Size': Cluster_Size_Set,
     "Contamination_Level": Contamination_Lvl_Ser,
     "PerReject": Second_Set
    })



Data_Frame_Final.to_csv('Sampling_Plan_Sens.csv')

