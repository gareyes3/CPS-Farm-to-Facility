# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 11:13:00 2022

@author: gareyes3
"""

#Marginal Gain analysis sampling, creation of the funtion. 

#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\reyes\Documents\GitHub\CPS-Farm-to-Facility\Model')

# %%
import random
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
import csv
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import sys
#import Trial_MainLoop_PH
import scipy.stats as st


reload(MainModel3z)
reload(Inputz)
reload(SCInputz)
reload(ScenCondz)

#%% MAIN FUCTION THE RUNS THE MODEL WITH THE INTERVENTIONS LISTES
#Options
     #Washing
     #Holding
     #Pre-cooling
     #Pre-Wash
     #PHS4d
     #PHS4h
     #PHSInt
     #HTrad
     #RTrad
     #FPTrad
     
def scenario_function(Cont_Scen_no,
                      #Intervention Strategies
                      Washing = False,
                      Washing_Optimized = False,
                      Holding = False,
                      Pre_Cooling = False,
                      PreS_Wash = False,
                      Sanitation = False,
                      
                      #Sampling Strategies.
                      PHS4d = False,
                      PHS4h= False,
                      PHSInt =False,
                      HSTrad = False,
                      RSTrad = False,
                      FPSTrad = False,
                      CSampling = False
                      ):
    ''' Docstring
    Select the scenarios from the arguements to run the scenario function. 
    Args: 
        
        Initial contamination (num): Number of CFU,
    
    Returns: 
        
        1. Outputs DF
        2. Progression DF
        3. Proportion Progression DF
    
    '''
    reload(MainModel3z)
    reload(Inputz)
    reload(SCInputz)
    reload(ScenCondz)
    reload(ContCondz)
    
    ScenCondz.Contamination_Scenario = Cont_Scen_no
        #Sampling Condition
    # Sampling Conditions, Baseline all conditions are off
    
    ScenCondz.Baseline_Sampling = 0  # all others must be 0 if this one is 1
    #PHS4d
    if PHS4d == True:
        ScenCondz.PH_Sampling = 1
        ScenCondz.PHS_4d = 1
        
    if PHS4h == True:
        ScenCondz.PH_Sampling = 1
        ScenCondz.PHS_4h = 1
    
    if PHSInt == True:
        ScenCondz.PH_Sampling = 1
        ScenCondz.PHS_Int = 1
    
    if HSTrad == True:
        ScenCondz.H_Sampling = 1
        ScenCondz.HS_Trad = 1
    
    if RSTrad == True:
        ScenCondz.R_Sampling = 1
    
    if FPSTrad == True:
        ScenCondz.FP_Sampling = 1
        ScenCondz.FPS_Trad = 1
    
    if CSampling == True:
        ScenCondz.C_Sampling = 1
    
    reload(SCInputz)
    
    
    
    #Management of System Control: 
    
    #Holding Time:
        #Potentially this scenario will only be applicable to the munre and irrigation scenario. 
    ScenCondz.Holding_Time= Holding #if we have 0-8 days or 2-8 days. 
    
    #Turning off Pre-Cooling
        #Precooling yes or not. 
    SCInputz.Pre_CoolingYN = Pre_Cooling 
    
    # Turning of Washing. 
        #Washing process yes or not. 
    SCInputz.Washing_YN = Washing
    
    SCInputz.Washing_Optimized =  Washing_Optimized
    
    #Harvest Pre-Wash: 
        #Harvest Pre-Wash yes or not
    SCInputz.Spray_WashYN =PreS_Wash
    
    #Sanitation:
    SCInputz.Sanitation_YN = Sanitation
    
        
    #Running The Model.
    Main_Mod_Outs = MainModel3z.F_MainLoop()
    
    #Getting the outputs from the function.
    OutputDF = Main_Mod_Outs[1]
    ProgDF = Main_Mod_Outs[0]
    PropProgDF = Main_Mod_Outs[2]

    return [OutputDF,ProgDF,PropProgDF]

#%% Important Functions for the Analysys

#COMPARING THE mean and 95% CI Reduction between two groups. 
import numpy as np
import math
import scipy.stats as st
import statsmodels.stats.api as sms

def mean_CI_ONE(Array):
    mean = Array.mean()
    CI = sms.DescrStatsW(Array).tconfint_mean()
    #log_mean = mean.log10()
    #CI_log = CI.log10()
    return [mean,CI]


def Calc_red(meanCI, treatments):
    list_1=[]
    for i in list(range(0,treatments)):
        mean = meanCI[i][0]
        list_1.append(((meanCI[0][0]-mean)/meanCI[0][0]))
    return list_1

def F_Sampling_Power (DF, Step_CFU_Acc, Step_CFU_Rej):
    Total_Contam = len(DF[(DF[Step_CFU_Acc]>0) | (DF[Step_CFU_Rej]>0) ])
    Total_Contam_Rej = len(DF[ (DF[Step_CFU_Rej]>0) ])
    Power = Total_Contam_Rej/Total_Contam
    return [Total_Contam, Total_Contam_Rej, Power]

def  CFU_Sampling_Stage (DF,Step_col_Name):
    Before_PHSInt= DF[Step_col_Name]
    return Before_PHSInt[Before_PHSInt>0].describe()

def Mean_Quantiles(df,columnname,q1,q2) : 
    mean_1 = df[columnname].mean()*100
    quantiles_1= df[columnname].quantile([q1,q2])
    return(mean_1, quantiles_1)

def Remove_Outlier_Indices(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    trueList = ~((df < (Q1 - 1.5 * IQR)) |(df > (Q3 + 1.5 * IQR)))
    return trueList


def Progression_DF_Melt(List_of_Outs):
    Column_Names = "BaselineNI BaselineAI Holding Precooling Washing PreSpray_Wash Sanitation".split()
    
    Index_1 = 0
    List_dfs = []
    for  i in List_of_Outs:
        #Progression Data. 
        i[1]["Type"] = Column_Names[Index_1]
        Melted_BNI_1 = i[1].melt(id_vars=['Type'])
        Index_1 = Index_1+1
        List_dfs.append(Melted_BNI_1)
    
    df = pd.concat(List_dfs)
    return df
    

def F_Outputs_Table(List_of_Outputs):
    rep = 0
    Columns_Final_Outs = [
        "Final_CFU_Acc_Portion_mean",
        "Final_CFU_Acc_Portion_5CI",
        "Final_CFU_Acc_Portion_95CI",
        "MeanComparison",
        "Final_CFU_Rej_Portion_mean",
        "Final_CFU_Rej_Portion_5CI",
        "Final_CFU_Rej_Portion_95CI",
        "Prevalence_Acc_Mean",
        "Prevalence_Acc_5CI",
        "Prevalence_Acc_95CI",
        "Prevalence_Comparison",
        "Prevalence_Rej_Mean",
        "Prevalence_Rej_5CI",
        "Prevalence_Rej_95CI",
        "Ratio_Product_accepted",
        "Pooled_CFU_g_mean",
        "Pooled_CFU_g_5CI",
        "Pooled_CFU_g_5CI",    
    ]
    
    Outputs_Df =pd.DataFrame(np.NaN, index= range(len(List_of_Outputs)), columns =Columns_Final_Outs)
    for i in List_of_Outputs:
        #Subseting Dataframe for the accepted
        Subset_Acc_NI_PHS4d= i[0][i[0]['C_Wei_Acc'] == 100_000].index
        Subset_Rej_NI_PHS4d= i[0][i[0]['C_Wei_Acc'] != 100_000].index
        
        
        #Final CFUs based on if accepted or rejected
            #Accepted
        Final_CFU_Acc_Portion_mean=i[1]["After CS Samp"].sum()
        Final_CFU_Acc_Portion_90CI=i[1]["After CS Samp"].quantile([0.05,0.95])
            #Rejected
        Final_CFU_Rej_Portion_mean=i[1][i[1].index.isin(Subset_Rej_NI_PHS4d)]["After CS Samp"].mean()
        Final_CFU_Rej_Portion_90CI=i[1][i[1].index.isin(Subset_Rej_NI_PHS4d)]["After CS Samp"].quantile([0.05,0.95])
        
        #Prevalence of contaminated packages 
            #Accepted
        Prevalence_Acc_Mean=i[2]["PropCont_A_CS_Whole"].mean()
        Prevalence_Acc_90CI=i[2]["PropCont_A_CS_Whole"].quantile([0.05,0.95])
            #Rejected
        Prevalence_Rej_Mean=i[2][i[2].index.isin(Subset_Rej_NI_PHS4d)]["PropCont_A_CS_Whole"].mean()
        Prevalence_Rej_90CI=i[2][i[2].index.isin(Subset_Rej_NI_PHS4d)]["PropCont_A_CS_Whole"].quantile([0.05,0.95])
            #Pooled Stuff
        Pooled_CFU_g_mean= ((i[1][i[1].index.isin(Subset_Acc_NI_PHS4d)]["After CS Samp"])/ (i[0][i[0].index.isin(Subset_Acc_NI_PHS4d)]["C_Wei_Acc"]*454)).mean()
        Pooled_CFU_g_90CI= ((i[1][i[1].index.isin(Subset_Acc_NI_PHS4d)]["After CS Samp"])/ (i[0][i[0].index.isin(Subset_Acc_NI_PHS4d)]["C_Wei_Acc"]*454)).quantile([0.05,0.95])
        
        #Ratio of product accepted all iterations (weight)
        i[0][i[0]['C_Wei_Acc'] == 50] = 0
        Ratio_Product_accepted=i[0]['C_Wei_Acc'].sum()/(100*100_000)
        
        Outputs_Df.at[rep,"Final_CFU_Acc_Portion_mean"] = Final_CFU_Acc_Portion_mean
        Outputs_Df.at[rep,"Final_CFU_Acc_Portion_5CI"] = Final_CFU_Acc_Portion_90CI.to_list()[0]
        Outputs_Df.at[rep,"Final_CFU_Acc_Portion_95CI"] = Final_CFU_Acc_Portion_90CI.to_list()[1]
        
        Outputs_Df.at[rep,"MeanComparison"] = Outputs_Df.at[rep,"Final_CFU_Acc_Portion_mean"] / Outputs_Df.at[0,"Final_CFU_Acc_Portion_mean"] 
        
        
        Outputs_Df.at[rep,"Final_CFU_Rej_Portion_mean"] = Final_CFU_Rej_Portion_mean
        Outputs_Df.at[rep,"Final_CFU_Rej_Portion_5CI"] = Final_CFU_Rej_Portion_90CI.to_list()[0]
        Outputs_Df.at[rep,"Final_CFU_Rej_Portion_95CI"] = Final_CFU_Rej_Portion_90CI.to_list()[1]
        
        Outputs_Df.at[rep,"Prevalence_Acc_Mean"] = Prevalence_Acc_Mean
        Outputs_Df.at[rep,"Prevalence_Acc_5CI"] = Prevalence_Acc_90CI.to_list()[0]
        Outputs_Df.at[rep,"Prevalence_Acc_95CI"] = Prevalence_Acc_90CI.to_list()[1]
        
        Outputs_Df.at[rep,"Prevalence_Comparison"] = Outputs_Df.at[rep,"Prevalence_Acc_Mean"] / Outputs_Df.at[0,"Prevalence_Acc_Mean"]
        
        Outputs_Df.at[rep,"Prevalence_Rej_Mean"] = Prevalence_Rej_Mean
        Outputs_Df.at[rep,"Prevalence_Rej_5CI"] = Prevalence_Rej_90CI.to_list()[0]
        Outputs_Df.at[rep,"Prevalence_Rej_95CI"] = Prevalence_Rej_90CI.to_list()[1]
        
        Outputs_Df.at[rep,"Ratio_Product_accepted"] = Ratio_Product_accepted
        
        Outputs_Df.at[rep,"Pooled_CFU_g_mean"] = Pooled_CFU_g_mean
        Outputs_Df.at[rep,"Pooled_CFU_g_5CI"] =  Pooled_CFU_g_90CI.to_list()[0]
        Outputs_Df.at[rep,"Pooled_CFU_g_95CI"] =  Pooled_CFU_g_90CI.to_list()[1]
        
        
        
        rep=rep+1
    return Outputs_Df
#%% Effect of Individual Interventions
'''
This chunk of code runs 6 systems talked on the effect of individual interventions section. 
'''
#Baseline No Intervention
Baseline_NI_1 =  scenario_function(Cont_Scen_no=1)
Baseline_NI_2 =  scenario_function(Cont_Scen_no=2)
Baseline_NI_3 =  scenario_function(Cont_Scen_no=3)

#Baseline All Intervention
Baseline_AI_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_AI_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_AI_3 =  scenario_function(Cont_Scen_no=3,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)

#Effect of Holding
Baseline_NI_Holding_1 =  scenario_function(Cont_Scen_no=1,Holding=True)
Baseline_NI_Holding_2 =  scenario_function(Cont_Scen_no=2,Holding=True)
Baseline_NI_Holding_3 =  scenario_function(Cont_Scen_no=3,Holding=True)

#Effect of PreCooling
Baseline_NI_Precooling_1 =  scenario_function(Cont_Scen_no=1,Pre_Cooling=True)
Baseline_NI_Precooling_2 =  scenario_function(Cont_Scen_no=2,Pre_Cooling=True)
Baseline_NI_Precooling_3 =  scenario_function(Cont_Scen_no=3,Pre_Cooling=True)

#Effect of Wash
Baseline_NI_Wash_1 =  scenario_function(Cont_Scen_no=1,Washing=True)
Baseline_NI_Wash_2 =  scenario_function(Cont_Scen_no=2,Washing=True)
Baseline_NI_Wash_3 =  scenario_function(Cont_Scen_no=3,Washing=True)

#Baseline_NI_Wash_Opt_1 =  scenario_function(Cont_Scen_no=1,Washing=True,Washing_Optimized =True)
#Baseline_NI_Wash_Opt_2 =  scenario_function(Cont_Scen_no=2,Washing=True,Washing_Optimized =True)
#Baseline_NI_Wash_Opt_3 =  scenario_function(Cont_Scen_no=3,Washing=True,Washing_Optimized =True)

#Harvest Wash
Baseline_NI_Sp_Wash_1 =  scenario_function(Cont_Scen_no=1,PreS_Wash=True)
Baseline_NI_Sp_Wash_2 =  scenario_function(Cont_Scen_no=2,PreS_Wash=True)
Baseline_NI_Sp_Wash_3 =  scenario_function(Cont_Scen_no=3,PreS_Wash=True)

#Baseline NI + Processing line Sanitation. 
Baseline_NI_PLS_1 = scenario_function(Cont_Scen_no=1,Sanitation=True)
Baseline_NI_PLS_2 = scenario_function(Cont_Scen_no=2,Sanitation=True)
Baseline_NI_PLS_3 = scenario_function(Cont_Scen_no=3,Sanitation=True)


#%%%

'''
This Chunk conducts data analysis for the effect of individual interventions
'''

# Data Analysis for getting the diffences between  individual interventions

#Scenario 1 Uniform Contamination
List_of_Outs_Ints_1 =    [Baseline_NI_1,
                         Baseline_AI_1,
                         Baseline_NI_Holding_1,
                         Baseline_NI_Precooling_1,
                         Baseline_NI_Wash_1,
                         Baseline_NI_Sp_Wash_1,
                         Baseline_NI_PLS_1
                         ]

Outputdf_INT_1 = F_Outputs_Table(List_of_Outs_Ints_1)


#Scenario 2 10% clustered Contamination
List_of_Outs_Ints_2 = [Baseline_NI_2,
                       Baseline_AI_2,
                       Baseline_NI_Holding_2,
                       Baseline_NI_Precooling_2,
                       Baseline_NI_Wash_2,
                       Baseline_NI_Sp_Wash_2,
                       Baseline_NI_PLS_2
                     ]

Outputdf_INT_2 = F_Outputs_Table(List_of_Outs_Ints_2)

#Scenario 1% cluster contamination.
List_of_Outs_Ints_3 = [Baseline_NI_3,
                       Baseline_AI_3,
                     Baseline_NI_Holding_3,
                     Baseline_NI_Precooling_3,
                     Baseline_NI_Wash_3,
                     Baseline_NI_Sp_Wash_3,
                     Baseline_NI_PLS_3
                     ]

Outputdf_INT_3 = F_Outputs_Table(List_of_Outs_Ints_3)


Melted_Prog_DF_NI_1 = Progression_DF_Melt(List_of_Outs = List_of_Outs_Ints_1)
Melted_Prog_DF_NI_2 = Progression_DF_Melt(List_of_Outs = List_of_Outs_Ints_2)
Melted_Prog_DF_NI_3 = Progression_DF_Melt(List_of_Outs = List_of_Outs_Ints_3)

Melted_Prog_DF_NI_1["Scenario"] = "Uniform"
Melted_Prog_DF_NI_2["Scenario"] = "1% Cluster"
Melted_Prog_DF_NI_3["Scenario"] = "10% Cluster"

##Exporting to CSV
Melted_Prog_DF_NI_1.to_csv(path_or_buf = "C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\Melted_Prog_DF_NI_1-R2T.csv")
Melted_Prog_DF_NI_2.to_csv(path_or_buf = "C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\Melted_Prog_DF_NI_2-R2T.csv")
Melted_Prog_DF_NI_3.to_csv(path_or_buf = "C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\Melted_Prog_DF_NI_3-R2T.csv")

#To generate plots see the Final Plots in R file.
#use section names Contamination progression.


########### --------Relative difference ----------- ################
Column_Names_INT = "BaselineNI BaselineAI Holding Precooling Washing PreSpray_Wash Sanitation".split()

Outputdf_INT_1["ScenarioN"] = Column_Names_INT
Outputdf_INT_2["ScenarioN"] = Column_Names_INT
Outputdf_INT_3["ScenarioN"] = Column_Names_INT

Outputdf_INT_1["Cont_Spread"] = "Uniform"
Outputdf_INT_2["Cont_Spread"] = "1% Cluster"
Outputdf_INT_3["Cont_Spread"] = "10% Cluster"

INT_Combined=pd.concat([Outputdf_INT_1,Outputdf_INT_3,Outputdf_INT_2])

#Compating two scenarios
def logredcomp  (df1,dfbaseline):
    return np.log10(df1[1]["After CS Samp"].mean() /dfbaseline[1]["After CS Samp"].mean())

def logredcomp_1  (df1,dfbaseline):
    a=  dfbaseline[1]["After CS Samp"]- df1[1]["After CS Samp"] 
    return a

def logredcomp_q  (df1,dfbaseline):
    q90= np.log10(df1[1]["After CS Samp"].quantile (0.95) /dfbaseline[1]["After CS Samp"].mean())
    q10= np.log10(df1[1]["After CS Samp"].quantile (0.15) /dfbaseline[1]["After CS Samp"].mean())
    return [q90,q10]

def red_CFU (df1,dfbaseline):
    sub =(df1[1]["After CS Samp"]-dfbaseline[1]["After CS Samp"])
    return [sub.mean(), sub.quantile(0.95), sub.quantile(0.05)]
    

#spray wash
logredcomp(Baseline_NI_Sp_Wash_1, Baseline_NI_1)
logredcomp(Baseline_NI_Sp_Wash_3, Baseline_NI_3)
logredcomp(Baseline_NI_Sp_Wash_2, Baseline_NI_2)

#wash
logredcomp(Baseline_NI_Wash_1, Baseline_NI_1)
logredcomp(Baseline_NI_Wash_3, Baseline_NI_3)
logredcomp(Baseline_NI_Wash_2, Baseline_NI_2)

#holding
logredcomp(Baseline_NI_Holding_1, Baseline_NI_1)
logredcomp(Baseline_NI_Holding_3, Baseline_NI_3)
logredcomp(Baseline_NI_Holding_2, Baseline_NI_2)

#precool
logredcomp(Baseline_NI_Precooling_1, Baseline_NI_1)
logredcomp(Baseline_NI_Precooling_3, Baseline_NI_3)
logredcomp(Baseline_NI_Precooling_2, Baseline_NI_2)

#Sanitation
logredcomp(Baseline_NI_PLS_1, Baseline_NI_1)
logredcomp(Baseline_NI_PLS_3, Baseline_NI_3)
logredcomp(Baseline_NI_PLS_2, Baseline_NI_2)

#All Interventions
logredcomp(Baseline_AI_1, Baseline_NI_1)
logredcomp(Baseline_AI_3, Baseline_NI_3)
logredcomp(Baseline_AI_2, Baseline_NI_2)


#%% Running the scenario Analysis 

### BASELINE  NO INTERVENTION ####

#Baseline Scenario No Intervention. 
Baseline_NI_1 =  scenario_function(Cont_Scen_no=1)
Baseline_NI_2 =  scenario_function(Cont_Scen_no=2)
Baseline_NI_3 =  scenario_function(Cont_Scen_no=3)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True)
Baseline_NI_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True)
Baseline_NI_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True)
Baseline_NI_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True)
Baseline_NI_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True)
Baseline_NI_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True)
Baseline_NI_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_HS_1=  scenario_function(Cont_Scen_no=1,HSTrad = True)
Baseline_NI_HS_2=  scenario_function(Cont_Scen_no=2,HSTrad = True)
Baseline_NI_HS_3=  scenario_function(Cont_Scen_no=3,HSTrad = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True)
Baseline_NI_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True)
Baseline_NI_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True)
Baseline_NI_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True)
Baseline_NI_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True)
Baseline_NI_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True)
Baseline_NI_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True)


### BASELINE  ALL INTERVENTIONS###

### Pre-Harvest Sampling 4 days.
#Baseline Scenario All Interventions.
Baseline_AI_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_AI_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_AI_3 =  scenario_function(Cont_Scen_no=3,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
 

#Baseline with intervention. 4 days pre-harvest sampling

Baseline_AI_PHS4d_1 =  scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4d = True)
Baseline_AI_PHS4d_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4d = True)
Baseline_AI_PHS4d_3 =  scenario_function(Cont_Scen_no=3,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4d = True)


### Pre-Harvest Sampling 4h
#Baseline with intervention. 4 hours pre-harvest sampling
Baseline_AI_PHS4h_1 =  scenario_function(Cont_Scen_no=1,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4h = True)
Baseline_AI_PHS4h_2 =  scenario_function(Cont_Scen_no=2,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4h = True)
Baseline_AI_PHS4h_3 =  scenario_function(Cont_Scen_no=3,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHS4h = True)

### Pre-Harvest Sampling Intense
#Baseline with intervention. Intense pre-harvest sampling
Baseline_AI_PHSInt_1 =  scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHSInt = True)
Baseline_AI_PHSInt_2 =  scenario_function(Cont_Scen_no=2,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHSInt = True)
Baseline_AI_PHSInt_3 =  scenario_function(Cont_Scen_no=3,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, PHSInt = True)

### Harvest Sampling Intense
#Baseline with intervention. Intense Harvest sampling
Baseline_AI_H_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, HSTrad= True)
Baseline_AI_H_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, HSTrad= True)
Baseline_AI_H_3 =  scenario_function(Cont_Scen_no=3,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, HSTrad= True)

### Receiving Sampling Intense
#Baseline with intervention. Intense Receiving sampling
Baseline_AI_R_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, RSTrad = True)
Baseline_AI_R_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, RSTrad = True)
Baseline_AI_R_3 =  scenario_function(Cont_Scen_no=3,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, RSTrad = True)

### FPS Sampling Intense
#Baseline with intervention. Intense Receiving sampling
Baseline_AI_FP_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, FPSTrad = True)
Baseline_AI_FP_2 =  scenario_function(Cont_Scen_no=2,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, FPSTrad = True)
Baseline_AI_FP_3 =  scenario_function(Cont_Scen_no=3,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, FPSTrad = True)

### FPS Sampling Intense
#Baseline with intervention. Intense Receiving sampling
Baseline_AI_CS_1 =  scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, CSampling = True)
Baseline_AI_CS_2 =  scenario_function(Cont_Scen_no=2,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, CSampling = True)
Baseline_AI_CS_3 =  scenario_function(Cont_Scen_no=3,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True, CSampling = True)


###BASELINE NI - WASHING  #####################
#Baseline Scenario No Intervention. 
Baseline_NI_W_1 =  scenario_function(Cont_Scen_no=1,Washing = True)
Baseline_NI_W_2 =  scenario_function(Cont_Scen_no=2,Washing = True)
Baseline_NI_W_3 =  scenario_function(Cont_Scen_no=3,Washing = True)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_W_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True,Washing = True)
Baseline_NI_W_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True,Washing = True)
Baseline_NI_W_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True,Washing = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_W_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True,Washing = True)
Baseline_NI_W_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True,Washing = True)
Baseline_NI_W_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True,Washing = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_W_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True,Washing = True)
Baseline_NI_W_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True,Washing = True)
Baseline_NI_W_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True,Washing = True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_W_H_1=  scenario_function(Cont_Scen_no=1,HSTrad = True,Washing = True)
Baseline_NI_W_H_2=  scenario_function(Cont_Scen_no=2,HSTrad = True,Washing = True)
Baseline_NI_W_H_3=  scenario_function(Cont_Scen_no=3,HSTrad = True,Washing = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_W_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True,Washing = True)
Baseline_NI_W_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True,Washing = True)
Baseline_NI_W_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True,Washing = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_W_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True,Washing = True)
Baseline_NI_W_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True,Washing = True)
Baseline_NI_W_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True,Washing = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_W_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True,Washing = True)
Baseline_NI_W_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True,Washing = True)
Baseline_NI_W_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True,Washing = True)


###########BASELINE NI -Holding  #####################
#Baseline Scenario No Intervention. 
Baseline_NI_H_1 =  scenario_function(Cont_Scen_no=1,Holding =  True)
Baseline_NI_H_2 =  scenario_function(Cont_Scen_no=2,Holding = True)
Baseline_NI_H_3 =  scenario_function(Cont_Scen_no=3,Holding =  True)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_H_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True,Holding = True)
Baseline_NI_H_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True,Holding = True)
Baseline_NI_H_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True,Holding = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_H_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True,Holding = True)
Baseline_NI_H_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True,Holding = True)
Baseline_NI_H_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True,Holding = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_H_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True,Holding = True)
Baseline_NI_H_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True,Holding = True)
Baseline_NI_H_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True,Holding = True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_H_H_1=  scenario_function(Cont_Scen_no=1,HSTrad = True,Holding = True)
Baseline_NI_H_H_2=  scenario_function(Cont_Scen_no=2,HSTrad = True,Holding = True)
Baseline_NI_H_H_3=  scenario_function(Cont_Scen_no=3,HSTrad = True,Holding = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_H_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True,Holding = True)
Baseline_NI_H_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True,Holding = True)
Baseline_NI_H_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True,Holding = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_H_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True,Holding = True)
Baseline_NI_H_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True,Holding = True)
Baseline_NI_H_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True,Holding = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_H_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True,Holding = True)
Baseline_NI_H_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True,Holding = True)
Baseline_NI_H_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True,Holding = True)

###############BASELINE NI -PreWash #####################
#Baseline Scenario No Intervention. 
Baseline_NI_PW_1 =  scenario_function(Cont_Scen_no=1, PreS_Wash=True)
Baseline_NI_PW_2 =  scenario_function(Cont_Scen_no=2, PreS_Wash=True)
Baseline_NI_PW_3 =  scenario_function(Cont_Scen_no=3, PreS_Wash=True)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_PW_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True, PreS_Wash=True)
Baseline_NI_PW_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True, PreS_Wash=True)
Baseline_NI_PW_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True, PreS_Wash=True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PW_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True, PreS_Wash=True)
Baseline_NI_PW_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True, PreS_Wash=True)
Baseline_NI_PW_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True, PreS_Wash=True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PW_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True, PreS_Wash=True)
Baseline_NI_PW_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True, PreS_Wash=True)
Baseline_NI_PW_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True, PreS_Wash=True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_PW_H_1=  scenario_function(Cont_Scen_no=1,HSTrad = True, PreS_Wash=True)
Baseline_NI_PW_H_2=  scenario_function(Cont_Scen_no=2,HSTrad = True, PreS_Wash=True)
Baseline_NI_PW_H_3=  scenario_function(Cont_Scen_no=3,HSTrad = True, PreS_Wash=True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PW_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True, PreS_Wash=True)
Baseline_NI_PW_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True, PreS_Wash=True)
Baseline_NI_PW_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True, PreS_Wash=True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PW_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True, PreS_Wash=True)
Baseline_NI_PW_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True, PreS_Wash=True)
Baseline_NI_PW_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True, PreS_Wash=True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PW_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True, PreS_Wash=True)
Baseline_NI_PW_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True, PreS_Wash=True)
Baseline_NI_PW_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True, PreS_Wash=True)

#############BASELINE NI -SANITATION  #####################
#Baseline Scenario No Intervention. 
Baseline_NI_S_1 =  scenario_function(Cont_Scen_no=1, Sanitation = True)
Baseline_NI_S_2 =  scenario_function(Cont_Scen_no=2, Sanitation = True)
Baseline_NI_S_3 =  scenario_function(Cont_Scen_no=3, Sanitation = True)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_S_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True, Sanitation = True)
Baseline_NI_S_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True, Sanitation = True)
Baseline_NI_S_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True, Sanitation = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_S_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True, Sanitation = True)
Baseline_NI_S_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True, Sanitation = True)
Baseline_NI_S_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True, Sanitation = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_S_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True, Sanitation = True)
Baseline_NI_S_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True, Sanitation = True)
Baseline_NI_S_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True, Sanitation = True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_S_H_1=  scenario_function(Cont_Scen_no=1,HSTrad = True, Sanitation = True)
Baseline_NI_S_H_2=  scenario_function(Cont_Scen_no=2,HSTrad = True, Sanitation = True)
Baseline_NI_S_H_3=  scenario_function(Cont_Scen_no=3,HSTrad = True, Sanitation = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_S_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True, Sanitation = True)
Baseline_NI_S_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True, Sanitation = True)
Baseline_NI_S_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True, Sanitation = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_S_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True, Sanitation = True)
Baseline_NI_S_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True, Sanitation = True)
Baseline_NI_S_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True, Sanitation = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_S_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True, Sanitation = True)
Baseline_NI_S_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True, Sanitation = True)
Baseline_NI_S_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True, Sanitation = True)


##############BASELINE NI -Precool#####################
#Baseline Scenario No Intervention. 
Baseline_NI_PC_1 =  scenario_function(Cont_Scen_no=1, Pre_Cooling = True)
Baseline_NI_PC_2 =  scenario_function(Cont_Scen_no=2, Pre_Cooling = True)
Baseline_NI_PC_3 =  scenario_function(Cont_Scen_no=3, Pre_Cooling = True)

### Pre-Harvest Sampling 4 days.
#Baseline no intervention. 4 days preharvest sampling
Baseline_NI_PC_PHS4d_1 =  scenario_function(Cont_Scen_no=1,PHS4d = True, Pre_Cooling = True)
Baseline_NI_PC_PHS4d_2 =  scenario_function(Cont_Scen_no=2,PHS4d = True, Pre_Cooling = True)
Baseline_NI_PC_PHS4d_3 =  scenario_function(Cont_Scen_no=3,PHS4d = True, Pre_Cooling = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PC_PHS4h_1 =  scenario_function(Cont_Scen_no=1,PHS4h = True, Pre_Cooling = True)
Baseline_NI_PC_PHS4h_2 =  scenario_function(Cont_Scen_no=2,PHS4h = True, Pre_Cooling = True)
Baseline_NI_PC_PHS4h_3 =  scenario_function(Cont_Scen_no=3,PHS4h = True, Pre_Cooling = True)

#Baseline no intervention. 4 hours preharvest sampling
Baseline_NI_PC_PHSInt_1 =  scenario_function(Cont_Scen_no=1,PHSInt = True, Pre_Cooling = True)
Baseline_NI_PC_PHSInt_2 =  scenario_function(Cont_Scen_no=2,PHSInt = True, Pre_Cooling = True)
Baseline_NI_PC_PHSInt_3 =  scenario_function(Cont_Scen_no=3,PHSInt = True, Pre_Cooling = True)

#Baseline no intervention Harvest Sampling Traditional
Baseline_NI_PC_H_1=  scenario_function(Cont_Scen_no=1,HSTrad = True, Pre_Cooling = True)
Baseline_NI_PC_H_2=  scenario_function(Cont_Scen_no=2,HSTrad = True, Pre_Cooling = True)
Baseline_NI_PC_H_3=  scenario_function(Cont_Scen_no=3,HSTrad = True, Pre_Cooling = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PC_R_1=  scenario_function(Cont_Scen_no=1,RSTrad =True, Pre_Cooling = True)
Baseline_NI_PC_R_2=  scenario_function(Cont_Scen_no=2,RSTrad =True, Pre_Cooling = True)
Baseline_NI_PC_R_3=  scenario_function(Cont_Scen_no=3,RSTrad =True, Pre_Cooling = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PC_FP_1=  scenario_function(Cont_Scen_no=1,FPSTrad =True, Pre_Cooling = True)
Baseline_NI_PC_FP_2=  scenario_function(Cont_Scen_no=2,FPSTrad =True, Pre_Cooling = True)
Baseline_NI_PC_FP_3=  scenario_function(Cont_Scen_no=3,FPSTrad =True, Pre_Cooling = True)

#Baseline no intervention Receiving Samplgin Traditional
Baseline_NI_PC_CS_1=  scenario_function(Cont_Scen_no=1,CSampling =True, Pre_Cooling = True)
Baseline_NI_PC_CS_2=  scenario_function(Cont_Scen_no=2,CSampling =True, Pre_Cooling = True)
Baseline_NI_PC_CS_3=  scenario_function(Cont_Scen_no=3,CSampling =True, Pre_Cooling = True)

#%% Data analysis

#step 1: Creation of output dataframes ###########
#No Intervention
List_of_Outs_NI_1 = [Baseline_NI_1,
                     Baseline_NI_PHS4d_1,
                     Baseline_NI_PHS4h_1,
                     Baseline_NI_PHSInt_1,
                     Baseline_NI_HS_1,
                     Baseline_NI_R_1,
                     Baseline_NI_FP_1,
                     Baseline_NI_CS_1
                     ]

List_of_Outs_NI_2 = [Baseline_NI_2,
                     Baseline_NI_PHS4d_2,
                     Baseline_NI_PHS4h_2,
                     Baseline_NI_PHSInt_2,
                     Baseline_NI_HS_2,
                     Baseline_NI_R_2,
                     Baseline_NI_FP_2,
                     Baseline_NI_CS_2
                     ]

List_of_Outs_NI_3 = [Baseline_NI_3,
                     Baseline_NI_PHS4d_3,
                     Baseline_NI_PHS4h_3,
                     Baseline_NI_PHSInt_3,
                     Baseline_NI_HS_3,
                     Baseline_NI_R_3,
                     Baseline_NI_FP_3,
                     Baseline_NI_CS_3
                     ]

Outputs_Df_NI_1=F_Outputs_Table(List_of_Outs_NI_1)  
Outputs_Df_NI_2=F_Outputs_Table(List_of_Outs_NI_2)
Outputs_Df_NI_3=F_Outputs_Table(List_of_Outs_NI_3)   

#All Interventions
List_of_Outs_AI_1 = [Baseline_AI_1,
                     Baseline_AI_PHS4d_1,
                     Baseline_AI_PHS4h_1,
                     Baseline_AI_PHSInt_1,
                     Baseline_AI_H_1,
                     Baseline_AI_R_1,
                     Baseline_AI_FP_1,
                     Baseline_AI_CS_1
                     ]

List_of_Outs_AI_2 = [Baseline_AI_2,
                     Baseline_AI_PHS4d_2,
                     Baseline_AI_PHS4h_2,
                     Baseline_AI_PHSInt_2,
                     Baseline_AI_H_2,
                     Baseline_AI_R_2,
                     Baseline_AI_FP_2,
                     Baseline_AI_CS_2
                     ]

List_of_Outs_AI_3 = [Baseline_AI_3,
                     Baseline_AI_PHS4d_3,
                     Baseline_AI_PHS4h_3,
                     Baseline_AI_PHSInt_3,
                     Baseline_AI_H_3,
                     Baseline_AI_R_3,
                     Baseline_AI_FP_3,
                     Baseline_AI_CS_3
                     ]

Outputs_Df_AI_1=F_Outputs_Table(List_of_Outs_AI_1)  
Outputs_Df_AI_2=F_Outputs_Table(List_of_Outs_AI_2)
Outputs_Df_AI_3=F_Outputs_Table(List_of_Outs_AI_3)  

#NI + WASHING
List_of_Outs_NI_W_1 = [Baseline_NI_W_1,
                     Baseline_NI_W_PHS4d_1,
                     Baseline_NI_W_PHS4h_1,
                     Baseline_NI_W_PHSInt_1,
                     Baseline_NI_W_H_1,
                     Baseline_NI_W_R_1,
                     Baseline_NI_W_FP_1,
                     Baseline_NI_W_CS_1
                     ]

List_of_Outs_NI_W_2 = [Baseline_NI_W_2,
                     Baseline_NI_W_PHS4d_2,
                     Baseline_NI_W_PHS4h_2,
                     Baseline_NI_W_PHSInt_2,
                     Baseline_NI_W_H_2,
                     Baseline_NI_W_R_2,
                     Baseline_NI_W_FP_2,
                     Baseline_NI_W_CS_2
                     ]

List_of_Outs_NI_W_3 = [Baseline_NI_W_3,
                     Baseline_NI_W_PHS4d_3,
                     Baseline_NI_W_PHS4h_3,
                     Baseline_NI_W_PHSInt_3,
                     Baseline_NI_W_H_3,
                     Baseline_NI_W_R_3,
                     Baseline_NI_W_FP_3,
                     Baseline_NI_W_CS_3
                     ]

Outputs_Df_NI_W_1=F_Outputs_Table(List_of_Outs_NI_W_1)  
Outputs_Df_NI_W_2=F_Outputs_Table(List_of_Outs_NI_W_2)
Outputs_Df_NI_W_3=F_Outputs_Table(List_of_Outs_NI_W_3)   



#No Intervention + Spray Washing
List_of_Outs_NI_PW_1 = [Baseline_NI_PW_1,
                     Baseline_NI_PW_PHS4d_1,
                     Baseline_NI_PW_PHS4h_1,
                     Baseline_NI_PW_PHSInt_1,
                     Baseline_NI_PW_H_1,
                     Baseline_NI_PW_R_1,
                     Baseline_NI_PW_FP_1,
                     Baseline_NI_PW_CS_1
                     ]

List_of_Outs_NI_PW_2 = [Baseline_NI_PW_2,
                     Baseline_NI_PW_PHS4d_2,
                     Baseline_NI_PW_PHS4h_2,
                     Baseline_NI_PW_PHSInt_2,
                     Baseline_NI_PW_H_2,
                     Baseline_NI_PW_R_2,
                     Baseline_NI_PW_FP_2,
                     Baseline_NI_PW_CS_2
                     ]

List_of_Outs_NI_PW_3 = [Baseline_NI_PW_3,
                     Baseline_NI_PW_PHS4d_3,
                     Baseline_NI_PW_PHS4h_3,
                     Baseline_NI_PW_PHSInt_3,
                     Baseline_NI_PW_H_3,
                     Baseline_NI_PW_R_3,
                     Baseline_NI_PW_FP_3,
                     Baseline_NI_PW_CS_3
                     ]

Outputs_Df_NI_PW_1=F_Outputs_Table(List_of_Outs_NI_PW_1)  
Outputs_Df_NI_PW_2=F_Outputs_Table(List_of_Outs_NI_PW_2)
Outputs_Df_NI_PW_3=F_Outputs_Table(List_of_Outs_NI_PW_3) 

#No Intervention + Holding
List_of_Outs_NI_H_1 = [Baseline_NI_H_1,
                     Baseline_NI_H_PHS4d_1,
                     Baseline_NI_H_PHS4h_1,
                     Baseline_NI_H_PHSInt_1,
                     Baseline_NI_H_H_1,
                     Baseline_NI_H_R_1,
                     Baseline_NI_H_FP_1,
                     Baseline_NI_H_CS_1
                     ]

List_of_Outs_NI_H_2 = [Baseline_NI_H_2,
                     Baseline_NI_H_PHS4d_2,
                     Baseline_NI_H_PHS4h_2,
                     Baseline_NI_H_PHSInt_2,
                     Baseline_NI_H_H_2,
                     Baseline_NI_H_R_2,
                     Baseline_NI_H_FP_2,
                     Baseline_NI_H_CS_2
                     ]

List_of_Outs_NI_H_3 = [Baseline_NI_H_3,
                     Baseline_NI_H_PHS4d_3,
                     Baseline_NI_H_PHS4h_3,
                     Baseline_NI_H_PHSInt_3,
                     Baseline_NI_H_H_3,
                     Baseline_NI_H_R_3,
                     Baseline_NI_H_FP_3,
                     Baseline_NI_H_CS_3
                     ]

Outputs_Df_NI_H_1=F_Outputs_Table(List_of_Outs_NI_H_1)  
Outputs_Df_NI_H_2=F_Outputs_Table(List_of_Outs_NI_H_2)
Outputs_Df_NI_H_3=F_Outputs_Table(List_of_Outs_NI_H_3)

#No Intervention + Sanitation
List_of_Outs_NI_S_1 = [Baseline_NI_S_1,
                     Baseline_NI_S_PHS4d_1,
                     Baseline_NI_S_PHS4h_1,
                     Baseline_NI_S_PHSInt_1,
                     Baseline_NI_S_H_1,
                     Baseline_NI_S_R_1,
                     Baseline_NI_S_FP_1,
                     Baseline_NI_S_CS_1
                     ]

List_of_Outs_NI_S_2 = [Baseline_NI_S_2,
                     Baseline_NI_S_PHS4d_2,
                     Baseline_NI_S_PHS4h_2,
                     Baseline_NI_S_PHSInt_2,
                     Baseline_NI_S_H_2,
                     Baseline_NI_S_R_2,
                     Baseline_NI_S_FP_2,
                     Baseline_NI_S_CS_2
                     ]

List_of_Outs_NI_S_3 = [Baseline_NI_S_3,
                     Baseline_NI_S_PHS4d_3,
                     Baseline_NI_S_PHS4h_3,
                     Baseline_NI_S_PHSInt_3,
                     Baseline_NI_S_H_3,
                     Baseline_NI_S_R_3,
                     Baseline_NI_S_FP_3,
                     Baseline_NI_S_CS_3
                     ]

Outputs_Df_NI_S_1=F_Outputs_Table(List_of_Outs_NI_S_1)  
Outputs_Df_NI_S_2=F_Outputs_Table(List_of_Outs_NI_S_2)
Outputs_Df_NI_S_3=F_Outputs_Table(List_of_Outs_NI_S_3)

#No Intervention + PRecooling
List_of_Outs_NI_PC_1 = [Baseline_NI_PC_1,
                     Baseline_NI_PC_PHS4d_1,
                     Baseline_NI_PC_PHS4h_1,
                     Baseline_NI_PC_PHSInt_1,
                     Baseline_NI_PC_H_1,
                     Baseline_NI_PC_R_1,
                     Baseline_NI_PC_FP_1,
                     Baseline_NI_PC_CS_1
                     ]

List_of_Outs_NI_PC_2 = [Baseline_NI_PC_2,
                     Baseline_NI_PC_PHS4d_2,
                     Baseline_NI_PC_PHS4h_2,
                     Baseline_NI_PC_PHSInt_2,
                     Baseline_NI_PC_H_2,
                     Baseline_NI_PC_R_2,
                     Baseline_NI_PC_FP_2,
                     Baseline_NI_PC_CS_2
                     ]

List_of_Outs_NI_PC_3 = [Baseline_NI_PC_3,
                     Baseline_NI_PC_PHS4d_3,
                     Baseline_NI_PC_PHS4h_3,
                     Baseline_NI_PC_PHSInt_3,
                     Baseline_NI_PC_H_3,
                     Baseline_NI_PC_R_3,
                     Baseline_NI_PC_FP_3,
                     Baseline_NI_PC_CS_3
                     ]

Outputs_Df_NI_PC_1=F_Outputs_Table(List_of_Outs_NI_PC_1)  
Outputs_Df_NI_PC_2=F_Outputs_Table(List_of_Outs_NI_PC_2)
Outputs_Df_NI_PC_3=F_Outputs_Table(List_of_Outs_NI_PC_3)



#%% Sampling Plan Power
    #To get power
def sampling_power(df,Step_Acc):
    return len(df[0][df[0][Step_Acc] ==0])/1000
    #to get number of detected
def sampling_power_2(df,Step_Acc):
    return len(df[0][df[0][Step_Acc] ==0])



Step_Acc_List = "PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc H_Wei_Acc R_Wei_Acc FP_Wei_Acc C_Wei_Acc".split()

######## No Intervention Powers ----------
Powers_NI = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_1[i],Step_Acc_List[i])
    Powers_NI.append(Power_1)   
Powers_NI

Powers_NI_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_1[i],Step_Acc_List[i])
    Powers_NI_N.append(Power_1)   
Powers_NI_N

#Scenario 2:  
Powers_NI_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_2[i],Step_Acc_List[i])
    Powers_NI_2.append(Power_1)   
Powers_NI_2

Powers_NI_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_2[i],Step_Acc_List[i])
    Powers_NI_N_2.append(Power_1)   
Powers_NI_N_2

#Scenario 3

Powers_NI_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_3[i],Step_Acc_List[i])
    Powers_NI_3.append(Power_1)
Powers_NI_3

Powers_NI_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_3[i],Step_Acc_List[i])
    Powers_NI_N_3.append(Power_1)
Powers_NI_N_3

###########All Intervention Powers: -------
 
Powers_AI = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_1[i],Step_Acc_List[i])
    Powers_AI.append(Power_1) 
Powers_AI


Powers_AI_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_1[i],Step_Acc_List[i])
    Powers_AI_N.append(Power_1)   
Powers_AI_N


#Scenario 2:    
Powers_AI_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_2[i],Step_Acc_List[i])
    Powers_AI_2.append(Power_1)   
Powers_AI_2


Powers_AI_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_2[i],Step_Acc_List[i])
    Powers_AI_N_2.append(Power_1)   
Powers_AI_N_2

#Scenario 3

Powers_AI_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_3[i],Step_Acc_List[i])
    Powers_AI_3.append(Power_1)   
Powers_AI_3

Powers_AI_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_3[i],Step_Acc_List[i])
    Powers_AI_N_3.append(Power_1)  
Powers_AI_N_3


######## No Intervention + WASHING ----------
Powers_NI_W = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_W_1[i],Step_Acc_List[i])
    Powers_NI_W.append(Power_1)   
Powers_NI_W

Powers_NI_W_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_W_1[i],Step_Acc_List[i])
    Powers_NI_W_N.append(Power_1)   
Powers_NI_W_N

#Scenario 2:  
Powers_NI_W_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_W_2[i],Step_Acc_List[i])
    Powers_NI_W_2.append(Power_1)   
Powers_NI_W_2

Powers_NI_W_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_W_2[i],Step_Acc_List[i])
    Powers_NI_W_N_2.append(Power_1)   
Powers_NI_W_N_2

#Scenario 3

Powers_NI_W_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_W_3[i],Step_Acc_List[i])
    Powers_NI_W_3.append(Power_1)
Powers_NI_W_3

Powers_NI_W_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_W_3[i],Step_Acc_List[i])
    Powers_NI_W_N_3.append(Power_1)
Powers_NI_W_N_3

######## No Intervention + Spray WASHING ----------
Powers_NI_PW = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PW_1[i],Step_Acc_List[i])
    Powers_NI_PW.append(Power_1)   
Powers_NI_PW

Powers_NI_PW_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PW_1[i],Step_Acc_List[i])
    Powers_NI_PW_N.append(Power_1)   
Powers_NI_PW_N

#Scenario 2:  
Powers_NI_PW_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PW_2[i],Step_Acc_List[i])
    Powers_NI_PW_2.append(Power_1)   
Powers_NI_PW_2

Powers_NI_PW_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PW_2[i],Step_Acc_List[i])
    Powers_NI_PW_N_2.append(Power_1)   
Powers_NI_PW_N_2

#Scenario 3

Powers_NI_PW_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PW_3[i],Step_Acc_List[i])
    Powers_NI_PW_3.append(Power_1)
Powers_NI_PW_3

Powers_NI_PW_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PW_3[i],Step_Acc_List[i])
    Powers_NI_PW_N_3.append(Power_1)
Powers_NI_PW_N_3

#NO Intervention + Holding #####################

Powers_NI_H = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_H_1[i],Step_Acc_List[i])
    Powers_NI_H.append(Power_1)   
Powers_NI_H
Powers_NI_H[0] = 0
Powers_NI_H_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_H_1[i],Step_Acc_List[i])
    Powers_NI_H_N.append(Power_1)   
Powers_NI_H_N

#Scenario 2:  
Powers_NI_H_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_H_2[i],Step_Acc_List[i])
    Powers_NI_H_2.append(Power_1)   
Powers_NI_H_2
Powers_NI_H_2[0] = 0

Powers_NI_H_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_H_2[i],Step_Acc_List[i])
    Powers_NI_H_N_2.append(Power_1)   
Powers_NI_H_N_2

#Scenario 3

Powers_NI_H_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_H_3[i],Step_Acc_List[i])
    Powers_NI_H_3.append(Power_1)
Powers_NI_H_3
Powers_NI_H_3[0] = 0

Powers_NI_H_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_H_3[i],Step_Acc_List[i])
    Powers_NI_H_N_3.append(Power_1)
Powers_NI_H_N_3


## No intervention + Sanitation ################


Powers_NI_S = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_S_1[i],Step_Acc_List[i])
    Powers_NI_S.append(Power_1)   
Powers_NI_S

Powers_NI_S_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_S_1[i],Step_Acc_List[i])
    Powers_NI_S_N.append(Power_1)   
Powers_NI_S_N

#Scenario 2:  
Powers_NI_S_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_S_2[i],Step_Acc_List[i])
    Powers_NI_S_2.append(Power_1)   
Powers_NI_S_2

Powers_NI_S_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_S_2[i],Step_Acc_List[i])
    Powers_NI_S_N_2.append(Power_1)   
Powers_NI_S_N_2

#Scenario 3

Powers_NI_S_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_S_3[i],Step_Acc_List[i])
    Powers_NI_S_3.append(Power_1)
Powers_NI_S_3

Powers_NI_S_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_S_3[i],Step_Acc_List[i])
    Powers_NI_S_N_3.append(Power_1)
Powers_NI_S_N_3

## No intervention + Precooling ################


Powers_NI_PC = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PC_1[i],Step_Acc_List[i])
    Powers_NI_PC.append(Power_1)   
Powers_NI_PC

Powers_NI_PC_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PC_1[i],Step_Acc_List[i])
    Powers_NI_PC_N.append(Power_1)   
Powers_NI_PC_N

#Scenario 2:  
Powers_NI_PC_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PC_2[i],Step_Acc_List[i])
    Powers_NI_PC_2.append(Power_1)   
Powers_NI_PC_2

Powers_NI_PC_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PC_2[i],Step_Acc_List[i])
    Powers_NI_PC_N_2.append(Power_1)   
Powers_NI_PC_N_2

#Scenario 3

Powers_NI_PC_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_PC_3[i],Step_Acc_List[i])
    Powers_NI_PC_3.append(Power_1)
Powers_NI_PC_3

Powers_NI_PC_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_PC_3[i],Step_Acc_List[i])
    Powers_NI_PC_N_3.append(Power_1)
Powers_NI_PC_N_3


Steps_For_Power = "None PHS4D PHS4H PHSInt HS RS FPS CS".split()
#Creating the output DF for powers
powers_df_1 = pd.DataFrame({
    "No Intervention": Powers_NI, #100%
    "All Intervention": Powers_AI, #100%
    "Washing Only": Powers_NI_W, #100%
    "Prewash Only": Powers_NI_PW, #100%
    "Holding Only": Powers_NI_H, #100%
    "Sanitation Only": Powers_NI_S, #100%
    "PreCooling Only": Powers_NI_PC, #100%
    "Type": Steps_For_Power 
    })

powers_df_1=powers_df_1.melt(id_vars=['Type'])
powers_df_1["ContS"] = "Random Uniform" 

powers_df_2 = pd.DataFrame({
   "No Intervention": Powers_NI_3, #10%
    "All Intervention": Powers_AI_3, #10%
    "Washing Only": Powers_NI_W_3, #100%
    "Prewash Only": Powers_NI_PW_3, #100%
    "Holding Only": Powers_NI_H_3, #100%
    "Sanitation Only": Powers_NI_S_3, #100%
    "PreCooling Only": Powers_NI_PC_3, #100%
    "Type": Steps_For_Power 
    })

powers_df_2=powers_df_2.melt(id_vars=['Type'])
powers_df_2["ContS"] = "10% Cluster" 


powers_df_3 = pd.DataFrame({   
    "No Intervention": Powers_NI_2, #1%
    "All Intervention": Powers_AI_2, #1%
    "Washing Only": Powers_NI_W_2, #100%
    "Prewash Only": Powers_NI_PW_2, #100%
    "Holding Only": Powers_NI_H_2, #100%
    "Sanitation Only": Powers_NI_S_2, #100%
    "PreCooling Only": Powers_NI_PC_2, #100%
    "Type": Steps_For_Power 
    })

powers_df_3=powers_df_3.melt(id_vars=['Type'])
powers_df_3["ContS"] = "1% Cluster" 

Power_For_Plot = pd.concat([powers_df_1,powers_df_2,powers_df_3])

Power_For_Plot.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\Powers.csv")


#%%
###### Contamination levels at sampling point ########


#Plot of contamination levels at sampling stage.
def Sampling_Loc_DF(list_1):
    a=list_1[1][1]["Bef Pre-Harvest Samp"]
    b=list_1[2][1]["Bef Pre-Harvest Samp"]
    c=list_1[3][1]["Bef Pre-Harvest Samp"]
    d=list_1[4][1]["Bef Harvest Samp"]
    e=list_1[5][1]["Bef Receiving Samp"]
    f=list_1[6][1]["Bef Final Prod S"]
    g=list_1[7][1]["Bef CS Samp"]
    
    df = pd.concat([a,b,c,d,e,f,g], axis = 1)
    Col_Names = ("Before PHS4d", "Before PHS4h","Before PHS Int","Before HS","Before RS","Before FPS","Before CS")
    df.columns = Col_Names
    
    df_melted = pd.melt(df)
    return df_melted

#No Intervention --- 
Sampling_Loc_NI_1 =Sampling_Loc_DF(List_of_Outs_NI_1)
Sampling_Loc_NI_2 =Sampling_Loc_DF(List_of_Outs_NI_2)
Sampling_Loc_NI_3 =Sampling_Loc_DF(List_of_Outs_NI_3)

Sampling_Loc_NI_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt = pd.concat([Sampling_Loc_NI_1,Sampling_Loc_NI_2,Sampling_Loc_NI_3])
Sampling_Loc_Melt["Process"] = "No Intervention"

Sampling_Loc_Melt.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\ContSamp-Final.csv")


#All Intervention ---
#Plot of contamination levels at sampling stage.
Sampling_Loc_AI_1 =Sampling_Loc_DF(List_of_Outs_AI_1)
Sampling_Loc_AI_2 =Sampling_Loc_DF(List_of_Outs_AI_2)
Sampling_Loc_AI_3 =Sampling_Loc_DF(List_of_Outs_AI_3)

Sampling_Loc_AI_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_AI_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_AI_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_AI = pd.concat([Sampling_Loc_AI_1,Sampling_Loc_AI_2,Sampling_Loc_AI_3])
Sampling_Loc_Melt_AI["Process"] = "All Intervention"

#Export CSV
Sampling_Loc_Melt_AI.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\ContSamp_AI-Final.csv")


#Washing Only---- 
Sampling_Loc_NI_W_1 =Sampling_Loc_DF(List_of_Outs_NI_W_1)
Sampling_Loc_NI_W_2 =Sampling_Loc_DF(List_of_Outs_NI_W_2)
Sampling_Loc_NI_W_3 =Sampling_Loc_DF(List_of_Outs_NI_W_3)

Sampling_Loc_NI_W_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_W_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_W_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_W = pd.concat([Sampling_Loc_NI_W_1,Sampling_Loc_NI_W_2,Sampling_Loc_NI_W_3])
Sampling_Loc_Melt_W["Process"] = "Washing Only"

#Prewas Only---- 
Sampling_Loc_NI_PW_1 =Sampling_Loc_DF(List_of_Outs_NI_PW_1)
Sampling_Loc_NI_PW_2 =Sampling_Loc_DF(List_of_Outs_NI_PW_2)
Sampling_Loc_NI_PW_3 =Sampling_Loc_DF(List_of_Outs_NI_PW_3)

Sampling_Loc_NI_PW_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_PW_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_PW_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_PW = pd.concat([Sampling_Loc_NI_PW_1,Sampling_Loc_NI_PW_2,Sampling_Loc_NI_PW_3])
Sampling_Loc_Melt_PW["Process"] = "PreWash Only"


#Holding Only Only---- 
Sampling_Loc_NI_H_1 =Sampling_Loc_DF(List_of_Outs_NI_H_1)
Sampling_Loc_NI_H_2 =Sampling_Loc_DF(List_of_Outs_NI_H_2)
Sampling_Loc_NI_H_3 =Sampling_Loc_DF(List_of_Outs_NI_H_3)

Sampling_Loc_NI_H_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_H_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_H_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_H = pd.concat([Sampling_Loc_NI_H_1,Sampling_Loc_NI_H_2,Sampling_Loc_NI_H_3])
Sampling_Loc_Melt_H["Process"] = "Holding Only"


#Sanitation Only Only---- 
Sampling_Loc_NI_S_1 =Sampling_Loc_DF(List_of_Outs_NI_S_1)
Sampling_Loc_NI_S_2 =Sampling_Loc_DF(List_of_Outs_NI_S_2)
Sampling_Loc_NI_S_3 =Sampling_Loc_DF(List_of_Outs_NI_S_3)

Sampling_Loc_NI_S_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_S_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_S_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_S = pd.concat([Sampling_Loc_NI_S_1,Sampling_Loc_NI_S_2,Sampling_Loc_NI_S_3])
Sampling_Loc_Melt_S["Process"] = "Sanitation Only"


#Precool Only Only---- 
Sampling_Loc_NI_PC_1 =Sampling_Loc_DF(List_of_Outs_NI_PC_1)
Sampling_Loc_NI_PC_2 =Sampling_Loc_DF(List_of_Outs_NI_PC_2)
Sampling_Loc_NI_PC_3 =Sampling_Loc_DF(List_of_Outs_NI_PC_3)

Sampling_Loc_NI_PC_1["Cont Scenario"] = "Random Uniform"
Sampling_Loc_NI_PC_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_PC_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_PC = pd.concat([Sampling_Loc_NI_PC_1,Sampling_Loc_NI_PC_2,Sampling_Loc_NI_PC_3])
Sampling_Loc_Melt_PC["Process"] = "PreCooling Only"

Contam_Samp_Point =pd.concat([Sampling_Loc_Melt,Sampling_Loc_Melt_AI,Sampling_Loc_Melt_W,
           Sampling_Loc_Melt_PW,Sampling_Loc_Melt_H,Sampling_Loc_Melt_S,Sampling_Loc_Melt_PC])

Contam_Samp_Point.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\ContSampPoint.csv")


#%% Consumer Exposure

#AI -----
Outputs_Df_AI_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_1["MeanComparison"]
Outputs_Df_AI_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_1["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_1["Cont_Spread"] = "Random"
A1=Outputs_Df_AI_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_AI_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_2["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_2["Cont_Spread"] = "1% Cluster"
A2=Outputs_Df_AI_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison", "MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_AI_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_3["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_3["Cont_Spread"] = "10% Cluster"
A3=Outputs_Df_AI_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI ------
Outputs_Df_NI_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_1["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_1["Cont_Spread"] = "Random"
N1=Outputs_Df_NI_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_2["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_2["Cont_Spread"] = "1% Cluster"
N2=Outputs_Df_NI_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_3["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_3["Cont_Spread"] = "10% Cluster"
N3=Outputs_Df_NI_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Wash
Outputs_Df_NI_W_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_W_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_W_1["Baseline_Scenario"] = "Wash Only"
Outputs_Df_NI_W_1["Cont_Spread"] = "Random"
N1_W=Outputs_Df_NI_W_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_W_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_W_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_W_2["Baseline_Scenario"] = "Wash Only"
Outputs_Df_NI_W_2["Cont_Spread"] = "1% Cluster"
N2_W=Outputs_Df_NI_W_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_W_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_W_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_W_3["Baseline_Scenario"] = "Wash Only"
Outputs_Df_NI_W_3["Cont_Spread"] = "10% Cluster"
N3_W=Outputs_Df_NI_W_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Prewash
Outputs_Df_NI_PW_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_1["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_1["Cont_Spread"] = "Random"
N1_W=Outputs_Df_NI_PW_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PW_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_2["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_2["Cont_Spread"] = "1% Cluster"
N2_W=Outputs_Df_NI_PW_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PW_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_3["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_3["Cont_Spread"] = "10% Cluster"
N3_W=Outputs_Df_NI_PW_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Prewash
Outputs_Df_NI_PW_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_1["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_1["Cont_Spread"] = "Random"
N1_PW=Outputs_Df_NI_PW_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PW_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_2["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_2["Cont_Spread"] = "1% Cluster"
N2_PW=Outputs_Df_NI_PW_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PW_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PW_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PW_3["Baseline_Scenario"] = "Prewash Only"
Outputs_Df_NI_PW_3["Cont_Spread"] = "10% Cluster"
N3_PW=Outputs_Df_NI_PW_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Holding
Outputs_Df_NI_H_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_H_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_H_1["Baseline_Scenario"] = "Holding Only"
Outputs_Df_NI_H_1["Cont_Spread"] = "Random"
N1_H=Outputs_Df_NI_H_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_H_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_H_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_H_2["Baseline_Scenario"] = "Holding Only"
Outputs_Df_NI_H_2["Cont_Spread"] = "1% Cluster"
N2_H=Outputs_Df_NI_H_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_H_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_H_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_H_3["Baseline_Scenario"] = "Holding Only"
Outputs_Df_NI_H_3["Cont_Spread"] = "10% Cluster"
N3_H=Outputs_Df_NI_H_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Sanitation
Outputs_Df_NI_S_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_S_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_S_1["Baseline_Scenario"] = "Sanitation Only"
Outputs_Df_NI_S_1["Cont_Spread"] = "Random"
N1_S=Outputs_Df_NI_S_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_S_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_S_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_S_2["Baseline_Scenario"] = "Sanitation Only"
Outputs_Df_NI_S_2["Cont_Spread"] = "1% Cluster"
N2_S=Outputs_Df_NI_S_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_S_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_S_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_S_3["Baseline_Scenario"] = "Sanitation Only"
Outputs_Df_NI_S_3["Cont_Spread"] = "10% Cluster"
N3_S=Outputs_Df_NI_S_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI _ Precooling
Outputs_Df_NI_PC_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PC_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PC_1["Baseline_Scenario"] = "PreCooling Only"
Outputs_Df_NI_PC_1["Cont_Spread"] = "Random"
N1_PC=Outputs_Df_NI_PC_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PC_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PC_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PC_2["Baseline_Scenario"] = "PreCooling Only"
Outputs_Df_NI_PC_2["Cont_Spread"] = "1% Cluster"
N2_PC=Outputs_Df_NI_PC_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_PC_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_PC_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_PC_3["Baseline_Scenario"] = "PreCooling Only"
Outputs_Df_NI_PC_3["Cont_Spread"] = "10% Cluster"
N3_PC=Outputs_Df_NI_PC_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]



Exposure_Chart_df= pd.concat([A1,A2,A3,N1,N2,N3,N1_W,N2_W,N3_W,N1_PW,N2_PW,N3_PW,N1_H,N2_H,N3_H,N1_S,N2_S,N3_S,N1_PC,N2_PC,N3_PC ])

Exposure_Chart_df.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\Exposure_Chart_df-Final2.csv")

#%%END
#%%
#%%
















################### Anything below this is for reference ############################
























#%%




#%%
#BASELINE NO INTERVETIONS -------------------------------------------------------

List_of_Outs_NI_1 = [Baseline_NI_1,
                     Baseline_NI_PHS4d_1,
                     Baseline_NI_PHS4h_1,
                     Baseline_NI_PHSInt_1,
                     Baseline_NI_H_1,
                     Baseline_NI_R_1,
                     Baseline_NI_FP_1,
                     Baseline_NI_CS_1
                     ]

List_of_Outs_NI_2 = [Baseline_NI_2,
                     Baseline_NI_PHS4d_2,
                     Baseline_NI_PHS4h_2,
                     Baseline_NI_PHSInt_2,
                     Baseline_NI_H_2,
                     Baseline_NI_R_2,
                     Baseline_NI_FP_2,
                     Baseline_NI_CS_2
                     ]

List_of_Outs_NI_3 = [Baseline_NI_3,
                     Baseline_NI_PHS4d_3,
                     Baseline_NI_PHS4h_3,
                     Baseline_NI_PHSInt_3,
                     Baseline_NI_H_3,
                     Baseline_NI_R_3,
                     Baseline_NI_FP_3,
                     Baseline_NI_CS_3
                     ]

#Names
Col_Names_NI = "BaselineNI PHS4D_NI PHS4H_NI PHSInt_NI HTrad_NI RSTrad_NI FPSTrad_NI CS_NI".split()
#Using function   
     
Outputs_Df_NI_1=F_Outputs_Table(List_of_Outs_NI_1)  
Outputs_Df_NI_2=F_Outputs_Table(List_of_Outs_NI_2)
Outputs_Df_NI_3=F_Outputs_Table(List_of_Outs_NI_3)    

#Sampling Results:

    #To get power
def sampling_power(df,Step_Acc):
    return len(df[0][df[0][Step_Acc] ==0])/10000
    #to get number of detected
def sampling_power_2(df,Step_Acc):
    return len(df[0][df[0][Step_Acc] ==0])

Step_Acc_List = "PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc H_Wei_Acc R_Wei_Acc FP_Wei_Acc C_Wei_Acc".split()

Powers_NI = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_1[i],Step_Acc_List[i])
    Powers_NI.append(Power_1)
    
Powers_NI

Powers_NI_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_1[i],Step_Acc_List[i])
    Powers_NI_N.append(Power_1)
    
Powers_NI_N

#Scenario 2: 
    
Powers_NI_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_2[i],Step_Acc_List[i])
    Powers_NI_2.append(Power_1)
    
Powers_NI_2

Powers_NI_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_2[i],Step_Acc_List[i])
    Powers_NI_N_2.append(Power_1)
    
Powers_NI_N_2

#Scenario 3

Powers_NI_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_NI_3[i],Step_Acc_List[i])
    Powers_NI_3.append(Power_1)
    
Powers_NI_3

Powers_NI_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_NI_3[i],Step_Acc_List[i])
    Powers_NI_N_3.append(Power_1)
    
Powers_NI_N_3


#Plot of contamination levels at sampling stage.
def Sampling_Loc_DF(list_1):
    a=list_1[1][1]["Bef Pre-Harvest Samp"]
    b=list_1[2][1]["Bef Pre-Harvest Samp"]
    c=list_1[3][1]["Bef Pre-Harvest Samp"]
    d=list_1[4][1]["Bef Harvest Samp"]
    e=list_1[5][1]["Bef Receiving Samp"]
    f=list_1[6][1]["Bef Final Prod S"]
    g=list_1[7][1]["Bef CS Samp"]
    
    df = pd.concat([a,b,c,d,e,f,g], axis = 1)
    Col_Names = ("Before PHS4d", "Before PHS4h","Before PHS Int","Before HS","Before RS","Before FPS","Before CS")
    df.columns = Col_Names
    
    df_melted = pd.melt(df)
    return df_melted

Sampling_Loc_NI_1 =Sampling_Loc_DF(List_of_Outs_NI_1)
Sampling_Loc_NI_2 =Sampling_Loc_DF(List_of_Outs_NI_2)
Sampling_Loc_NI_3 =Sampling_Loc_DF(List_of_Outs_NI_3)

Sampling_Loc_NI_1["Cont Scenario"] = "Random"
Sampling_Loc_NI_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_NI_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt = pd.concat([Sampling_Loc_NI_1,Sampling_Loc_NI_2,Sampling_Loc_NI_3])

#Exporting to CSV
Sampling_Loc_Melt.to_csv("C:\\Users\\Gustavo Reyes\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\ContSamp-Final.csv")

H=sns.catplot(x="variable", y="value", hue = "Cont Scenario", kind = "bar" ,
data=Sampling_Loc_Melt)
plt.xlabel("Sampling Stage")
plt.ylabel("Total CFUs in System")
plt.yscale('log')
plt.yticks([1,10,100,1000,10000,100000])
plt.title("No Intervention: System contamination before sampling steps")
plt.xticks(rotation=-90)

Outputs_Df_NI_1["ScenarioN"] = Col_Names_NI
Outputs_Df_NI_2["ScenarioN"] = Col_Names_NI
Outputs_Df_NI_3["ScenarioN"] = Col_Names_NI

Outputs_Df_NI_1["Cont_Spread"] = "Uniform"
Outputs_Df_NI_2["Cont_Spread"] = "1% Cluster"
Outputs_Df_NI_3["Cont_Spread"] = "10% Cluster"



#Bar Chart for Relative Difference: 
sns.barplot(x ="MeanComparison" , y = "ScenarioN", data = Outputs_Df_NI_1)
plt.xlabel("Relative Difference")
plt.title("No-Intervention: Uniform")

#Bar Chart for Relative Difference: 
sns.barplot(x ="MeanComparison" , y = "ScenarioN", data = Outputs_Df_NI_3)
plt.xlabel("Relative Difference")
plt.title("No-Intervention: 10% Cluster")

NI_Combined=pd.concat([Outputs_Df_NI_1,Outputs_Df_NI_3,Outputs_Df_NI_2])


#Bar Chart for Relative Difference: 
sns.barplot(x ="MeanComparison" , y = "ScenarioN", col = "Cont_Spread", data = NI_Combined)
plt.xlabel("Relative Difference")
plt.title("No-Intervention: 10% Cluster")

H= sns.catplot(x ="MeanComparison", y = "ScenarioN", col="Cont_Spread",
                data=NI_Combined, kind="bar",
                height=4, aspect=1)
H.set_axis_labels("Relative Difference", "Sampling Plan Scenario")


#%% New All intervention analysis

List_of_Outs_AI_1 = [Baseline_AI_1,
                     Baseline_AI_PHS4d_1,
                     Baseline_AI_PHS4h_1,
                     Baseline_AI_PHSInt_1,
                     Baseline_AI_H_1,
                     Baseline_AI_R_1,
                     Baseline_AI_FP_1,
                     Baseline_AI_CS_1
                     ]

List_of_Outs_AI_2 = [Baseline_AI_2,
                     Baseline_AI_PHS4d_2,
                     Baseline_AI_PHS4h_2,
                     Baseline_AI_PHSInt_2,
                     Baseline_AI_H_2,
                     Baseline_AI_R_2,
                     Baseline_AI_FP_2,
                     Baseline_AI_CS_2
                     ]

List_of_Outs_AI_3 = [Baseline_AI_3,
                     Baseline_AI_PHS4d_3,
                     Baseline_AI_PHS4h_3,
                     Baseline_AI_PHSInt_3,
                     Baseline_AI_H_3,
                     Baseline_AI_R_3,
                     Baseline_AI_FP_3,
                     Baseline_AI_CS_3
                     ]

#Names
Col_Names_AI = "BaselineAI PHS4D_AI PHS4H_AI PHSInt_AI HTrad_AI RSTrad_AI FPSTrad_AI CS_AI".split()
#Using function   
     
Outputs_Df_AI_1=F_Outputs_Table(List_of_Outs_AI_1)  
Outputs_Df_AI_2=F_Outputs_Table(List_of_Outs_AI_2)
Outputs_Df_AI_3=F_Outputs_Table(List_of_Outs_AI_3)    


##For outputs for new chart
Outputs_Df_AI_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_1["MeanComparison"]
Outputs_Df_AI_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_1["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_1["Cont_Spread"] = "Random"
A1=Outputs_Df_AI_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_AI_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_2["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_2["Cont_Spread"] = "1% Cluster"
A2=Outputs_Df_AI_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison", "MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_AI_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_AI_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_AI_3["Baseline_Scenario"] = "All Intervention"
Outputs_Df_AI_3["Cont_Spread"] = "10% Cluster"
A3=Outputs_Df_AI_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

#NI
Outputs_Df_NI_1["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_1["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_1["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_1["Cont_Spread"] = "Random"
N1=Outputs_Df_NI_1[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_2["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_2["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_2["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_2["Cont_Spread"] = "1% Cluster"
N2=Outputs_Df_NI_2[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Outputs_Df_NI_3["Final_CFU_Acc_Portion_mean"]
Outputs_Df_NI_3["Sampling_Plan"] = "Baseline PHS4D PHS4H PHSInt HTrad RSTrad FPSTrad CS".split()
Outputs_Df_NI_3["Baseline_Scenario"] = "No Intervention"
Outputs_Df_NI_3["Cont_Spread"] = "10% Cluster"
N3=Outputs_Df_NI_3[["Final_CFU_Acc_Portion_mean","Prevalence_Acc_Mean","Prevalence_Comparison","MeanComparison","Cont_Spread","Sampling_Plan","Baseline_Scenario"]]

Exposure_Chart_df= pd.concat([A1,A2,A3,N1,N2,N3])

Exposure_Chart_df.to_csv("C:\\Users\\Gustavo Reyes\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Exposure_Chart_df-Final.csv")




#Sampling Results:

def sampling_power(df,Step_Acc):
    return len(df[0][df[0][Step_Acc] ==0])/10000

Step_Acc_List = "PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc PH_Wei_Acc H_Wei_Acc R_Wei_Acc FP_Wei_Acc C_Wei_Acc".split()

Powers_AI = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_1[i],Step_Acc_List[i])
    Powers_AI.append(Power_1)
    
Powers_AI


Powers_AI_N = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_1[i],Step_Acc_List[i])
    Powers_AI_N.append(Power_1)
    
Powers_AI_N


#Scenario 2: 
    
Powers_AI_2 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_2[i],Step_Acc_List[i])
    Powers_AI_2.append(Power_1)
    
Powers_AI_2


Powers_AI_N_2 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_2[i],Step_Acc_List[i])
    Powers_AI_N_2.append(Power_1)
    
Powers_AI_N_2

#Scenario 3

Powers_AI_3 = []
for i in list(range(8)):
    Power_1= sampling_power(List_of_Outs_AI_3[i],Step_Acc_List[i])
    Powers_AI_3.append(Power_1)
    
Powers_AI_3


Powers_AI_N_3 = []
for i in list(range(8)):
    Power_1= sampling_power_2(List_of_Outs_AI_3[i],Step_Acc_List[i])
    Powers_AI_N_3.append(Power_1)
    
Powers_AI_N_3

powers_df = pd.DataFrame({
    "1NI": Powers_NI,
    "1AI": Powers_AI,
    "2NI": Powers_NI_3,
    "2AI": Powers_AI_3,
    "3NI": Powers_NI_2,
    "3AI": Powers_AI_2
    })


#Plot of contamination levels at sampling stage.
Sampling_Loc_AI_1 =Sampling_Loc_DF(List_of_Outs_AI_1)
Sampling_Loc_AI_2 =Sampling_Loc_DF(List_of_Outs_AI_2)
Sampling_Loc_AI_3 =Sampling_Loc_DF(List_of_Outs_AI_3)

Sampling_Loc_AI_1["Cont Scenario"] = "Random"
Sampling_Loc_AI_2["Cont Scenario"] = "1% Cluster"
Sampling_Loc_AI_3["Cont Scenario"] = "10% Cluster"

Sampling_Loc_Melt_AI = pd.concat([Sampling_Loc_AI_1,Sampling_Loc_AI_2,Sampling_Loc_AI_3])

#Export CSV
Sampling_Loc_Melt_AI.to_csv("C:\\Users\\Gustavo Reyes\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\ContSamp_AI-Final.csv")

H=sns.catplot(x="variable", y="value", hue = "Cont Scenario", kind = "bar" ,
data=Sampling_Loc_Melt)
plt.xlabel("Sampling Stage")
plt.ylabel("Total CFUs in System")
plt.yscale('log')
plt.yticks([1,10,100,1000,10000,100000])
plt.title("All Interventions: System contamination before sampling steps")
plt.xticks(rotation=-90)



#Noramlized Contamination
Power_AI_PHS4d = 1-0.127 
Power_AI_FP = 1-0

Subset_Acc_AI_PHS4d= Baseline_AI_PHS4d_1[0][Baseline_AI_PHS4d_1[0]['PH_Wei_Acc'] == 100_000].index
Subset_Rej_AI_PHS4d= Baseline_AI_PHS4d_1[0][Baseline_AI_PHS4d_1[0]['PH_Wei_Acc'] != 100_000].index


Subset_Acc_AI_FP= Baseline_AI_FP_1[0][Baseline_AI_FP_1[0]['FP_Wei_Acc'] == 100_000].index
Subset_Rej_AI_FP= Baseline_AI_FP_1[0][Baseline_AI_FP_1[0]['FP_Wei_Acc'] != 100_000].index

#Final CFUs based on if accepted or rejected

#Accepted
(Baseline_AI_PHS4d_1[1][Baseline_AI_PHS4d_1[1].index.isin(Subset_Acc_AI_PHS4d)]["After CS Samp"]*Power_AI_PHS4d).mean()

(Baseline_AI_FP_1[1][Baseline_AI_FP_1[1].index.isin(Subset_Acc_AI_FP)]["After CS Samp"]*Power_AI_FP).mean()


Final_CFU_Acc_Portion_90CI=i[1]["After CS Samp"].quantile([0.05,0.95])


###
Outputs_Df_AI_1["ScenarioN"] = Col_Names_AI
Outputs_Df_AI_2["ScenarioN"] = Col_Names_AI
Outputs_Df_AI_3["ScenarioN"] = Col_Names_AI

Outputs_Df_AI_1["Cont_Spread"] = "Uniform"
Outputs_Df_AI_2["Cont_Spread"] = "1% Cluster"
Outputs_Df_AI_3["Cont_Spread"] = "10% Cluster"

AI_Combined=pd.concat([Outputs_Df_AI_1,Outputs_Df_AI_3,Outputs_Df_AI_2])



H= sns.catplot(x ="MeanComparison", y = "ScenarioN", col="Cont_Spread",
                data=AI_Combined, kind="bar",
                height=4, aspect=1)
H.set_axis_labels("Relative Difference", "Sampling Plan Scenario")


#%%    
Outputs_Df.columns
Outputs_Df["Treatments"] = Col_Names_NI

#Exploring the relationships
sns.scatterplot(data=Outputs_Df, x="Final_CFU_Acc_Portion_mean", y="Ratio_Product_accepted", hue="Treatments")

#Pooled CFU_g, 







#Creating list of Final Contmainations
List_Final_Cont_NI = [x["Final Product Facility"] for x in List_of_Progs_NI]
Col_Names_NI = "BaselineNI PHS4D_NI PHS4H_NI PHSInt_NI HTrad_NI RSTrad_NI FPSTrad_NI".split()
#Creating and Melting the Dataframe of Final Contminations
DF_Final_Cont_NI = pd.concat(List_Final_Cont_NI, axis = 1)
DF_Final_Cont_NI.columns = Col_Names_NI
DF_Final_Cont_NI_melted = DF_Final_Cont_NI.melt()

#Creating Cat plot bar or box, change arguement
H=sns.catplot(x="variable", y="value", kind = "box" ,
            data=DF_Final_Cont_NI_melted)
plt.xlabel("Intervention")
plt.ylabel("Total CFUs at Finished Product")
plt.yscale('log')
plt.title("CFU Final Contamination: Baseline No INterventions")
plt.xticks(rotation=70)



#Means compared to each other. 
Final_Cont_NI_MeanCI= [mean_CI_ONE(x) for x in List_Final_Cont_NI]

Reduction_Cont_NI = Calc_red(meanCI =Final_Cont_NI_MeanCI,treatments= 7)

DF_Reduction_Cont_NI = pd.DataFrame({"Treatment": Col_Names_NI,
                             "Reduction": Reduction_Cont_NI} )

DF_Reduction_Cont_NI=DF_Reduction_Cont_NI.sort_values('Reduction',ascending=False).reset_index()

chart = sns.barplot(data = DF_Reduction_Cont_NI, x = "Treatment", y = "Reduction", order=DF_Reduction_Cont_NI['Treatment'])
chart.bar_label(chart.containers[0])
plt.xlabel("Intervention")
plt.ylabel("Percent Reduction from Baseline NI")



#Sampling Results: 

#PHS4d
F_Sampling_Power(Baseline_NI_PHS4d_1[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_PHS4d_1[1],"Bef Pre-Harvest Samp")
#PHS4h
F_Sampling_Power(Baseline_NI_PHS4h[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_PHS4h[1],"Bef Pre-Harvest Samp")
#PHSInt    
F_Sampling_Power(Baseline_NI_PHSInt[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_PHSInt[1],"Bef Pre-Harvest Samp")
#HS   
F_Sampling_Power(Baseline_NI_H[0],"H_CFU_Acc","H_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_H[1],"Bef Harvest Samp")
#RS
F_Sampling_Power(Baseline_NI_R[0],"R_CFU_Acc","R_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_R[1],"Bef Receiving Samp")
#FPS
F_Sampling_Power(Baseline_NI_FP[0],"FP_CFU_Acc","FP_CFU_Rej")
CFU_Sampling_Stage (Baseline_NI_FP[1],'Bef Final Prod S')

#Final CFUs for all. 
#prints the mean total CFUs and the CIs
Final_Cont_NI_MeanCI

#Getting the Prevalence of Contaminated Packages.

Mean_Quantiles(df = Baseline_NI[2],columnname ="PropCont_A_FP_Whole" ,q1 = 0.05,q2 = 0.95)  

#Proportion of packages contaminated
Intervention_Props_NI = [Baseline_NI[2],
                     Baseline_NI_PHS4d[2],
                     Baseline_NI_PHS4h[2],
                     Baseline_NI_PHSInt[2],
                     Baseline_NI_H[2],
                     Baseline_NI_R[2],
                     Baseline_NI_FP[2]
                     ]

#Creating dataframe of final contamination for every intervention. 
Intervention_PropCont_List_NI = [x["PropCont_A_FP_Whole"] for x in Intervention_Props_NI]

Prop_Cont_NI = pd.concat(Intervention_PropCont_List_NI, axis = 1)
Prop_Cont_NI.columns = Col_Names_NI
Prop_Cont_NI_melted = Prop_Cont_NI.melt()

#Plotting the bar graph or boxplot of the differences. 
H=sns.catplot(x="variable", y="value", kind = "box" ,
            data=Prop_Cont_NI_melted)
plt.xlabel("Intervention")
plt.ylabel("Proportion Contaminated at Final Product")
plt.xticks(rotation=70)

Mean_Props_NI= [mean_CI_ONE(x) for x in Intervention_PropCont_List_NI]

#Histplot
h=sns.displot( data =Prop_Cont_NI_melted, 
            x = "value" , 
            col = "variable", 
            col_wrap=3,
             stat = "count",
             bins = 30,
            facet_kws=dict(sharey=False,sharex= False))
plt.suptitle("Proportion of Paclages Contminated",) 

def specs(x, **kwargs):
    plt.axvline(x.mean(), c='red', ls='-', lw=2.5)
    plt.axvline(x.median(), c='orange', ls='--', lw=2.5)

h.map(specs,"value" )

 

### Baseline All Interventions, good system.

#Creating list of contamination progression
List_of_Progs_AI = [Baseline_AI[1],
                     Baseline_AI_PHS4d[1],
                     Baseline_AI_PHS4h[1],
                     Baseline_AI_PHSInt[1],
                     Baseline_AI_H[1],
                     Baseline_AI_R[1],
                     Baseline_AI_FP[1]
                     ]



#Creating list of Final Contmainations
List_Final_Cont_AI = [x["Final Product Facility"] for x in List_of_Progs_AI]
Col_Names_AI = "BaselineAI PHS4D_AI PHS4H_AI PHSInt_AI HTrad_AI RSTrad_AI FPSTrad_AI".split()
#Creating and Melting the Dataframe of Final Contminations
DF_Final_Cont_AI = pd.concat(List_Final_Cont_AI, axis = 1)
DF_Final_Cont_AI.columns = Col_Names_AI
DF_Final_Cont_AI_melted = DF_Final_Cont_AI.melt()

#Creating Cat plot bar or box, change arguement
H=sns.catplot(x="variable", y="value", kind = "box" ,
            data=DF_Final_Cont_AI_melted)
plt.xlabel("Intervention")
plt.ylabel("Total CFUs at Finished Product")
plt.yscale('log')
plt.title("CFU Final Contamination: Baseline All Interventions")
plt.xticks(rotation=70)

#Means compared to each other. 
Final_Cont_AI_MeanCI= [mean_CI_ONE(x) for x in List_Final_Cont_AI]

for i in range(7):
    print(Final_Cont_AI_MeanCI[i][0]/Final_Cont_AI_MeanCI[0][0])


Reduction_Cont_AI = Calc_red(meanCI =Final_Cont_AI_MeanCI,treatments= 7)

DF_Reduction_Cont_AI = pd.DataFrame({"Treatment": Col_Names_AI,
                             "Reduction": Reduction_Cont_AI} )

DF_Reduction_Cont_AI=DF_Reduction_Cont_AI.sort_values('Reduction',ascending=False).reset_index()

chart = sns.barplot(data = DF_Reduction_Cont_AI, x = "Treatment", y = "Reduction", order=DF_Reduction_Cont_AI['Treatment'])
chart.bar_label(chart.containers[0])
plt.xlabel("Intervention")
plt.ylabel("Percent Reduction from BaselineAI")



#Sampling Results: 

#PHS4d
F_Sampling_Power(Baseline_AI_PHS4d[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_PHS4d[1],"Bef Pre-Harvest Samp")
#PHS4h
F_Sampling_Power(Baseline_AI_PHS4h[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_PHS4h[1],"Bef Pre-Harvest Samp")
#PHSInt    
F_Sampling_Power(Baseline_AI_PHSInt[0],"PH_CFU_Acc","PH_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_PHSInt[1],"Bef Pre-Harvest Samp")
#HS   
F_Sampling_Power(Baseline_AI_H[0],"H_CFU_Acc","H_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_H[1],"Bef Harvest Samp")
#RS
F_Sampling_Power(Baseline_AI_R[0],"R_CFU_Acc","R_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_R[1],"Bef Receiving Samp")
#FPS
F_Sampling_Power(Baseline_AI_FP[0],"FP_CFU_Acc","FP_CFU_Rej")
CFU_Sampling_Stage (Baseline_AI_FP[1],'Bef Final Prod S')

#Final CFUs for all. 
#prints the mean total CFUs and the CIs
Final_Cont_AI_MeanCI

#Getting the Prevalence of Contaminated Packages.

Mean_Quantiles(df = Baseline_AI[2],columnname ="PropCont_A_FP_Whole" ,q1 = 0.05,q2 = 0.95)  

#Proportion of packages contaminated
Intervention_Props_AI = [Baseline_AI[2],
                     Baseline_AI_PHS4d[2],
                     Baseline_AI_PHS4h[2],
                     Baseline_AI_PHSInt[2],
                     Baseline_AI_H[2],
                     Baseline_AI_R[2],
                     Baseline_AI_FP[2]
                     ]

#Creating dataframe of final contamination for every intervention. 
Intervention_PropCont_List_AI = [x["PropCont_A_FP_Whole"] for x in Intervention_Props_AI]

Prop_Cont_AI = pd.concat(Intervention_PropCont_List_AI, axis = 1)
Prop_Cont_AI.columns = Col_Names_AI
Prop_Cont_AI_melted = Prop_Cont_AI.melt()

#Plotting the bar graph or boxplot of the differences. 
H=sns.catplot(x="variable", y="value", kind = "bar" ,
            data=Prop_Cont_AI_melted)
plt.xlabel("Intervention")
plt.ylabel("Proportion Contaminated at Final Product")
plt.xticks(rotation=70)

Mean_Props_AI= [mean_CI_ONE(x) for x in Intervention_PropCont_List_AI]

for i in range(7):
    print(Mean_Props_AI[i][0]/Mean_Props_AI[0][0])

#Histplot
h=sns.displot( data =Prop_Cont_AI_melted, 
            x = "value" , 
            col = "variable", 
            col_wrap=3,
             stat = "count",
             bins = 30,
            facet_kws=dict(sharey=False,sharex= False))
plt.suptitle("Proportion of Paclages Contminated",) 

def specs(x, **kwargs):
    plt.axvline(x.mean(), c='red', ls='-', lw=2.5)
    plt.axvline(x.median(), c='orange', ls='--', lw=2.5)

h.map(specs,"value" )





#%%

sns.barplot(data = Baseline_NI[1])
plt.xticks(rotation=-80)
plt.yscale('log')
plt.yticks([1,10,100,1000,10000,100000])
plt.title("Contamination Progression")
plt.ylabel("CFUs in System")



sns.barplot(data = Baseline_AI[1])
plt.xticks(rotation=-80)
plt.yticks([1,10,100,1000,10000,100000])
plt.yscale('log')
plt.title("Contamination Progression")
plt.ylabel("CFUs in System")



sns.barplot(data = Baseline_AI_PHS4d[1])
plt.xticks(rotation=-80)
plt.yticks([1,10,100,1000,10000,100000])
plt.yscale('log')
plt.title("Contamination Progression")
plt.ylabel("CFUs in System")



dfBaseline_NI_melted= pd.melt(Baseline_NI[1]["Final Product Facility"])
dfBaseline_NI_melted["Type"] = "BaselineNI"



dfBaseline_NI_melted= Baseline_NI[1]["Final Product Facility"]
dfBaseline_NI_melted["Type"] = "BaselineNI"

dfBaseline_NI_PH4D_melted= Baseline_NI_PHS4d[1]["Final Product Facility"]
dfBaseline_NI_PH4D_melted["Type"] = "BaselineNI_PH4d"

NI_compared = pd.concat([dfBaseline_NI_melted, dfBaseline_NI_PH4D_melted], 0)

sns.barplot(data = NI_compared, x  = "variable", y= "value", hue = "Type")
plt.xticks(rotation=-80)
plt.yticks([1,10,100,1000,10000,100000])
plt.yscale('log')
plt.title("Contamination Progression")
plt.ylabel("CFUs in System")


dfBaseline_AI_melted= pd.melt(Baseline_AI[1])
dfBaseline_AI_melted["Type"] = "BaselineAI"




dfBaseline_AI_PH4D_melted= pd.melt(Baseline_AI_PHS4d[1])
dfBaseline_AI_PH4D_melted["Type"] = "BaselineAI_PH4d"

AI_compared = pd.concat([dfBaseline_AI_melted, dfBaseline_AI_PH4D_melted], 0)

sns.barplot(data = AI_compared, x  = "variable", y= "value", hue = "Type")
plt.xticks(rotation=-80)
plt.yticks([1,10,100,1000,10000,100000])
plt.yscale('log')
plt.title("Contamination Progression")
plt.ylabel("CFUs in System")


df_comparisons = pd.DataFrame ({"Baseline NI": Baseline_NI[1]["Final Product Facility"],
                              "PH4d": Baseline_NI_PHS4d[1]["Final Product Facility"] })
df_comparisons_melt = pd.melt(df_comparisons)


sns.barplot(data =df_comparisons_melt, x  = "variable", y= "value")
plt.xticks(rotation=-80)
#plt.yticks([1,10,100,1000,10000,100000])
#plt.yscale('log')
plt.title("Contamination at Final Product")
plt.ylabel("CFUs ")


df_comparisons_AI = pd.DataFrame ({"Baseline AI": Baseline_AI[1]["Final Product Facility"],
                              "PH4d": Baseline_AI_PHS4d[1]["Final Product Facility"] })
df_comparisons_melt_AI = pd.melt(df_comparisons_AI)


sns.barplot(data =df_comparisons_melt_AI, x  = "variable", y= "value")
plt.xticks(rotation=-80)
#plt.yticks([1,10,100,1000,10000,100000])
#plt.yscale('log')
plt.title("Contamination at Final Product")
plt.ylabel("CFUs ")