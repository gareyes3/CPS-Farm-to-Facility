# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 10:21:20 2021

@author: gareyes3
"""

#Background contamination in our case won't be necessary. Efficacy is always 0. 
#%%
import sys, os
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

# %%
from importlib import reload
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
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append(
    'C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')


#%% Systematic Contamination
#---10% Clusters. 
#--- Initial Hazard level 50k CFU 
# %% Changing to Point Source Contamination
# Contamination Challenges
ContCondz.Background_C = True
ContCondz.Point_Source_C = False
ContCondz.Systematic_C = False

# Harvester Contamination
ContCondz.Crew_C = False
ContCondz.Harvester_C = False

# Processing equipment
ContCondz.PE_C = False
ContCondz.PE_Cont_Loc = False,  # 1,2,3,4,5
# 1 = Shredder, #2 = Belt, #3 = Washing, #4 Shaker, #5Centrifuge
ContCondz.Pack_C = False
#%% 4Days--------------------------------------------------

Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for i in Tuning_SampleSize:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Grab_number:  
        # Sampling Conditions, Baseline all conditions are off
        ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
        ScenCondz.PH_Sampling = 1
        ScenCondz.H_Sampling = 0
        ScenCondz.R_Sampling = 0
        ScenCondz.FP_Sampling = 0
        # Pre_Harvest 4 Days
        # Pre-Harvest
        ScenCondz.PHS_4d = 1  # Scenario 1
        ScenCondz.PHS_4h = 0  # Scenario 2
        ScenCondz.PHS_Int = 0  # Scenario 3
        
        reload(SCInputz)  # Reload Inputz
        reload(Listz)  # ReUPdate Lists
        
        
        SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        SCInputz.No_Grabs_PH = j
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["GrabNo"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)
#%%
Main_list = []     
for  i in Output_Collection_List  :   
    Combined_df = pd.concat(i)
    Main_list.append(Combined_df)

Main_Combined = pd.concat(Main_list)


Main_Combined["Efficacy"] = Main_Combined["PH_CFU_PerR"] -  Main_Combined["PH_Wei_PerR"]
Main_Combined["SampType"] = "PH4d"

#%%  4 HOURS--------------------------------------------------------------------

Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for i in Tuning_SampleSize:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Grab_number:  
        # Sampling Conditions, Baseline all conditions are off
        ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
        ScenCondz.PH_Sampling = 1
        ScenCondz.H_Sampling = 0
        ScenCondz.R_Sampling = 0
        ScenCondz.FP_Sampling = 0
        # Pre_Harvest 4 Days
        # Pre-Harvest
        ScenCondz.PHS_4d = 0  # Scenario 1
        ScenCondz.PHS_4h = 1  # Scenario 2
        ScenCondz.PHS_Int = 0  # Scenario 3
        
        reload(SCInputz)  # Reload Inputz
        reload(Listz)  # ReUPdate Lists
        
        
        SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        SCInputz.No_Grabs_PH = j
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["GrabNo"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)

#%%
Main_list = []     
for  i in Output_Collection_List  :   
    Combined_df = pd.concat(i)
    Main_list.append(Combined_df)

Main_Combined_4h = pd.concat(Main_list)


Main_Combined_4h["Efficacy"] = Main_Combined_4h["PH_CFU_PerR"] -  Main_Combined_4h["PH_Wei_PerR"]
Main_Combined_4h["SampType"] = "PH4h"



#%% Intense Sampling
Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for i in Tuning_SampleSize:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Grab_number:  
        # Sampling Conditions, Baseline all conditions are off
        ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
        ScenCondz.PH_Sampling = 1
        ScenCondz.H_Sampling = 0
        ScenCondz.R_Sampling = 0
        ScenCondz.FP_Sampling = 0
        # Pre_Harvest 4 Days
        # Pre-Harvest
        ScenCondz.PHS_4d = 0  # Scenario 1
        ScenCondz.PHS_4h = 0  # Scenario 2
        ScenCondz.PHS_Int = 1 # Scenario 3
        
        
        reload(SCInputz)  # Reload Inputz
        reload(Listz)  # ReUPdate Lists
        
        
        if ScenCondz.PHS_Int == 1:
            SCInputz.RR_PH_Int = "Sublot" 
        
        SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        SCInputz.No_Grabs_PH = j
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["GrabNo"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)

#%%
Main_list = []     
for  i in Output_Collection_List  :   
    Combined_df = pd.concat(i)
    Main_list.append(Combined_df)

Main_Combined_Int = pd.concat(Main_list)


Main_Combined_Int["Efficacy"] = Main_Combined_Int["PH_CFU_PerR"] -  Main_Combined_Int["PH_Wei_PerR"]
Main_Combined_Int["SampType"] = "Intense"


#%% HArvest Traditional
Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for i in Tuning_SampleSize:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Grab_number:  
        # Sampling Conditions, Baseline all conditions are off
        ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
        ScenCondz.PH_Sampling = 1
        ScenCondz.H_Sampling = 0
        ScenCondz.R_Sampling = 0
        ScenCondz.FP_Sampling = 0
        # Pre_Harvest 4 Days
        # Pre-Harvest
        ScenCondz.PHS_4d = 0  # Scenario 1
        ScenCondz.PHS_4h = 0  # Scenario 2
        ScenCondz.PHS_Int = 1 # Scenario 3
        
        
        reload(SCInputz)  # Reload Inputz
        reload(Listz)  # ReUPdate Lists
        
        
        if ScenCondz.PHS_Int == 1:
            SCInputz.RR_PH_Int = "Sublot" 
        
        SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        SCInputz.No_Grabs_PH = j
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["GrabNo"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)

#%%
Main_list = []     
for  i in Output_Collection_List  :   
    Combined_df = pd.concat(i)
    Main_list.append(Combined_df)

Main_Combined_HTrad = pd.concat(Main_list)


Main_Combined_HTrad["Efficacy"] = Main_Combined_HTrad["PH_CFU_PerR"] -  Main_Combined_HTrad["PH_Wei_PerR"]
Main_Combined_HTrad["SampType"] = "H Trad"






#%%
Frames = [Main_Combined,Main_Combined_4h,Main_Combined_Int,Main_Combined_HTrad]
Combined_all = pd.concat(Frames)
#%%
H = sns.catplot(x="GrabNo", y="PH_CFU_PerR" ,
                hue="Sample Mass", col="SampType",data=Combined_all, kind = "bar")
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "Per CFU Rejected (POWER)" ))
H.fig.subplots_adjust(top=0.8) # adjust the Figure in rp
H.fig.suptitle('Uniform Contamination, Total 50k CFUs')



#Efficacy Violin Plot
sns.set(font_scale = 1.4)
Combined_all = Combined_all.rename(columns={"SampSize":"Sample Mass"})
H = sns.catplot(x="GrabNo", y="Efficacy" ,
                col="Sample Mass", col_wrap=3,hue="SampType",data=Combined_all, kind = "bar")
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "Efficacy" ))
#plt.xlabel("Number  of Grabs")
#plt.ylabel("Efficacy")
#H.map(plt.axhline, y=1, ls='--', c='red',)
#H.fig.subplots_adjust(top=0.8) # adjust the Figure in rp
#H.fig.suptitle('Point Source-Contamination, Efficacy Comparison Different Sampling Methods')

sns.set(font_scale = 1.4)
Combined_all = Combined_all.rename(columns={"SampSize":"Sample Mass"})
H = sns.catplot(x="GrabNo", y="PH_CFU_PerR" ,
                col="Sample Mass", col_wrap=3,hue="SampType",data=Combined_all, kind = "bar")
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "Prportion Rejected" ))



sns.set(font_scale = 1.4)
Combined_all = Combined_all.rename(columns={"SampSize":"Sample Mass"})
H = sns.catplot(x="GrabNo", y="PH_CFU_PerR" ,
                col="Sample Mass", row="SampType",data=Combined_all, kind = "bar")
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "PH_CFU_PerR" ))
H.xlabel("Number  of Grabs")
H.ylabel("Percentage of CFU Rejected")
H.map(plt.axhline, y=1, ls='--', c='red',)
H.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
H.fig.suptitle('Point Source-Contamination, Efficacy Comparison Different Sampling Methods')







means_efficacy_all = Combined_all.groupby(by=['Sample Mass', "GrabNo","SampType"])['Efficacy','PH_CFU_PerR', 'PH_Wei_PerR'].mean()

means_CFU_R_all = Combined_all.groupby(by=['Sample Mass', "GrabNo","SampType"])['PH_CFU_PerR'].mean()



#%% SYSTEMATIC CONTAMINATION



