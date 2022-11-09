# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:37:03 2021

@author: gareyes3
"""

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
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append(
    'C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')
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
Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]
Tuning_Contamination_levels = [CFU_10000g,CFU_1000g,CFU_100g,CFU_10g,CFU_g,CFU_0_1g]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for k in Tuning_Contamination_levels:
    for i in Tuning_SampleSize:
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
            
            SCInputz.PSHazard_lvl = k
            SCInputz.sample_size_PH = i
            #SCInputz.BGHazard_lvl = j
            SCInputz.No_Grabs_PH = j
            
            Main_Mod_Outs = MainModel3z.F_MainLoop()
            OutputDF = Main_Mod_Outs[1]
            ProgDF = Main_Mod_Outs[0]
            
            DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
            #DF= pd.melt(DF)
            DF["ContLevel"] = k
            DF["GrabNo"] = j
            DF["SampSize"] = i
            Output_Collection_List.append(DF)
             

Combined_df_PS = pd.concat(Output_Collection_List)
#%%

Combined_df_PS["PH_CFU_PerA"] = 1-Combined_df_PS["PH_CFU_PerR"]
Combined_df_PS["ContLevel"]=np.log10(Combined_df_PS["ContLevel"])


g = sns.FacetGrid(Combined_df_PS, col="SampSize", hue = "GrabNo",col_wrap=3, height = 6)
g.map(sns.lineplot, "ContLevel","PH_CFU_PerA")
g.add_legend()
g= (g.set_axis_labels(x_var ="Initial Cont Level Log CFU",y_var = "% of CFU Acc (1-POWER)" ))
g.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
g.fig.suptitle('Point Source Contamination -Effect of  #Grabs | Sample Mass | Initial Contamination Levels')




sns.set_style("whitegrid")
sns.set(font_scale = 1.4)


Combined_df_PS["ContLevel"].replace({CFU_10000g: "1 CFU/10kg",
                                     CFU_1000g: "1 CFU/kg",
                                     CFU_100g: "1 CFU/100g", 
                                     CFU_10g: "1 CFU/10"}, inplace=True)



H = sns.catplot(x="GrabNo", y="PH_CFU_PerR" ,
                col="ContLevel", col_wrap=3,hue="SampSize",data=Combined_df_PS, kind = "point")
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "Proportion of CFU Rejected (POWER)" )) 
H.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
H.fig.suptitle('Point Source-Contamination: Composite Mass vs. No Grabs vs. Contamination Level')


sns.regplot(x="GrabNo", y="PH_CFU_PerR" ,)


#%% SYSTEMATIC

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
#%%
Tuning_SampleSize = [60,120,300,600,1200] #list(range(100,1200,200))
#Tuning_HazardLevel = [4536,45360,450360,4503600]#list(range(10000,130000,10000)) #1 CFU/10kg, 1CFU kg, 1CFU/100g, 1CFU/10g. 
Tuning_Grab_number = [1,30,60,120,320]
Tuning_Contamination_levels = [1000,10000,100000,1000000]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"""]
Output_Collection_List = [] #First Index
for k in Tuning_Contamination_levels:
    for i in Tuning_SampleSize:
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
            
            SCInputz.PSHazard_lvl = k
            SCInputz.sample_size_PH = i
            #SCInputz.BGHazard_lvl = j
            SCInputz.No_Grabs_PH = j
            
            Main_Mod_Outs = MainModel3z.F_MainLoop()
            OutputDF = Main_Mod_Outs[1]
            ProgDF = Main_Mod_Outs[0]
            
            DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
            #DF= pd.melt(DF)
            DF["ContLevel"] = k
            DF["GrabNo"] = j
            DF["SampSize"] = i
            Output_Collection_List.append(DF)
             

Combined_df_BG = pd.concat(Output_Collection_List)



H = sns.catplot(x="GrabNo", y="PH_CFU_PerR" ,
                col="ContLevel", col_wrap=3,hue="SampSize",data=Combined_df_BG, kind = "point" )
H= (H.set_axis_labels(x_var ="Number  of Grabs",y_var = "Proportion of CFU Rejected (POWER)" )) 
H.fig.subplots_adjust(top=0.9) # adjust the Figure in rp
H.fig.suptitle('Bakcground-Contamination: Composite Mass vs. No Grabs vs. Contamination Level')


