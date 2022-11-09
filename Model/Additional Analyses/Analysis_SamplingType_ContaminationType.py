# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 14:10:19 2021

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


# %% Pre-Harvest- Sampling Type and Contamination Type
#%%
#Tuning Parameter SampleSize and Hazard Level


#Only 2 can be tuned simulatenously
Tuning_Sampling_Type = ["4d", "4h", "Int"] #Tune sampling type
Tuning_Sampling_Step = ["PH"] #Tune Sampling Step
Tuning_Contamination_Type = ["BackGround", "PointSource", "Systematic"]

Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"]
Output_Collection_List = [] #First Index
for i in Tuning_Contamination_Type:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Sampling_Type:  
        
        #SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        #SCInputz.No_Grabs_PH = j
        Sampling_Step = "PH" #Sampling Step we are evaluating. 
        Sampling_Type =j
        Contamination_Type = i
        
        # Sampling Conditions, Baseline all conditions are off
        reload(ScenCondz) #Refreshing ScenCondz. 
        if Sampling_Step == "Baseline":
            ScenCondz.Baseline_Sampling = 1
        elif Sampling_Step == "PH":
            ScenCondz.PH_Sampling = 1
            if Sampling_Type == "4d":
                ScenCondz.PHS_4d = 1
            if Sampling_Type == "4h":
                ScenCondz.PHS_4h = 1
            if Sampling_Type == "Int":
                ScenCondz.PHS_Int = 1
        elif Sampling_Step == "H":
            ScenCondz.H_Sampling = 1
            if Sampling_Type == "Trad":
                ScenCondz.HS_Trad = 1
            if Sampling_Type == "Agg":
                ScenCondz.HS_Agg = 1
        elif Sampling_Step == "Receiving":
            ScenCondz.R_Sampling = 1
        elif Sampling_Step == "FP":
            ScenCondz.FP_Sampling = 1
            if Sampling_Type == "Trad":
                ScenCondz.FP_Trad = 1
            if Sampling_Type == "Agg":
                ScenCondz.FP_Agg = 1
                
        reload(ContCondz)        
        # Contamination Challenges
        if Contamination_Type=="BackGround":  
            ContCondz.Background_C = True
        if Contamination_Type=="PointSource":
            ContCondz.Point_Source_C = True
        if Contamination_Type=="Systematic":
            ContCondz.Systematic_C = True
        

        reload(SCInputz)  # Reload Inputz
        
        #REJECTION RULE AREA-----------------------------------------------------------------------------
        #Updating the Rejection Rule:
        if Sampling_Type == "Int":
            SCInputz.RR_PH_Int = "Sublot"
        
        #Reloading the list
        reload(Listz)  # ReUPdate Lists
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["HLev"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)

#%%
Main_list = []     
for  i in Output_Collection_List  :   
    Combined_df = pd.concat(i)
    Main_list.append(Combined_df)


Main_Combined_SS = pd.concat(Main_list)
#Main_Combined_PRej = Main_Combined.loc[Main_Combined['variable'] == "PH_CFU_PerR"]

Main_Combined_SS["Efficacy"] = Main_Combined_SS["PH_CFU_PerR"] -  Main_Combined_SS["PH_Wei_PerR"]


means_efficacy_ss = Main_Combined_SS.groupby(by=["HLev", "SampSize"])['Efficacy','PH_CFU_PerR','PH_Wei_PerR'].mean()

#Facet Line Plot
'''
g= sns.lmplot(y="PH_CFU_PerR", x="PH_Wei_PerR", hue="HLev",
              col = "SampSize", 
              data=Main_Combined, scatter=True,
              line_kws=dict(alpha=0.8),
              scatter_kws=dict(alpha = 0.3)).set_titles('{col_name}')
g= (g.set_axis_labels(x_var ="Percent Weight Rejected at Sampling Step",y_var = "Percent of CFU Rejected at Sampling Step" ,))
g= (g._legend.set_title("Sampling Type"))

#Lineplot
myPlot = sns.FacetGrid(col="SampSize", hue="HLev", data=Main_Combined_SS,size=5)
myPlot = myPlot.map(sns.regplot, "PH_Wei_PerR", "PH_CFU_PerR" ,line_kws=dict(alpha=1), scatter_kws=dict(alpha = 0.3))
myPlot = myPlot.map_dataframe(plt.plot, [1,0], [1,0], 'r-',color = "black",linestyle='dashed',linewidth=1, alpha = 0.3).add_legend().set_axis_labels("PH_Wei_PerR", "PH_CFU_PerR")
myPlot= (myPlot.set_axis_labels(x_var ="Percent Weight Rejected at Sampling Step",y_var = "Percent of CFU Rejected at Sampling Step" ,))
myPlot= (myPlot._legend.set_title("Sampling Type"))

plt.show()
'''

#Efficacy Violin Plot
H = sns.catplot(x="HLev", y="Efficacy",
                col="SampSize", data=Main_Combined_SS, kind = "bar")
H= (H.set_axis_labels(x_var ="Sampling Type",y_var = "Efficacy" )
    .set_xticklabels(["4 Day", "4 Hour", "Intense"]))
#H.map(plt.axhline, y=1, ls='--', c='red',)
H.fig.subplots_adjust(top=0.8) # adjust the Figure in rp
H.fig.suptitle('Contamination | Weight Rejected Efficacy: Sublot Rejection Rule ')

#Percentage Rejected Plot
H = sns.catplot(x="HLev", y="PH_CFU_PerR",
                col="SampSize", data=Main_Combined_SS, kind = "bar")
H= (H.set_axis_labels(x_var ="Sampling Type",y_var = "Proportion Rejected" )
    .set_xticklabels(["4 Day", "4 Hour", "Intense"]))
#H.map(plt.axhline, y=1, ls='--', c='red',)
H.fig.subplots_adjust(top=0.8) # adjust the Figure in rp
H.fig.suptitle('Proportion Rejected by Sampling Type')



#Histogrm of percentage rejected
myPlot2 = sns.FacetGrid(col="SampSize", hue="HLev", data=Main_Combined_SS,size=5)
myPlot2 = myPlot2.map(sns.histplot, "PH_CFU_PerR", binwidth=.1, alpha=0.4)



# %% Harvest- Sampling Type and Contamination Type


#Only 2 can be tuned simulatenously
Tuning_Sampling_Type = ["Trad", "Agg"] #Tune sampling type
Tuning_Sampling_Step = ["H"] #Tune Sampling Step
Tuning_Contamination_Type = ["BackGround", "PointSource", "Systematic"]



Desired_Outputs = ["PH_CFU_PerR", "PH_Wei_PerR"]
Output_Collection_List = [] #First Index
for i in Tuning_Contamination_Type:
    Hazard_Level_list = [] #SEcond INdex
    for j in Tuning_Sampling_Type:  
        
        #SCInputz.sample_size_PH = i
        #SCInputz.BGHazard_lvl = j
        #SCInputz.No_Grabs_PH = j
        Sampling_Step = "PH" #Sampling Step we are evaluating. 
        Sampling_Type =j
        Contamination_Type = i
        
        # Sampling Conditions, Baseline all conditions are off
        reload(ScenCondz) #Refreshing ScenCondz. 
        if Sampling_Step == "Baseline":
            ScenCondz.Baseline_Sampling = 1
        elif Sampling_Step == "PH":
            ScenCondz.PH_Sampling = 1
            if Sampling_Type == "4d":
                ScenCondz.PHS_4d = 1
            if Sampling_Type == "4h":
                ScenCondz.PHS_4h = 1
            if Sampling_Type == "Int":
                ScenCondz.PHS_Int = 1
        elif Sampling_Step == "H":
            ScenCondz.H_Sampling = 1
            if Sampling_Type == "Trad":
                ScenCondz.HS_Trad = 1
            if Sampling_Type == "Agg":
                ScenCondz.HS_Agg = 1
        elif Sampling_Step == "Receiving":
            ScenCondz.R_Sampling = 1
        elif Sampling_Step == "FP":
            ScenCondz.FP_Sampling = 1
            if Sampling_Type == "Trad":
                ScenCondz.FP_Trad = 1
            if Sampling_Type == "Agg":
                ScenCondz.FP_Agg = 1
                
        reload(ContCondz)        
        # Contamination Challenges
        if Contamination_Type=="BackGround":  
            ContCondz.Background_C = True
        if Contamination_Type=="PointSource":
            ContCondz.Point_Source_C = True
        if Contamination_Type=="Systematic":
            ContCondz.Systematic_C = True
        

        reload(SCInputz)  # Reload Inputz
        
        #REJECTION RULE AREA-----------------------------------------------------------------------------
        #Updating the Rejection Rule:
        if Sampling_Type == "Int":
            SCInputz.RR_PH_Int = "Sublot"
        
        
        #Reloading the list
        reload(Listz)  # ReUPdate Lists
        
        Main_Mod_Outs = MainModel3z.F_MainLoop()
        OutputDF = Main_Mod_Outs[1]
        ProgDF = Main_Mod_Outs[0]
        
        DF= OutFunz.F_Output_get_cols(Outdf = OutputDF , ColNames = Desired_Outputs)
        #DF= pd.melt(DF)
        DF["HLev"] = j
        DF["SampSize"] = i
        Hazard_Level_list.append(DF)
    Output_Collection_List.append(Hazard_Level_list)
    
means_efficacy_ss = Main_Combined_SS.groupby(by=["HLev", "SampSize"])['Efficacy','PH_CFU_PerR','PH_Wei_PerR'].mean()

means_CFU_R_ss = Combined_all.groupby(by=['Sample Mass', "GrabNo","SampType"])['PH_CFU_PerR'].mean()





