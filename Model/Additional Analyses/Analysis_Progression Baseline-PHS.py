# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 08:21:31 2021

@author: Gustavo Reyes
"""

#Analysis Compating Sampling Step vs Ending Step

#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
#sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

# %%
from importlib import reload
import Listz
import pandas as pd
import MainModel3z
import Inputz
import SCInputz
import ContCondz
import ScenCondz
import InFunz
import OutFunz
import ContScen
import Funz
from matplotlib import pyplot as plt
import seaborn as sns
import sys
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append(
    'C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

# %% Libraries

# Own Libraries

# Main Model Loops

# %%
# Progression DataFrame.
#Progression_DFS = []
# %% BAseline Sampling:
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

# %% Baseline Sampling
reload(ScenCondz)
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = True  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = False
ScenCondz.H_Sampling = False
ScenCondz.R_Sampling = False
ScenCondz.FP_Sampling = False
ScenCondz.Field_Pack = False

reload(SCInputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists
Main_Mod_Outs = MainModel3z.F_MainLoop()
D_Baseline = Main_Mod_Outs[1]
DProg_Baseline = Main_Mod_Outs[0]

# %% Pre-Harvest
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

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_PH4d = Main_Mod_Outs[1]
DProg_PH4d = Main_Mod_Outs[0]

 #%%
def melting_type(df, typename):
    df_melt = pd.melt(df)
    df_melt["type"] = typename
    return df_melt

#PRogression
DProg_Baseline_melt = melting_type(DProg_Baseline, "Baseline")
DProg_PH4d_melt = melting_type(DProg_PH4d, "PH 4d")

all_list = [DProg_Baseline_melt, DProg_PH4d_melt]
Frame_all = pd.concat(all_list)

sns.catplot(x="variable", y="value", hue="type", kind="box",
            data=Frame_all, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression Through System")
plt.xticks(rotation=70)

#%% 

Final_Base = DProg_Baseline["Final Product Facility"]
Final_PH4d = DProg_PH4d["Final Product Facility"]

Final_Compared=pd.concat([Final_Base, Final_PH4d], axis=1, ignore_index=True)
Final_Compared.rename(columns={0: 'Final Baseline', 1: 'Final PH4d'}, inplace=True)
Final_Compared_melt=pd.melt(Final_Compared)

H=sns.catplot(x="variable", y="value", kind="box",
            data=Final_Compared_melt)
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")
plt.title("CFUs in Final Product")
#plt.xticks(rotation=70)


medians = Final_Compared_melt.groupby(['variable'])['value'].median()



#%%

Frame_beg = Frame_all.loc[(Frame_all["variable"] == "Bef Pre-Harvest Samp") | (Frame_all["variable"] == "Aft Pre-Harvest Samp")]
sns.catplot(x="variable", y="value", hue="type", kind="box",
            data=Frame_beg)
plt.xlabel("Sampling Scenario")
plt.ylabel("Pre-Harvest Sampling Step")
plt.title("Contamination Progression Through System")
plt.xticks(rotation=70)

medians = Frame_beg.groupby(['variable', 'type'])['value'].median()
