# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 09:26:18 2021

@author: gareyes3
"""

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
ContCondz.Systematic_C = True

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

# Pre_Harvest 4 hr
# Pre-Harvest

ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 1  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_PH4h = Main_Mod_Outs[1]
DProg_PH4h = Main_Mod_Outs[0]

# Pre_Harvest Int
# Pre-Harvest


ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 1  # Scenario 3

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists
reload(SCInputz)  # ReUPdate Lists
Main_Mod_Outs = MainModel3z.F_MainLoop()
D_PHInt = Main_Mod_Outs[1]
DProg_PHInt = Main_Mod_Outs[0]


# %% Harvest
reload(ScenCondz)
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 1
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0

# Harvest Traditional


ScenCondz.HS_Trad = 1  # Scenario 5
ScenCondz.HS_Agg = 0  # Scenario 6

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_HTrad = Main_Mod_Outs[1]
DProg_HTrad = Main_Mod_Outs[0]

# Harvest Aggregative
ScenCondz.HS_Trad = 0  # Scenario 5
ScenCondz.HS_Agg = 1  # Scenario 6

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists
reload(SCInputz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_HAgg = Main_Mod_Outs[1]
DProg_HAgg = Main_Mod_Outs[0]

# %% Receiving
reload(ScenCondz)
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 1
ScenCondz.FP_Sampling = 0

# Harvest Traditional
reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_R = Main_Mod_Outs[1]
DProg_R = Main_Mod_Outs[0]

# %% Final Product
reload(ScenCondz)
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 1

# Harvest Traditional

ScenCondz.FPS_Trad = 1  # Scenario 5
ScenCondz.FPS_Agg = 0  # Scenario 6

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_FPSTrad = Main_Mod_Outs[1]
DProg_FPSTrad = Main_Mod_Outs[0]

# Harvest Aggregative

ScenCondz.FPS_Trad = 0  # Scenario 5
ScenCondz.FPS_Agg = 1  # Scenario 6

reload(Inputz)  # Reload Inputz
reload(Listz)  # ReUPdate Lists

Main_Mod_Outs = MainModel3z.F_MainLoop()
D_FPSAgg = Main_Mod_Outs[1]
DProg_FPSAgg = Main_Mod_Outs[0]

#%%

##CFU that leak through the system
L_Final_Conts = {
    'Baseline': D_Baseline["FP_CFU_Acc"],
    'PH4d': D_PH4d["FP_CFU_Acc"],
    'PH4h': D_PH4h["FP_CFU_Acc"],
    'PHInt': D_PHInt["FP_CFU_Acc"],
    'HTrad': D_HTrad["FP_CFU_Acc"],
    #HAgg': D_HAgg["FP_CFU_Acc"],
    'R': D_R["FP_CFU_Acc"],
    'FPTrad': D_FPSTrad["FP_CFU_Acc"],
    #'FPAgg': D_FPSAgg["FP_CFU_Acc"],
}

df_L_Final_Conts = pd.DataFrame(L_Final_Conts)
df_L_Final_Conts_melted = pd.melt(df_L_Final_Conts)

Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_Conts_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Total CFU in Final Product")
plt.title("Total CFU that leak through system")

#Sampling POWER
L_Final_Conts = {
    'Baseline': D_Baseline["PH_CFU_PerR"],
    'PH4d': D_PH4d["PH_CFU_PerR"],
    'PH4h': D_PH4h["PH_CFU_PerR"],
    'PHInt': D_PHInt["PH_CFU_PerR"],
    'HTrad': D_HTrad["H_CFU_PerR"],
    #HAgg': D_HAgg["FP_CFU_Acc"],
    'R': D_R["R_CFU_PerR"],
    'FPTrad': D_FPSTrad["FP_CFU_PerR"],
    #'FPAgg': D_FPSAgg["FP_CFU_Acc"],
}

df_L_Final_Conts = pd.DataFrame(L_Final_Conts)
df_L_Final_Conts_melted = pd.melt(df_L_Final_Conts)

Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_Conts_melted)
plt.xlabel("Sampling Type")
plt.ylabel("% of CFU Rejected by Sampling Plan")
plt.title("Power of Sampling Plans")


def melting_type(df, typename):
    df_melt = pd.melt(df)
    df_melt["type"] = typename
    return df_melt

#PRogression
DProg_Baseline_melt = melting_type(DProg_Baseline, "Baseline")
DProg_PH4d_melt = melting_type(DProg_PH4d, "PH 4d")
DProg_PH4h_melt = melting_type(DProg_PH4h, "PH 4h")
DProg_PHInt_melt = melting_type(DProg_PHInt, "PH Int")
DProg_HTrad_melt = melting_type(DProg_HTrad, "H Trad")
DProg_R_melt = melting_type(DProg_R, "Receiving")
DProg_FPSTrad_melt = melting_type(DProg_FPSTrad, "FPSTrad")

all_list = [DProg_Baseline_melt, DProg_PH4d_melt, DProg_PH4h_melt, DProg_PHInt_melt,
            DProg_HTrad_melt, DProg_R_melt, DProg_FPSTrad_melt]
Frame_all = pd.concat(all_list)

sns.catplot(x="variable", y="value", hue="type", kind="bar",
            data=Frame_all, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression Through System")
plt.xticks(rotation=70)




# %%Final graphs

L_Final_Conts = {
    'Baseline': D_Baseline["Total_CFU_A"],
    'PH4d': D_PH4d["Total_CFU_A"],
    'PH4h': D_PH4h["Total_CFU_A"],
    'PHInt': D_PHInt["Total_CFU_A"],
    'HTrad': D_HTrad["Total_CFU_A"],
    'HAgg': D_HAgg["Total_CFU_A"],
    'R': D_R["Total_CFU_A"],
    'FPTrad': D_FPSTrad["Total_CFU_A"],
    'FPAgg': D_FPSAgg["Total_CFU_A"],
}

df_L_Final_Conts = pd.DataFrame(L_Final_Conts)
df_L_Final_Conts_melted = pd.melt(df_L_Final_Conts)


#Scenariosplot =sns.catplot(x="variable", y="value", data=df_L_Final_Conts_melted ,kind="bar",  height=4, aspect=2 )
Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_Conts_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Total CFU in Final Product")


L_Final_Per = {
    'PH4d': D_PH4d["PerRejected at PH"],
    'PH4h': D_PH4h["PerRejected at PH"],
    'PHInt': D_PHInt["PerRejected at PH"],
    'HTrad': D_HTrad["PerRejected at H"],
    'HAgg': D_HAgg["PerRejected at H"],
    'R': D_R["PerRejected at R"],
    'FPTrad': D_FPSTrad["PerRejected at FP"],
    'FPAgg': D_FPSAgg["PerRejected at FP"],
}

df_L_Final_Per = pd.DataFrame(L_Final_Per)
df_L_Final_Per_melted = pd.melt(df_L_Final_Per)


#Scenariosplot =sns.catplot(x="variable", y="value", data=df_L_Final_Per_melted ,kind="bar",  height=4, aspect=2 )
Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_Per_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Percentage CFU Rejected at Sampling Step")


# Weight Rejected, Weight Accepted.
L_Final_PerRejWeight = {
    'PH4d': D_PH4d["PerRejectedWeight"],
    'PH4h': D_PH4h["PerRejectedWeight"],
    'PHInt': D_PHInt["PerRejectedWeight"],
    'HTrad': D_HTrad["PerRejectedWeight"],
    'HAgg': D_HAgg["PerRejectedWeight"],
    'R': D_R["PerRejectedWeight"],
    'FPTrad': D_FPSTrad["PerRejectedWeight"],
    'FPAgg': D_FPSAgg["PerRejectedWeight"],
}

df_L_Final_PerRejWeight = pd.DataFrame(L_Final_PerRejWeight)
df_L_Final_PerRejWeight_melted = pd.melt(df_L_Final_PerRejWeight)


#Scenariosplot =sns.catplot(x="variable", y="value", data=df_L_Final_Per_melted ,kind="bar",  height=4, aspect=2 )
Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_PerRejWeight_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Percentage total lb. Rejected at Sampling Step")

# TotalWeightRejected
# Weight Rejected, Weight Accepted.
L_Final_RejWeight = {
    'PH4d': D_PH4d["Total_Weight_R"],
    'PH4h': D_PH4h["Total_Weight_R"],
    'PHInt': D_PHInt["Total_Weight_R"],
    'HTrad': D_HTrad["Total_Weight_R"],
    'HAgg': D_HAgg["Total_Weight_R"],
    'R': D_R["Total_Weight_R"],
    'FPTrad': D_FPSTrad["Total_Weight_R"],
    'FPAgg': D_FPSAgg["Total_Weight_R"],
}

df_L_Final_RejWeight = pd.DataFrame(L_Final_RejWeight)
df_L_Final_RejWeight_melted = pd.melt(df_L_Final_RejWeight)


#Scenariosplot =sns.catplot(x="variable", y="value", data=df_L_Final_Per_melted ,kind="bar",  height=4, aspect=2 )
Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_RejWeight_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Total lb. Rejected due to sampling")


# Weight Rejected, Weight Accepted Compariason
L_Final_PerRatio = {
    'PH4d': D_PH4d["PerRejectedWeight"] / D_PH4d["PerRejected at PH"],
    'PH4h': D_PH4h["PerRejectedWeight"] / D_PH4h["PerRejected at PH"],
    'PHInt': D_PHInt["PerRejectedWeight"] / D_PHInt["PerRejected at PH"],
    'HTrad': D_HTrad["PerRejectedWeight"] / D_HTrad["PerRejected at H"],
    'HAgg': D_HAgg["PerRejectedWeight"] / D_HAgg["PerRejected at H"],
    'R': D_R["PerRejectedWeight"] / D_R["PerRejected at R"],
    'FPTrad': D_FPSTrad["PerRejectedWeight"] / D_FPSTrad["PerRejected at FP"],
    'FPAgg': D_FPSAgg["PerRejectedWeight"] / D_FPSAgg["PerRejected at FP"],
}

df_L_Final_PerRatio = pd.DataFrame(L_Final_PerRatio)
df_L_Final_PerRatio_melted = pd.melt(df_L_Final_PerRatio)


#Scenariosplot =sns.catplot(x="variable", y="value", data=df_L_Final_Per_melted ,kind="bar",  height=4, aspect=2 )
Scenariosplot = sns.boxplot(x="variable", y="value",
                            data=df_L_Final_PerRatio_melted)
plt.xlabel("Sampling Type")
plt.ylabel("Percentage lb Rejected/  Percentage CFU Rejected")


# Combined Data PEr Rejected and per Lb Rejected.
Percentage_Rejected = df_L_Final_Per_melted
Percentage_Rejected['Type'] = "Proportion of CFU Rejected"

Percentage_lb_Rejected = df_L_Final_PerRejWeight_melted
Percentage_lb_Rejected['Type'] = "Proportion of lb Rejected"

Combined_Percentages = pd.concat([Percentage_Rejected,
                                  Percentage_lb_Rejected])
# Proportion of CFU rejected Paired Proportion of LB Rejected
sns.set(rc={'figure.figsize': (8, 6)})
sns.set(style="darkgrid")
Scenariosplot = sns.lineplot(
    x="variable", y="value", hue="Type", data=Combined_Percentages)
plt.xlabel("Sampling Scenario")
plt.ylabel("Proportion Rejected")
plt.title("Testing Plan Efficacy")
Scenariosplot.legend_.set_title('Type')
# Put the legend out of the figure
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Boxplot
Scenariosplot = sns.boxplot(x="variable", y="value",
                            hue="Type", data=Combined_Percentages)
plt.xlabel("Sampling Scenario")
plt.ylabel("Proportion Rejection")
Scenariosplot.legend_.set_title('Type')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


Combined_Percentages_Dif = Percentage_Rejected["value"] - \
    Percentage_lb_Rejected["value"]
Combined_Percentages_Dif = pd.DataFrame(Combined_Percentages_Dif)
Combined_Percentages_Dif["Type"] = Percentage_Rejected['variable']

sns.set(style="darkgrid")
Scenariosplot = sns.lineplot(
    x="Type", y="value",  data=Combined_Percentages_Dif)
plt.xlabel("Sampling Scenario")
plt.ylabel("Change Prop CFU Rej & Prop Lb Rejc ")


Scenariosplot = sns.boxplot(x="variable", y="value",
                            hue="Type", data=Combined_Percentages)
plt.xlabel("Sampling Scenario")
plt.ylabel("Proportion Rejection")
Scenariosplot.legend_.set_title('Type')


# Conbined Per Rejected CFU and lb Rejectet

Percentage_Rejected = df_L_Final_Per_melted
Percentage_Rejected['Type'] = "Percentage CFU Rejected"

Total_lb_Rejected = df_L_Final_RejWeight_melted
Total_lb_Rejected['Type'] = "Total lb Rejected"

Combined_PerCFUTotalWeight = Percentage_Rejected
Combined_PerCFUTotalWeight['value2'] = Total_lb_Rejected['value']


# Contamination Progression

def melting_type(df, typename):
    df_melt = pd.melt(df)
    df_melt["type"] = typename
    return df_melt


DProg_Baseline_melt = melting_type(DProg_Baseline, "Baseline")
DProg_PH4d_melt = melting_type(DProg_PH4d, "PH 4d")
DProg_PH4h_melt = melting_type(DProg_PH4h, "PH 4h")
DProg_PHInt_melt = melting_type(DProg_PHInt, "PH Int")
DProg_HTrad_melt = melting_type(DProg_HTrad, "H Trad")
DProg_HAgg_melt = melting_type(DProg_HAgg, "H Agg")
DProg_R_melt = melting_type(DProg_R, "Receiving")
DProg_FPSTrad_melt = melting_type(DProg_FPSTrad, "FPSTrad")
DProg_FPSAgg_melt = melting_type(DProg_FPSAgg, "FPS Agg")

all_list = [DProg_Baseline_melt, DProg_PH4d_melt, DProg_PH4h_melt, DProg_PHInt_melt,
            DProg_HTrad_melt, DProg_HAgg_melt, DProg_R_melt, DProg_FPSTrad_melt, DProg_FPSAgg_melt]
Frame_all = pd.concat(all_list)

sns.catplot(x="variable", y="value", hue="type", kind="bar",
            data=Frame_all, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression")
plt.xticks(rotation=70)


def melting_type(df, typename):
    df_melt = pd.melt(df)
    df_melt["type"] = typename
    return df_melt


DProg_Baseline_melt = melting_type(DProg_Baseline, "Baseline")
DProg_PH4d_melt = melting_type(DProg_PH4d, "PH 4d")
DProg_PH4h_melt = melting_type(DProg_PH4h, "PH 4h")
DProg_PHInt_melt = melting_type(DProg_PHInt, "PH Int")

DProg_HTrad_melt = melting_type(DProg_HTrad, "H Trad")
DProg_HAgg_melt = melting_type(DProg_HAgg, "H Agg")


DProg_R_melt = melting_type(DProg_R, "Receiving")
DProg_FPSTrad_melt = melting_type(DProg_FPSTrad, "FPSTrad")
DProg_FPSAgg_melt = melting_type(DProg_FPSAgg, "FPS Agg")

all_list = [DProg_Baseline_melt, DProg_PH4d_melt, DProg_PH4h_melt, DProg_PHInt_melt,DProg_HTrad_melt,
            DProg_HAgg_melt,DProg_R_melt,DProg_FPSTrad_melt,DProg_FPSAgg_melt]
Frame_all = pd.concat(all_list)

sns.catplot(x="variable", y="value", hue="type", kind="bar",
            data=Frame_all, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression")
plt.xticks(rotation=70)



# %%

# Additional Analysis

Weight_Rejected = {
    'PH4d': D_PH4d["PerRejectedWeight"],
    'PH4h': D_PH4h["PerRejectedWeight"],
    'PHInt': D_PHInt["PerRejectedWeight"],
}

df_Weight_Rejected = pd.DataFrame(Weight_Rejected)
df_Weight_Rejected_melted = pd.melt(df_Weight_Rejected)

Weight_Rejected_CFU = {
    'PH4d': D_PH4d["PerRejected at PH"],
    'PH4h': D_PH4h["PerRejected at PH"],
    'PHInt': D_PHInt["PerRejected at PH"],
}

df_CFU_Rejected = pd.DataFrame(Weight_Rejected_CFU)
df_CFU_Rejected_melted = pd.melt(df_CFU_Rejected)
CFU_Per = df_CFU_Rejected_melted.value


df_Weight_Rejected_melted["CFU"] = CFU_Per

sns.set(style="darkgrid")
Scenariosplot = sns.regplot(
    x="value", y="CFU", hue="variable", data=df_Weight_Rejected_melted)
plt.xlabel("Weight Rejected %")
plt.ylabel("CFU Rejected %")
plt.title("Testing Plan Efficacy")


g = sns.FacetGrid(df_Weight_Rejected_melted, col="variable",
                  margin_titles=True, height=4)
g.map(sns.regplot, "value", "CFU")
