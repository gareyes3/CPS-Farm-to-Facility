# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 16:00:57 2021

@author: Gustavo Reyes
"""

#Sensityvity Analysis

#Inputs to Analyze: 
    #1: Initial Contamination
    #2: Contamination Clustering Level
    #3: Sample Size
    #4: Nuber of Grabs
    
    #5: Washing Step 1 or 0
    #6: 

#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
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

#reload(MainModel3z)

#%%

#%% Sensitivity. 

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Changing Sens_analysis inputs. 
#Turning on the Sensitivity Analysis


#Sampling Condition

SCInputz.Sensitivity_Analysis = True


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF = Main_Mod_Outs[1]
ProgDF = Main_Mod_Outs[0]
PropProgDF = Main_Mod_Outs[2]    
Sensitivity_Df = Main_Mod_Outs[5]    
    

#%%
Sensitivity_Df.to_csv("C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\SensitivityDF.csv")

#Sensitivity_Df.to_csv("C:\\Users\\Gustavo Reyes\\Documents\\GitHubFiles\\CPS-Farm-to-Facility\\Model\\SensitivityDF2.csv")
Sensitivity_Df.to_csv("C:\\Users\\gareyes3\\Documents\\GitHub\\CPS-Farm-to-Facility\\Model\\SensitivityDF2.csv")
