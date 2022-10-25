# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:47:47 2022

@author: gareyes3
"""

#%% Running a Scenario API


#%%
'''
Change these sources depending on where the model is saved in your computer.
You can still run the lines below, but remember to add you own path in order to 
be able to load the necessary libraries
'''
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\reyes\Documents\GitHub\CPS-Farm-to-Facility\Model')

#%%
'''
This chunk loads the necesary libraries to run the model. Not all of them are
necessary, but they are loaded in case they may be needed for other processes.
Ignore the triengle with exclamation marks, they just mean that the library may not
be in use. 
'''
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


#%% MAIN FUCTION THE RUNS THE MODEL WITH THE INTERVENTIONS LISTES

'''
This chunk loads the main function that makes running the model possible without
having to dig through multiple code files. The inputs can always be updated inside the 
function. 
The function allow you to select any of the interventions, as well as to select aany
of the sampling plans
'''
'''
Args: 
    Cont_Scen_no: The contamination scenario
        1: Uniform Contamination
        2:  1% cluster
        3: 10% cluster
        To change contamination levels or any other clustering access the Inputz.py
        file directly. Change the inputs direclty from there. 
    Total_Iterations: The total fields you would like to run
'''
     
def scenario_function(
                      #Model setup inputs
                      Cont_Scen_no,
                      Total_Iterations,
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
    
    #this here selects the contamination scenario
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
    
    #Management of static inputs
    SCInputz.N_Iterations =  Total_Iterations
    
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


#%%
'''
This section evaluated the important ourputs of the model that were shown in the analysis. 
While the model is capable of doing a lot more, this are the outputs that were caputured
for the analysis. 

As you can se below you can run many scenarios based on the system you would like to represent. 

For the purpose of this model we can compare 2 scenarios. 
'''
My_Scenario = scenario_function(Cont_Scen_no= 1,
                                Total_Iterations= 100)

My_Senario_2 = scenario_function(Cont_Scen_no= 1,
                                Total_Iterations= 100, 
                                Washing= True,
                                PHS4d=True)
