# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:43:28 2022

@author: gareyes3
"""

####Factor Sensitivity Analysis

#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\reyes\Documents\GitHub\CPS-Farm-to-Facility')

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
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import sys
#import Trial_MainLoop_PH


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
                      Always_Washing_Optimized = False,
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
    
    Cont_Scen_no = np.random.choice([1,2,3])
    
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
        ScenCondz.C_Sampling == 1
    
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
    
    SCInputz.Always_Washing_Optimized = Always_Washing_Optimized
    
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

Baseline = scenario_function(Cont_Scen_no=1)
#Baseline_washing_opt = scenario_function(Cont_Scen_no=1,Washing = True,Always_Washing_Optimized= True)
Baseline_PHS4d = scenario_function(Cont_Scen_no=1, PHS4d = True) #
Baseline_PHS4h = scenario_function(Cont_Scen_no=1, PHS4h = True) #
Baseline_PHS4Int = scenario_function(Cont_Scen_no=1, PHSInt = True)#
Baseline_HS = scenario_function(Cont_Scen_no=1, HSTrad = True)#
Baseline_RS = scenario_function(Cont_Scen_no=1, RSTrad = True)
Baseline_FPSTrad = scenario_function(Cont_Scen_no=1, FPSTrad = True)
Baseline_CS = scenario_function(Cont_Scen_no=1, CSampling = True)

Baseline_washing = scenario_function(Cont_Scen_no=1, Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True )
Baseline_washing_PHS4D = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,PHS4d = True)
Baseline_washing_PHS4h = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,PHS4h = True)
Baseline_washing_PHSInt = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,PHSInt = True)
Baseline_washing_HS = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,HSTrad  = True)
Baseline_washing_RS = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,RSTrad = True)
Baseline_washing_FPS = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,FPSTrad = True)
Baseline_washing_CS = scenario_function(Cont_Scen_no=1,Holding=True, Pre_Cooling = True, PreS_Wash=True, Sanitation = True ,CSampling = True)

Baseline_Holding = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_Holding_PHS4D = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHS4d = True)
Baseline_Holding_PHS4h = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHS4h = True)
Baseline_Holding_PHSInt = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHSInt = True)
Baseline_Holding_HS = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,HSTrad  = True)
Baseline_Holding_RS = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,RSTrad = True)
Baseline_Holding_FPS = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,FPSTrad = True)
Baseline_Holding_CS = scenario_function(Cont_Scen_no=1,Washing = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,CSampling = True)

Baseline_PreS_Wash = scenario_function(Cont_Scen_no=1, Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True) #
Baseline_PreS_Wash_PHS4D = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,PHS4d = True)
Baseline_PreS_Wash_PHS4h = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,PHS4h = True)
Baseline_PreS_Wash_PHSInt = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,PHSInt = True)
Baseline_PreS_Wash_HS = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,HSTrad  = True)
Baseline_PreS_Wash_RS = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,RSTrad = True)
Baseline_PreS_Wash_FPS = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,FPSTrad = True)
Baseline_PreS_Wash_CS = scenario_function(Cont_Scen_no=1,Holding = True, Washing = True,Pre_Cooling = True, Sanitation = True,CSampling = True)

Baseline_Sanitation = scenario_function(Cont_Scen_no=1, Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True) #
Baseline_Sanitation_PHS4D = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,PHS4d = True)
Baseline_Sanitation_PHS4h = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,PHS4h = True)
Baseline_Sanitation_PHSInt = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,PHSInt = True)
Baseline_Sanitation_HS = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,HSTrad  = True)
Baseline_Sanitation_RS = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,RSTrad = True)
Baseline_Sanitation_FPS = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,FPSTrad = True)
Baseline_Sanitation_CS = scenario_function(Cont_Scen_no=1,Washing = True,Holding = True,Pre_Cooling = True, PreS_Wash=True,CSampling = True)

Baseline_Pre_Cooling = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True)
Baseline_Pre_Cooling_PHS4D = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,PHS4d = True)
Baseline_Pre_Cooling_PHS4h = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,PHS4h = True)
Baseline_Pre_Cooling_PHSInt = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,PHSInt = True)
Baseline_Pre_Cooling_HS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,HSTrad  = True)
Baseline_Pre_Cooling_RS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,RSTrad = True)
Baseline_Pre_Cooling_FPS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,FPSTrad = True)
Baseline_Pre_Cooling_CS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True, PreS_Wash=True, Sanitation = True,CSampling = True)


Baseline_AI = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True)
Baseline_AI_PHS4d = scenario_function(Cont_Scen_no=1,Washing = True,  Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHS4d = True)
Baseline_AI_PHS4h = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHS4h = True)
Baseline_AI_PHSInt = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,PHSInt =True)
Baseline_AI_HSTrad = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,HSTrad = True)
Baseline_AI_RS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,RSTrad = True)
Baseline_AI_FPSTrad = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,FPSTrad = True)
Baseline_AI_CS = scenario_function(Cont_Scen_no=1,Washing = True, Holding = True,Pre_Cooling = True, PreS_Wash=True, Sanitation = True,CSampling = True)





List_of_Names = ["Baseline NI", "Baseline NI PHS4D", "Baseline NI PHS4H", "Baseline NI PHS Int", "Baseline NI HS", "Baseline NI RS", "Baseline NI FPS", "Baseline NI CS",
                 "No Washing", "No Washing PHS4D", "No Washing PHS4H", "No Washing PHS Int", "No Washing HS", "No Washing RS", "No Washing FPS", "No Washing CS",
                 "No Holding", "No Holding PHS4D", "No Holding PHS4H", "No Holding PHS Int", "No Holding HS", "No Holding RS", "No Holding FPS", "No Holding CS",
                 "No PreWash", "No PreWash PHS4D", "No PreWash PHS4H", "No PreWash PHS Int", "No PreWash HS", "No PreWash RS", "No PreWash FPS", "No PreWash CS",
                 "No Precooling", "No Precooling PHS4D", "No Precooling PHS4H", "No Precooling PHS Int", "No Precooling HS", "No Precooling RS", "No Precooling FPS", "No Precooling CS",
                 "No Sanitation", "No Sanitation PHS4D", "No Sanitation PHS4H", "No Sanitation PHS Int", "No Sanitation HS", "No Sanitation RS", "No Sanitation FPS", "No Sanitation CS",                 
                 "Baseline  AI", "Baseline  AI PHS4D", "Baseline  AI PHS4H", "Baseline  AI PHS Int", "Baseline  AI HS", "Baseline  AI RS", "Baseline  AI FPS", "Baseline  AI CS",
    ]

List_of_systems = ["Baseline NI","Baseline NI","Baseline NI","Baseline NI","Baseline NI","Baseline NI","Baseline NI","Baseline NI",
                   "No Washing","No Washing","No Washing","No Washing","No Washing","No Washing","No Washing","No Washing",
                   "No Holding","No Holding","No Holding","No Holding","No Holding","No Holding","No Holding","No Holding",
                   "No PreWash", "No PreWash", "No PreWash", "No PreWash", "No PreWash", "No PreWash", "No PreWash", "No PreWash", 
                   "No Precooling","No Precooling","No Precooling","No Precooling","No Precooling","No Precooling","No Precooling","No Precooling",
                   "No Sanitation","No Sanitation","No Sanitation","No Sanitation","No Sanitation","No Sanitation","No Sanitation","No Sanitation",
                   "Baseline  AI", "Baseline  AI", "Baseline  AI", "Baseline  AI", "Baseline  AI", "Baseline  AI", "Baseline  AI", "Baseline  AI"
                   
                   ]

List_of_FS = [Baseline,Baseline_PHS4d,Baseline_PHS4h,Baseline_PHS4Int,Baseline_HS,Baseline_RS,Baseline_FPSTrad,Baseline_CS,
              Baseline_washing,Baseline_washing_PHS4D,Baseline_washing_PHS4h,Baseline_washing_PHSInt,Baseline_washing_HS,Baseline_washing_RS,Baseline_washing_FPS,Baseline_washing_CS,
              Baseline_Holding,Baseline_Holding_PHS4D,Baseline_Holding_PHS4h,Baseline_Holding_PHSInt,Baseline_Holding_HS,Baseline_Holding_RS,Baseline_Holding_FPS,Baseline_Holding_CS,
              Baseline_PreS_Wash,Baseline_PreS_Wash_PHS4D,Baseline_PreS_Wash_PHS4h,Baseline_PreS_Wash_PHSInt,Baseline_PreS_Wash_HS,Baseline_PreS_Wash_RS,Baseline_PreS_Wash_FPS,Baseline_PreS_Wash_CS,
              Baseline_Sanitation,Baseline_Sanitation_PHS4D,Baseline_Sanitation_PHS4h,Baseline_Sanitation_PHSInt,Baseline_Sanitation_HS,Baseline_Sanitation_RS,Baseline_Sanitation_FPS,Baseline_Sanitation_CS,
              Baseline_Pre_Cooling,Baseline_Pre_Cooling_PHS4D,Baseline_Pre_Cooling_PHS4h,Baseline_Pre_Cooling_PHSInt,Baseline_Pre_Cooling_HS,Baseline_Pre_Cooling_RS,Baseline_Pre_Cooling_FPS,Baseline_Pre_Cooling_CS,
              Baseline_AI,Baseline_AI_PHS4d, Baseline_AI_PHS4h,Baseline_AI_PHSInt,Baseline_AI_HSTrad,Baseline_AI_RS,Baseline_AI_FPSTrad,Baseline_AI_CS
              
    ]




List_of_FS = [
    Baseline,
    Baseline_washing,
    #Baseline_washing_opt,
    Baseline_Holding,
    Baseline_Pre_Cooling,
    Baseline_PreS_Wash,
    Baseline_Sanitation,
    Baseline_PHS4d,
    Baseline_PHS4h,
    Baseline_PHS4Int,
    Baseline_HS,
    Baseline_RS,
    Baseline_FPSTrad,
    Baseline_CS,
    Baseline_AI,
    Baseline_AI_PHS4d,
    Baseline_AI_PHS4h,
    Baseline_AI_PHSInt,
    Baseline_AI_HSTrad,
    Baseline_AI_RS,
    Baseline_AI_FPSTrad,
    Baseline_AI_CS
    ]



def get_FS(List):
    vector_of_FS = []
    for i in List:
        factorin = i[1]["After CS Samp"].sum()
        baselinein=List[0][1]["After CS Samp"].sum()
        FS = np.log10(factorin/baselinein)
        vector_of_FS.append(FS)
    return vector_of_FS

FS = get_FS(List_of_FS)

FSdf= pd.DataFrame({
    "FS": FS,
    "Name": List_of_Names,
    "System":List_of_systems
    })


FSdf.to_csv(path_or_buf = "C:\\Users\\gareyes3\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\Review 2\\FSdf-Final.csv") 
#FSdf.to_csv(path_or_buf = "C:\\Users\\gareyes3\\Box Sync\\CPS Project- Farm to Facility\\Papers\\CSV Data\\FSdf-Final.csv") 
    

