# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 09:12:08 2021

@author: gareyes3
"""
#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

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
import Trial_MainLoop_PH


reload(MainModel3z)
reload(Inputz)
reload(SCInputz)
reload(ScenCondz)
#%%

#Scenarios Creation
#Systematic Failure from irrigation. 100% of field affected at 1CFU/g 

# 1 CFU/lb Contmination Level
One_CFU_lb = 100000 #CFUs

#Uniform Contmaination in the Field. 
# Contamination Challenges
ContCondz.Background_C = False
ContCondz.Point_Source_C =True
ContCondz.Systematic_C = False

# Harvester Contamination
ContCondz.Crew_C = False
ContCondz.Harvester_C = False

# Processing equipment
ContCondz.PE_C = False
ContCondz.PE_Cont_Loc = False,  # 1,2,3,4,5
# 1 = Shredder, #2 = Belt, #3 = Washing, #4 Shaker, #5Centrifuge
ContCondz.Pack_C = False


#%% Scenario 1. No intervention Strategies, Natural Field and No Washing. 

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Management of Control: 

#Holding Time: 
ScenCondz.Holding_Time= False #0,4,8 Holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False 

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False

#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF = Main_Mod_Outs[1]
ProgDF = Main_Mod_Outs[0]
PropProgDF = Main_Mod_Outs[2]
#FinalConts = Main_Mod_Outs[6]

 


#%% Scenario 2. No Sampling, Only Intervention Strategy is Washing

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = True

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_Wash = Main_Mod_Outs[1]
ProgDF_Base_Wash = Main_Mod_Outs[0]
PropProgDF_Base_Wash = Main_Mod_Outs[2]
#FinalConts_Base_Wash = Main_Mod_Outs[6]



#%% Scenario 3 PH Sampling 4D, Washing off , and normal baseline contamination reduction
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot
    #Sample Mass: 365g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 1
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 1  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
   
#Only on when sampling,
SCInputz.test_unit_PH = "Sublot"
SCInputz.sample_size_PH = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_slot_PH = 1 # (Input) Samples per sublot of product
SCInputz.No_Grabs_PH = 60
SCInputz.Limit_PH = 0
SCInputz.RR_PH_Trad = "Lot" #Reject by Sublot


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_PHS4d = Main_Mod_Outs[1]
ProgDF_Base_PHS4d = Main_Mod_Outs[0]
PropProgDF_Base_PHS4d = Main_Mod_Outs[2]
#FinalConts_Base_PHS4d = Main_Mod_Outs[6]




#%% Scenario 4 PH Sampling 4h, Washing off , and normal baseline contamination reduction
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Sublot"
    #Sample Mass: 367g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 1
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 1  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
   
#Only on when sampling,
SCInputz.test_unit_PH = "Sublot"
SCInputz.sample_size_PH = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_slot_PH = 1 # (Input) Samples per sublot of product
SCInputz.No_Grabs_PH = 60
SCInputz.Limit_PH = 0
SCInputz.RR_PH_Trad = "Lot" #Reject by Sublot


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_PHS4h = Main_Mod_Outs[1]
ProgDF_Base_PHS4h = Main_Mod_Outs[0]
PropProgDF_Base_PHS4h = Main_Mod_Outs[2]
#FinalConts_Base_PHS4h = Main_Mod_Outs[6]




#%% Scenario 5 PH Sampling Int, Washing off , and normal baseline contamination reduction
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot"
    #Sample Mass: 367g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 1
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 1  # Scenario 3
   
#Only on when sampling,
SCInputz.test_unit_PH = "Lot"
SCInputz.sample_size_PH = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_lot_PH = 10 # (Input) Samples per sublot of product
SCInputz.No_Grabs_PH = 60
SCInputz.Limit_PH = 0
SCInputz.RR_PH_Int = "Lot" #Reject by Sublot


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_PHSInt = Main_Mod_Outs[1]
ProgDF_Base_PHSInt = Main_Mod_Outs[0]
PropProgDF_Base_PHSInt = Main_Mod_Outs[2]
#FinalConts_Base_PHSInt = Main_Mod_Outs[6]




#%% Harvest Sampling Traditional 
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot"
    #Sample Mass: 367g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 1
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
#
ScenCondz.HS_Trad = 1  # Scenario 1


   
#Only on when sampling,
SCInputz.test_unit_H = "Sublot"
SCInputz.sample_size_H = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_slot_H = 1 # (Input) Samples per sublot of product
SCInputz.No_Grabs_H = 60
SCInputz.Limit_H = 0
SCInputz.RR_H_Trad = "Lot" #Reject by Sublot



#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_HTrad = Main_Mod_Outs[1]
ProgDF_Base_HTrad = Main_Mod_Outs[0]
PropProgDF_Base_HTrad = Main_Mod_Outs[2]
#FinalConts_Base_HTrad = Main_Mod_Outs[6]


#%% Receiving Sampling
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot"
    #Sample Mass: 367g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 1
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
# HArvest Sampling
ScenCondz.HS_Trad = 0  # Scenario 1

#Receiving Sampling


   
#Only on when sampling,
SCInputz.test_unit_R = "PalletNo"
SCInputz.sample_size_R = 150 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_pallet = 1 # (Input) Samples per sublot of product
SCInputz.No_Grabs_R = 24
SCInputz.Limit_R = 0
SCInputz.RR_R_Trad = "Lot" #Reject by Sublot



#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_RTrad = Main_Mod_Outs[1]
ProgDF_Base_RTrad = Main_Mod_Outs[0]
PropProgDF_Base_RTrad = Main_Mod_Outs[2]
#FinalConts_Base_RTrad = Main_Mod_Outs[6]




#%%  Finsihed Product Sampling

    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot"
    #Sample Mass: 367g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False #Defaults

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 1
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
# HArvest Sampling
ScenCondz.HS_Trad = 0  # Scenario 1

#Receiving Sampling
#Finished Product
ScenCondz.FPS_Trad = 1


   
#Only on when sampling,
SCInputz.test_unit_FP = "Lot"
SCInputz.sample_size_FP = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_FP = 10 # (Input) Samples per sublot of product
SCInputz.N_Packages_Samples = 60
SCInputz.Limit_FP = 0
SCInputz.RR_R_Trad = "Lot" #Reject by Sublot



#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_FPTrad = Main_Mod_Outs[1]
ProgDF_Base_FPTrad = Main_Mod_Outs[0]
PropProgDF_Base_FPTrad = Main_Mod_Outs[2]
#FinalConts_Base_FPTrad = Main_Mod_Outs[6]



#%%
#Scenario 5 Holding Time on Everything Else Off

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= True #Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN = False 

# Turning of Washing. 
SCInputz.Washing_YN = False

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster


#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_Holding = Main_Mod_Outs[1]
ProgDF_Base_Holding = Main_Mod_Outs[0]
PropProgDF_Base_Holding = Main_Mod_Outs[2]
#FinalConts_Base_Holding = Main_Mod_Outs[6]



#%%
#Scenario #6 Only Intervention if Pre-Cooling

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False#Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN =True

# Turning of Washing. 
SCInputz.Washing_YN = False

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =False

#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster


#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_Precool = Main_Mod_Outs[1]
ProgDF_Base_Precool = Main_Mod_Outs[0]
PropProgDF_Base_Precool = Main_Mod_Outs[2]
#FinalConts_Base_Precool = Main_Mod_Outs[6]




#%%
#Scenario #7 Only Intervention if Pre-Harvest SprayWash

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= False#Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN =False
# Turning of Washing. |
SCInputz.Washing_YN = False

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =True


#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster


#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 0
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 0  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_Base_Prewash = Main_Mod_Outs[1]
ProgDF_Base_Prewash = Main_Mod_Outs[0]
PropProgDF_Base_Prewash = Main_Mod_Outs[2]
#FinalConts_Base_Prewash = Main_Mod_Outs[6]



#%% Scenario 8 PH Sampling 4D, Washing ON ,Pre-cooling on and everything else on and normal baseline contamination reduction
    #PH Sample: 1 Sample/Sublot, 
    #Rejection Rule: "Lot
    #Sample Mass: 365g per sublot
    # GRabs per Sublot = N60

reload(SCInputz)  # Reload Inputz
reload(Listz)  # Reload Lists

#Holding Time: 
ScenCondz.Holding_Time= True#Defaults to 1 Day holding Time

#Turning off Pre-Cooling
SCInputz.Pre_CoolingYN =True

# Turning of Washing. 
SCInputz.Washing_YN = True

#Harvest Pre-Wash: 
SCInputz.C_Spray_HYN =True

#SCInputz.BGHazard_lvl = One_CFU_lb

SCInputz.PSHazard_lvl = One_CFU_lb  #CFU # background contamination
SCInputz.PSCluster_Size = np.random.choice(np.arange(1000,100001,1000)) #lb 1%
SCInputz.PSNo_Cont_Clusters = 1 #4 1 1% Cluster

#Sampling Condition
# Sampling Conditions, Baseline all conditions are off
ScenCondz.Baseline_Sampling = 0  # all others must be 0if this one is 1
ScenCondz.PH_Sampling = 1
ScenCondz.H_Sampling = 0
ScenCondz.R_Sampling = 0
ScenCondz.FP_Sampling = 0
# Pre_Harvest 4 Days
ScenCondz.PHS_4d = 1  # Scenario 1
ScenCondz.PHS_4h = 0  # Scenario 2
ScenCondz.PHS_Int = 0  # Scenario 3
   
#Only on when sampling,
SCInputz.test_unit_PH = "Sublot"
SCInputz.sample_size_PH = 375 # (Input) g #Sample Size in grams for Pre Harvest
SCInputz.n_samples_slot_PH = 1 # (Input) Samples per sublot of product
SCInputz.No_Grabs_PH = 60
SCInputz.Limit_PH = 0
SCInputz.RR_PH_Trad = "Sublot" #Reject by Sublot


#Running The Model.
Main_Mod_Outs = MainModel3z.F_MainLoop()

#Getting the outputs from the function.
OutputDF_All = Main_Mod_Outs[1]
ProgDF_All = Main_Mod_Outs[0]
PropProgDF_All = Main_Mod_Outs[2]
#FinalConts_All = Main_Mod_Outs[6]



#%% Data Cleaning

#Plotting the current data. 
List_of_OutputDFs = [OutputDF,
                     OutputDF_Base_Wash,
                     OutputDF_Base_PHS4d,
                     OutputDF_Base_PHS4h,
                     OutputDF_Base_PHSInt,
                     OutputDF_Base_HTrad,
                     OutputDF_Base_RTrad,
                     OutputDF_Base_FPTrad,
                     OutputDF_Base_Holding,
                     OutputDF_Base_Precool,
                     OutputDF_Base_Prewash,
                     OutputDF_All]
List_of_ProgDfs = [ProgDF,
                     ProgDF_Base_Wash,
                     ProgDF_Base_PHS4d,
                     ProgDF_Base_PHS4h,
                     ProgDF_Base_PHSInt,
                     ProgDF_Base_HTrad,
                     ProgDF_Base_RTrad,
                     ProgDF_Base_FPTrad,
                     ProgDF_Base_Holding,
                     ProgDF_Base_Precool,
                     ProgDF_Base_Prewash,
                     ProgDF_All]


List_of_Propprog  = [PropProgDF,
                     PropProgDF_Base_Wash,
                     PropProgDF_Base_PHS4d,
                     PropProgDF_Base_PHS4h,
                     PropProgDF_Base_PHSInt,
                     PropProgDF_Base_HTrad,
                     PropProgDF_Base_RTrad,
                     PropProgDF_Base_FPTrad,
                     PropProgDF_Base_Holding,
                     PropProgDF_Base_Precool,
                     PropProgDF_Base_Prewash,
                     PropProgDF_All]
'''
List_Contaminated_Packages = [FinalConts,
                     FinalConts_Base_Wash,
                     FinalConts_Base_PHS4d,
                     FinalConts_Base_PHS4h,
                     FinalConts_Base_PHSInt,
                     FinalConts_Base_HTrad,
                     FinalConts_Base_RTrad,
                     FinalConts_Base_FPTrad,
                     FinalConts_Base_Holding,
                     FinalConts_Base_Precool,
                     FinalConts_Base_Prewash,
                     FinalConts_All]
'''

List_of_Final_Conts = [x["Final Product Facility"] for x in List_of_ProgDfs]
Column_Names = "Baseline Wash PH4D PH4h PHInt HTrad RTrad FPTrad Holding Precool PreWash All".split()

#Removing outliers from the contamination before.
def remove_outliers(list_a):
    new_list=[]
    for i in list_a:
        Q1 = i.quantile(0.25)
        Q3 = i.quantile(0.75)
        IQR = Q3 - Q1
        new = i[~((i < (Q1 - 1.5 * IQR)) |(i > (Q3 + 1.5 * IQR)))]
        new_list.append(new)
    return new_list

List_FinalCont_NoOuts = remove_outliers(list_a = List_of_Final_Conts)

#CREATING THE DATAFRAME for contaminations without outliers
Final_Contdf = pd.concat(List_FinalCont_NoOuts, axis = 1)
Final_Contdf.columns = Column_Names
#melting that dataframe
Final_Contdf_melted = Final_Contdf.melt()

'''
#Contaminated Packages
#strp 1, make them into a list of 
aaaaa=pd.concat(List_Contaminated_Packages[0], axis =1)

Non_Zeros=List_Contaminated_Packages[0][2][List_Contaminated_Packages[0][2]>0]
Prop_Cont = len(Non_Zeros)/len(List_Contaminated_Packages[0][2])
CFU_g_Cont=Non_Zeros/5
CFU_g_Cont.describe()
'''
    

#%% Data Analysis and Plotting

#Boxplot Y axis log scale adjusted. 
H=sns.catplot(x="variable", y="value", kind = "box" ,
            data=Final_Contdf_melted)
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")
plt.yscale('log')
plt.title("Log10 CFU Final Contaminations")
plt.xticks(rotation=70)

#same plot but instead we use the means
H=sns.catplot(x="variable", y="value" ,
            data=Final_Contdf_melted)
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")
plt.title("Log10 CFU Final Contaminations")
plt.xticks(rotation=70)


#%% DATA ANALYSIS

#COMPARING THE mean and 95% CI Reduction between two groups. 
import numpy as np
import math
import scipy.stats as st
import statsmodels.stats.api as sms
#Cmparing Pairs to Baseline Model


#Function with confidence interval.
def mean_CI_RED(a,b):
    Reductio_Change_Wash = a - b
    mean_r = Reductio_Change_Wash.mean()
    CI_r =st.t.interval(alpha=0.95, df=len(Reductio_Change_Wash)-1, loc=np.mean(Reductio_Change_Wash), scale=st.sem(Reductio_Change_Wash)) 
    outputs = [ mean_r,CI_r]
    return outputs

def mean_CI_RED_per(a,b):
    Reductio_Change_Wash = (a - b)/a
    mean_r = Reductio_Change_Wash.mean()
    CI_r =st.t.interval(alpha=0.95, df=len(Reductio_Change_Wash)-1, loc=np.mean(Reductio_Change_Wash), scale=st.sem(Reductio_Change_Wash)) 
    outputs = [ mean_r,CI_r]
    return outputs

def mean_CI_ONE(Array):
    mean = Array.mean()
    CI = sms.DescrStatsW(Array).tconfint_mean()
    #log_mean = mean.log10()
    #CI_log = CI.log10()
    return [mean,CI]

def mean_CI_ONE_log(Array):
    mean = Array.mean()
    CI = sms.DescrStatsW(Array).tconfint_mean()
    log_mean = np.log10(mean)
    CI_log = np.log10(CI)
    return [log_mean,CI_log]


#Getting the mean log reduction per Name
[mean_CI_ONE_log(x) for x in List_of_Final_Conts]
[mean_CI_ONE(x) for x in List_of_Final_Conts]
(7049-1828)/7049
2735.6133333333332
(7049-2735.6133333333332)/7049


[mean_CI_ONE(x)[0]/8776.7325 for x in List_of_Final_Conts]


#desnsity plot
h=sns.displot( data =Final_Contdf_melted, 
            x = "value" , 
            col = "variable", 
            col_wrap=3,
             stat = "probability",
             bins = 30,
            facet_kws=dict(sharey=False,sharex= False))
plt.suptitle("Distribution of CFUs in Final Product",) 

def specs(x, **kwargs):
    plt.axvline(x.mean(), c='red', ls='-', lw=2.5)
    plt.axvline(x.median(), c='orange', ls='--', lw=2.5)

h.map(specs,"value" )



#Statistical analysis. 
from scipy import stats
import scikit_posthocs as sp
stats.kruskal(List_of_Final_Conts)
Dunn_posthoc=sp.posthoc_dunn(List_of_Final_Conts, p_adjust = 'bonferroni')  











Reductio_Change_Wash = Final_CFU_Baseline["Final Product Facility"] - Final_CFU_Base_Wash["Final Product Facility"]
Reductio_Change_Wash.mean()
st.t.interval(alpha=0.95, df=len(Reductio_Change_Wash)-1, loc=np.mean(Reductio_Change_Wash), scale=st.sem(Reductio_Change_Wash)) 

st.ttest_ind(Final_CFU_Baseline["Final Product Facility"], Final_CFU_Base_Wash["Final Product Facility"], equal_var=False)

st.normaltest(Final_CFU_Baseline["Final Product Facility"])
st.normaltest(Final_CFU_Base_Wash["Final Product Facility"])


#Getting the mean. 
Final_Compared.groupby("Type").mean()


#Log Change in CFU population
Final_Compared["CFU Change"] = 100000- Final_Compared["Final Product Facility"]
Final_Compared["CFU Change Log"]=np.log10(Final_Compared["CFU Change"])
Final_Compared["CFU Change Per"]=Final_Compared["CFU Change"]/100000

#Confidence interval Creation. you can change the columns desired. 
import math


stats = Final_Compared.groupby(["Type"])["CFU Change Per"].agg(['mean', 'count', 'std'])
print(stats)
print('-'*30)

ci95_hi = []
ci95_lo = []

for i in stats.index:
    m, c, s = stats.loc[i]
    ci95_hi.append(m + 1.96*s/math.sqrt(c))
    ci95_lo.append(m - 1.96*s/math.sqrt(c))

stats['ci95_hi'] = ci95_hi
stats['ci95_lo'] = ci95_lo
print(stats)



##3compring means with CI

List_of_series = [x["Final Product Facility"] for x in List_of_dfs ]
List_of_Names = "Baseline Wash PH4d PH4h Holding Precool HWash All".split()
Final_Compared_Cols = pd.concat(List_of_series,axis=1)
Final_Compared_Cols.columns = List_of_Names

import numpy as np, statsmodels.stats.api as sms


cm = sms.CompareMeans(sms.DescrStatsW(Final_Compared_Cols['Baseline']), sms.DescrStatsW(Final_Compared_Cols['Wash']))
print (cm.tconfint_diff(usevar='unequal'))

#%% PLotting

#Initial levels represent the initial contamination level before or after PH

    #2-8 days before harvest. Triangular distribution.     
H=sns.catplot(x="Type", y="Final Product Facility", 
            data=Final_Compared)
#H.map(sns.swarmplot, 'Type', 'Final Product Facility', color='k')
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")
plt.title("CFUs Initial vs Strategies")
plt.xticks(rotation=70)

H=sns.catplot(x="Type", y="Final Product Facility", kind = "box" ,
            data=Final_Compared)
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")
plt.yscale('log')
plt.title("CFUs Initial vs Strategies")
plt.xticks(rotation=70)


#Log Change bar bar chart
H=sns.catplot(x="Type", y="CFU Change Log", kind = "box" ,
            data=Final_Compared)
plt.xlabel("Sampling Scenario")
plt.ylabel("Total CFUs")

plt.title("Log10 CFU Change per intervention")
plt.xticks(rotation=70)

##Progression plot for baseline scenario: 
    
ProgDF_melt = pd.melt(ProgDF)
sns.catplot(x="variable", y="value", kind="bar",
            data=ProgDF_melt, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression Through System")
plt.xticks(rotation=70)


sns.catplot(x="variable", y="value", kind="box",
            data=ProgDF_melt, height=4, aspect=12/4)
plt.xlabel("Sampling Scenario")
plt.ylabel("CFU in System")
plt.title("Contamination Progression Through System")
plt.xticks(rotation=70)




#Desnity Plots for Final Contamination. 
#g = sns.FacetGrid(Final_Compared, col="Type", col_wrap=3)
#g.map_dataframe(sns.histplot, x="Final Product Facility")


h=sns.displot( data =Final_Compared, 
            x = "Final Product Facility" , 
            col = "Type", 
            col_wrap=3,
             stat = "probability",
             bins = 30,
            facet_kws=dict(sharey=False,sharex= False))
plt.suptitle("Distribution of CFUs in Final Product",) 

def specs(x, **kwargs):
    plt.axvline(x.mean(), c='red', ls='-', lw=2.5)
    plt.axvline(x.median(), c='orange', ls='--', lw=2.5)

h.map(specs,"Final Product Facility" )


#sns.displot(Final_Compared, x="Final Product Facility", col="Type", kind="kde", common_norm=False, col_wrap=3)

#Statistical analysis. 
from scipy import stats
import scikit_posthocs as sp
stats.kruskal(A,B,C,G,D,E,F,Z)
data = [A,B,C,G,D,E,F,Z]
Dunn_posthoc=sp.posthoc_dunn(data, p_adjust = 'bonferroni')  

#%%
#Percent Contaminated at each stage. 
#Adding Types To merge
Final_Prop_Baseline = Final_Prop_Baseline.to_frame()
Final_Prop_Baseline["Type"] = "No Sampling, No Wash"
x_prp = Final_Prop_Baseline["PropCont_A_FP"]

#totalcontaminated pakages
x_totCont = PropProgDF["TotalCont_A_FP"].to_frame()
x_totCont["Type"] = "No Sampling, No Wash"

#-------------

Final_Prop_Base_Wash = Final_Prop_Base_Wash.to_frame()
Final_Prop_Base_Wash["Type"] = "No Sampling, Wash"
y_prp = PropProgDF_Base_Wash["PropCont_A_FP"]

y_totCont = PropProgDF_Base_Wash["TotalCont_A_FP"].to_frame()
y_totCont["Type"] = "No Sampling, Wash"


#-------------------
Final_Prop_Base_PHS4d=Final_Prop_Base_PHS4d.to_frame()
Final_Prop_Base_PHS4d["Type"] = "Sampling, No Wash"
z_prp = PropProgDF_Base_PHS4d["PropCont_A_FP"]

z_totCont = PropProgDF_Base_PHS4d["TotalCont_A_FP"].to_frame()
z_totCont["Type"] = "Sampling, No Wash"

#----------------------
Final_Prop_All = Final_Prop_All.to_frame()
Final_Prop_All["Type"] = "Sampling, Wash"
xy_prp = PropProgDF_All["PropCont_A_FP"]

xy_totCont = PropProgDF_All["TotalCont_A_FP"].to_frame()
xy_totCont["Type"] = "Sampling, Wash"



#Combining them for final proportions
Final_Compared_Prop=pd.concat([Final_Prop_Baseline, 
                          Final_Prop_Base_Wash,
                          Final_Prop_Base_PHS4d,
                          Final_Prop_All], 
                         axis=0, 
                         ignore_index=True)

#sum of all 

Final_Compared_TotCont=pd.concat([x_totCont, 
                             y_totCont,
                             z_totCont,
                             xy_totCont], 
                             axis=0, 
                             ignore_index=True)



#Desnity Plots for proportion of contaminated final packages
#g = sns.FacetGrid(Final_Compared_Prop, col="Type", col_wrap=3)
#g.map_dataframe(sns.histplot, x="PropCont_A_FP", binwidth=0.005, stat = "probability")


g=sns.displot( data =Final_Compared_Prop, 
            x = "PropCont_A_FP" , 
            col = "Type", 
            col_wrap=3,
             stat = "probability",
            facet_kws=dict(sharey=False,sharex= False))

def specs(x, **kwargs):
    plt.axvline(x.mean(), c='red', ls='-', lw=2.5)
    plt.axvline(x.median(), c='orange', ls='--', lw=2.5)

g.map(specs,"PropCont_A_FP" )




#Barplot for proportion of contaminated final packages. 
j = sns.barplot(data = Final_Compared_Prop, y = "Type", x = "PropCont_A_FP")

#Bar plot for total contaminated final packages.
j = sns.barplot(data = Final_Compared_TotCont, y = "Type", x = "TotalCont_A_FP")

#%%

#Total CFUs, for all finished product bags. 

#Baseline no Sampling
FinalConts  = [item for sublist in FinalConts for item in sublist]
dfFinalConts = pd.DataFrame(FinalConts, columns=["CFU"])
dfFinalConts["CFU_g"] = dfFinalConts["CFU"]/(5*454)
sns.displot(dfFinalConts["CFU_g"], bins=30,stat = "probability")
cont_Packages=sum(dfFinalConts["CFU_g"]>0)
Non_contPackages= sum(dfFinalConts["CFU_g"]==0)
RejectedPackages = ((SCInputz.Field_Weight/Inputz.Pack_Weight_FP)*SCInputz.N_Iterations -len(dfFinalConts["CFU_g"]))
cont_Packages/(Non_contPackages+RejectedPackages)*100


#Washing
FinalConts_Base_Wash  = [item for sublist in FinalConts_Base_Wash for item in sublist]
dfFinalConts_Base_Wash = pd.DataFrame(FinalConts_Base_Wash, columns=["CFU"])
dfFinalConts_Base_Wash["CFU_g"] = dfFinalConts_Base_Wash["CFU"]/(5*454)
sns.displot(dfFinalConts_Base_Wash["CFU_g"], bins=30,stat = "probability")
cont_Packages=sum(dfFinalConts_Base_Wash["CFU_g"]>0)
Non_contPackages= sum(dfFinalConts_Base_Wash["CFU_g"]==0)
RejectedPackages = ((SCInputz.Field_Weight/Inputz.Pack_Weight_FP)*SCInputz.N_Iterations -len(dfFinalConts_Base_Wash["CFU_g"]))
cont_Packages/(Non_contPackages+RejectedPackages)*100


#Sampling No Wash
FinalConts_Base_PHS4d  = [item for sublist in FinalConts_Base_PHS4d for item in sublist]
dfFinalConts_Base_PHS4d = pd.DataFrame(FinalConts_Base_PHS4d, columns=["CFU"])
dfFinalConts_Base_PHS4d["CFU_g"] = dfFinalConts_Base_PHS4d["CFU"]/(5*454)
sns.displot(dfFinalConts_Base_PHS4d["CFU_g"], bins=30,stat = "probability")
cont_Packages=sum(dfFinalConts_Base_PHS4d["CFU_g"]>0)
Non_contPackages= sum(dfFinalConts_Base_PHS4d["CFU_g"]==0)
RejectedPackages = ((SCInputz.Field_Weight/Inputz.Pack_Weight_FP)*SCInputz.N_Iterations -len(dfFinalConts_Base_PHS4d["CFU_g"]))
cont_Packages/(Non_contPackages+RejectedPackages)*100

#Sampling and Wash.
FinalConts_All  = [item for sublist in FinalConts_All for item in sublist]
dfFinalConts_All = pd.DataFrame(FinalConts_All, columns=["CFU"])
dfFinalConts_All["CFU_g"] = dfFinalConts_All["CFU"]/(5*454)
sns.displot(dfFinalConts_All["CFU_g"], bins=30,stat = "probability")
cont_Packages=sum(dfFinalConts_All["CFU_g"]>0)
Non_contPackages= sum(dfFinalConts_All["CFU_g"]==0)
RejectedPackages = ((SCInputz.Field_Weight/Inputz.Pack_Weight_FP)*SCInputz.N_Iterations -len(dfFinalConts_All["CFU_g"]))
cont_Packages/(Non_contPackages+RejectedPackages)*100

