# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 10:56:37 2021

@author: gareyes3
"""
import ScenCondz
import ContCondz
import Funz
import numpy as np    


Random_Contam = False

#Initial Inputs for Field Setup and Challenges
Sensitivity_Analysis = False


N_Iterations= 10000
Field_Weight= 100000 #total lb in field
slot_weight = 10000 #weight of the sublot in field.
Partition_Weight = 50 #Partition weight for contamination Events, Partition to have better definition of field. 

slot_number = int(Field_Weight/slot_weight) #Number of sublots per field
Partition_Units = int(Field_Weight/Partition_Weight) #Number of partition units per field




#Sampling Inputs. 

#Pre_Harvest   
test_unit_PH = "Sublot"
sample_size_PH = 37.5 # (Input) g #Sample Size in grams for Pre Harvest
n_samples_slot_PH = 1 # (Input) Samples per sublot of product
No_Grabs_PH = 6 #Number of graps per sample. stratified per sample. 
Limit_PH = 0
RR_PH_Trad = "Lot" #Reject by Sublot




#Intesne pre-harvest sampling
test_unit_PH_Int = "Lot"
sample_size_PH_Int = 375
n_samples_lot_PH_Int = 10 # (Input) Samples per lot of product'
No_Grabs_PH_Int = 60
RR_PH_Int = "Lot"
test_unit_PH_Int = "Lot"
    

#Harvest Inputs: 
test_unit_H = "Sublot"
sample_size_H = 37.5 #g #Sample Size in grams
n_samples_slot_H = 1 # Samples per lot of product
No_Grabs_H = 6 
Limit_H = 0
RR_H_Trad  = "Lot"
#aggregative_prehavest sampling
if ScenCondz.HS_Agg ==True:
    n_samples_slot_H = 10 # (Input) Samples per lot of product'
    RR_H_Agg  = "Sublot"


#Receiving sampling:
test_unit_R = "PalletNo"
n_samples_pallet = 1 #samples taken per pallet
sample_size_R = 15 #375 #g #Sample Size in grams Receiving
No_Grabs_R = np.random.choice([2,3])
Limit_R =0
RR_R_Trad= "Lot"


#Finished Product Sampling: 
test_unit_FP = "HourProd"
sample_size_FP = 200 #g #Sample Size in grams
n_samples_FP = 1 #number of samples per lot final product
N_Packages_Samples = 4
RR_FP_Trad = "Lot"
Limit_FP =0
if ScenCondz.FPS_Agg ==True:
    n_samples_FP = 1
    Limit_FP =0
    RR_FP_Agg= "Lot"

#Receiving Sampling Field Pack/Also Final Product Sampling. 
test_unit_R_FP = "Lot"
n_samples_R_FP = 1
sample_size_R_FP = 375
No_GRabs_R_FP = 60
Limit_R_FP = 0
RR_R_FP_Trad  = "Lot"

#Retail Sampling
test_unit_CS = "PalletNoFS"
n_samples_CS = 1
sample_size_CS = 25
No_GRabs_CS = 3
Limit_CS = 0
RR_CS  = "Lot"

#%% Toggles
Washing_YN =False

Washing_Optimized = False

Always_Washing_Optimized = False

Pre_CoolingYN = False

Spray_WashYN =False

Sanitation_YN = False



