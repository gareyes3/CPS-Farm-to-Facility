# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 08:32:41 2021

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
import random
import Dictionariez
import SCInputsValidation
sys.path
#sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append(
    'C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

#%%
def F_MainLoop_Validation():
    
    #DataCollection DataFrame for outputs.  
    df_Output_Contprog = Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Progression, SCInputz.N_Iterations) #Progression Dataframe
    df_Output_PH =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe pre harvest


    for  i in range(SCInputz.N_Iterations):
        Iteration_In = i
        print(i)
        reload(Inputz)
        
    
        #Adding Contmination to the Field
        
        #STEP 0 CONTAMINATION SCENARIOS  ----------------------------------------------------------------------------------------------------

        #Creation of the Data Frame to Track: 
        df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                          Field_Weight = SCInputz.Field_Weight, 
                          slot_number = SCInputz.slot_number)
            
            
        #Adding Contamination depending on challenge Systematic Sampling
        if ContCondz.Systematic_C == True:
            df = ContScen.F_systematic_C(df=df, Hazard_lvl= SCInputsValidation.Hazard_Level,
                                         No_Cont_Clusters =1,
                                         Cluster_Size= 100_000,
                                         Partition_Weight = SCInputz.Partition_Weight)
            
        

    
        #STEP 1 PREHARVEST ------------------------------------------------------------------------------------------------------------------
        
        #print("Dieoff",LV_Die_Off_CE_PHS)
            
        LO_Cont_B_PH = sum(df.CFU) #Contamination before rejection sampling
        #print("before",LO_Cont_B_PH)
        LO_Weight_B_PH = sum(df.Weight)
        
        #Contprog Before Pre-Harvest
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Bef Pre-Harvest Samp", 
                                                     i =Iteration_In )
        
        #Sampling at Pre-Harvest
        if ScenCondz.PH_Sampling ==True: #If function to turn off Pre-Harvest Sampling
            if ScenCondz.PHS_Int ==True: #Intense pre harvest sampling
                df = Funz.F_Sampling_2(df =df,Test_Unit ="Lot", 
                                              NSamp_Unit = 1, 
                                              Samp_Size =SCInputz.sample_size_PH, 
                                              Partition_Weight =SCInputz.Partition_Weight, 
                                              NoGrab =SCInputz.No_Grabs_PH )
            elif ScenCondz.PHS_4d==True or ScenCondz.PHS_4h == True :
            #Pre-Harvest Sampling, Traditional
                 df = Funz.F_Sampling_2(df =df,
                                        Test_Unit ="Lot", 
                                         NSamp_Unit = 1, 
                                         Samp_Size =SCInputz.sample_size_PH, 
                                         Partition_Weight =SCInputz.Partition_Weight, 
                                         NoGrab =SCInputz.No_Grabs_PH)
            
            
        
        Listz.List_BPHS_CFU.append( LO_Cont_B_PH) #List of contamination before sampling
        
        #Filtering out the Rejected lots, Pre-Harvest
        if ScenCondz.PHS_Int ==True: #Rejection intense
           df= Funz.F_Rejection_Rule3(df =df, Test_Unit = SCInputz.RR_PH_Int, limit = SCInputz.Limit_PH)  
        else:  #Rejection normal
            df=Funz.F_Rejection_Rule3(df =df, Test_Unit = SCInputz.RR_PH_Trad, limit = SCInputz.Limit_PH) 
           
        #print("PH", sum(df["CFU"]))
        #Contprog After Pre-Harvest
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Aft Pre-Harvest Samp", 
                                                     i =Iteration_In )
                           
    
        
        #Outputs Function, Instead of collecting outputs. 
        df_Output_PH = Dictionariez.Output_Collection_Final(df = df, 
                                                            outputDF = df_Output_PH, 
                                                            Step = "PH", 
                                                            Cont_Before = LO_Cont_B_PH, 
                                                            Weight_Before = LO_Weight_B_PH, 
                                                            i = Iteration_In, 
                                                            Niterations = SCInputz.N_Iterations)
        
    df_outputs = df_Output_PH
    
    
    outputs = [df_Output_Contprog, df_outputs]

        
    return outputs


