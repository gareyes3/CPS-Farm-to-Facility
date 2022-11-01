

#%%
import sys
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

#%% 
#Libraries, Modules
import pandas as pd  #for working
#import seaborn as sns
#from matplotlib import pyplot as plt
import Funz #Functions that we might need
import ContScen #Contamination Scenarios
import Listz #Empty lists that I need to things thoughut process, #Update
import InFunz #Input functions
import ScenCondz #Scenario Conditions
import ContCondz #contamination conditions
import Inputz #Random Inputz. 
import SCInputz #Inputz that can be changed later on through analysis
import Dictionariez #for dataframe creation
from importlib import reload  
import numpy as np
#import Sensitivity_analysis
import random

reload(Listz)

#%%

def F_MainLoop():
    '''
    This main model puts all the functions together into the Leafy green process. 
    
    Returns the 4 main outputs of the model
    0.the output dataframe
    1. 
    -------
    TYPE
        DESCRIPTION.

    '''
    
    #DataCollection DataFrame for outputs.
    #Progression Data Frames
    df_Output_Contprog = Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Progression, SCInputz.N_Iterations) #Progression Dataframe
    df_Output_Propprog = Dictionariez.Output_DF_Creation(Dictionariez.Column_Per_Contaminated, SCInputz.N_Iterations) #Progression Dataframe
    
    #Outputs dataframe for each of the sampling steps. 
    df_Output_PH =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe pre harvest
    df_Output_H =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe Harvest
    df_Output_R =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe Receiving
    df_Output_FP =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe Final Product
    df_Output_CS =  Dictionariez.Output_DF_Creation(Dictionariez.Column_Names_Outs, SCInputz.N_Iterations)#Main outputs Dataframe Final Product
    df_Sensitivity = Dictionariez.Output_DF_Creation(Dictionariez.Sensitivity_Analysis_Dic, SCInputz.N_Iterations)
    Series_Final_Conts = []

    
    for  i in range(SCInputz.N_Iterations):
        #i = 1
        Iteration_In = i
        print(Iteration_In,"iteration")
        np.random.seed(i)
        reload(Inputz)
        
        
        #PRCC SENSITYIVITY ANALYSIS SECTION Randomized Inputs only if Sens Analysis is ON. 
        if SCInputz.Sensitivity_Analysis == True:

            reload(ScenCondz)
            #Sampling Contamination

            Sampling_Type_SA = np.random.choice(["PH", "H", "R", "FP", "CS"])
            if Sampling_Type_SA == "PH":
                ScenCondz.PH_Sampling = True
                PH_Type_S = np.random.choice(["4d", "4h", "Int"])
                if PH_Type_S == "4d":
                    ScenCondz.PHS_4d = True
                elif PH_Type_S == "4h":
                    ScenCondz.PHS_4h= True
                elif PH_Type_S == "Int":
                    ScenCondz.PHS_Int= True
            if Sampling_Type_SA == "H":
                ScenCondz.H_Sampling = True
                ScenCondz.HS_Trad = True
            if Sampling_Type_SA == "R":
                ScenCondz.R_Sampling = True
            if Sampling_Type_SA == "FP":
                ScenCondz.FP_Sampling = True
                ScenCondz.FPS_Trad= True 
            if Sampling_Type_SA == "CS":
                ScenCondz.C_Sampling  = True
   
            
            ScenCondz.Contamination_Scenario = np.random.choice([1,2,3])
            #Receiving Factors
            SCInputz.Pre_CoolingYN = np.random.choice([True,False])
            #Processing Factors
            SCInputz.Sanitation_YN = np.random.choice([True,False])
            #7 Washing Random Choice
            SCInputz.Washing_YN = np.random.choice([True,False])
                       
            SCInputz.Spray_WashYN = np.random.choice([True,False])
        
        
        
        ##################BEGGINING OF THE MODEL ###########################
        
        #Adding Contmination to the Field if Contmination Event Occurs Before Pre-Harvest
        
        #STEP 0 CONTAMINATION SCENARIOS  ----------------------------------------------------------------------------------------------------

        
        #Creation of the main model dataframe. 
        df= InFunz.F_InDF(Partition_Units = SCInputz.Partition_Units,
                          Field_Weight = SCInputz.Field_Weight, 
                          slot_number = SCInputz.slot_number)
        #Adding Background contamination to the field. 
        df = ContScen.F_systematic_C(df=df, 
                             Hazard_lvl= 636,
                             No_Cont_Clusters =1,
                             Cluster_Size= 100_000,
                             Partition_Weight = SCInputz.Partition_Weight)
                
        
        #adding contamination to the dataframe is contamination happened before PH. 
        if Inputz.Time_CE_PHS>0: 
            

            #Adding Contamination depending on challenge Systematic Sampling
            
            #Pre-Harvest Contamination Scenarios:
            if ScenCondz.Contamination_Scenario in [1,2,3]:
                df = ContScen.F_systematic_C(df=df, 
                                             Hazard_lvl= Inputz.Hazard_Lvl,
                                             No_Cont_Clusters =Inputz.Cont_Cluster,
                                             Cluster_Size= Inputz.Cluster_Size,
                                             Partition_Weight = SCInputz.Partition_Weight)
                
                
            # Local Outputs: Initial Contamination     
            LV_Initial_CFU= sum(df.CFU) #Initial Contamination 
            Listz.List_Initial_CFU.append(LV_Initial_CFU) #Adding Initial Contamintion to List
            
            #Collection of progression
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                         outputDF = df_Output_Contprog,
                                                         Step_Column = "Contam Event Before PHS", 
                                                         i =Iteration_In )
            #Collection of proportsion contaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column = "PropCont_CE_B_PHS", 
                                                                   i = Iteration_In)
    
        #STEP 1 PREHARVEST ------------------------------------------------------------------------------------------------------------------
        
        #Die-off From Contamination Event to Pre-Havrvest
        #LV_Die_Off_CE_PHS =Funz.F_DieOff_IR_PH(Inputz.Time_CE_PHS,Inputz.Break_Point, Inputz.Dieoff1, Inputz.Dieoff2) #Die off rate from Irrigation to pre harvest sampling, Belias et al. 
        
        
        #DIE-OFF, if contamination before PH Occured.
        LV_Die_off_Total = Funz.F_Simple_DieOff(Inputz.Time_CE_H)
        LV_Die_Off_CE_PHS = ((Inputz.Time_CE_PHS)/Inputz.Time_CE_H)*LV_Die_off_Total

        df = Funz.Applying_dieoff(df=df, Dieoff=LV_Die_Off_CE_PHS ) #Applying Die off to CFU Column in the DF

            
        LO_Cont_B_PH = sum(df.CFU) #Contamination before rejection sampling

        LO_Weight_B_PH = sum(df.Weight)
        

        
        #Contprog Before Pre-Harvest Sampling
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Bef Pre-Harvest Samp", 
                                                     i =Iteration_In )
        
        #PropoContaminated
        df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                               outputDF =df_Output_Propprog, 
                                                               Step_Column =  "PropCont_B_PHS", 
                                                               i = Iteration_In)
        
        Listz.List_BPHS_CFU.append( LO_Cont_B_PH) #List of contamination before sampling

        
        #Sampling at Pre-Harvest
        if ScenCondz.PH_Sampling ==True: #If function to turn off Pre-Harvest Sampling
            if ScenCondz.PHS_Int ==True: #Intense pre harvest sampling
                df = Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_PH_Int, 
                                              NSamp_Unit = SCInputz.n_samples_lot_PH_Int, 
                                              Samp_Size =SCInputz.sample_size_PH_Int, 
                                              Partition_Weight =SCInputz.Partition_Weight, 
                                              NoGrab =SCInputz.No_Grabs_PH_Int)
            elif ScenCondz.PHS_4d==True or ScenCondz.PHS_4h == True :
            #Pre-Harvest Sampling, Traditional
                 df = Funz.F_Sampling_2(df =df,Test_Unit = SCInputz.test_unit_PH, 
                                           NSamp_Unit = SCInputz.n_samples_slot_PH, 
                                           Samp_Size =SCInputz.sample_size_PH, 
                                           Partition_Weight =SCInputz.Partition_Weight, 
                                           NoGrab =SCInputz.No_Grabs_PH)
            
        
        #Filtering out the Rejected lots, Pre-Harvest
        if ScenCondz.PHS_Int ==True: #Rejection intense
           df= Funz.F_Rejection_Rule3(df =df, Test_Unit = SCInputz.RR_PH_Int, limit = SCInputz.Limit_PH)  
        else:  #Rejection normal
            df=Funz.F_Rejection_Rule3(df =df, Test_Unit = SCInputz.RR_PH_Trad, limit = SCInputz.Limit_PH) 
           
        
        #Contprog After Pre-Harvest Sampling
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Aft Pre-Harvest Samp", 
                                                     i =Iteration_In )
        
        #PropoContaminated
        df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                               outputDF =df_Output_Propprog, 
                                                               Step_Column =  "PropCont_A_PHS", 
                                                               i = Iteration_In)
        
                           
        #Outputs Function, Instead of collecting outputs. 
        df_Output_PH = Dictionariez.Output_Collection_Final(df = df, 
                                                            outputDF = df_Output_PH, 
                                                            Step = "PH", 
                                                            Cont_Before = LO_Cont_B_PH, 
                                                            Weight_Before = LO_Weight_B_PH, 
                                                            i = Iteration_In, 
                                                            Niterations = SCInputz.N_Iterations)
        
        
        if Inputz.Time_CE_PHS==0:
        
            
            #Adding Contamination depending on challenge Pre-harvest challenges
            if df["Weight"].sum() != 50:
                if ScenCondz.Contamination_Scenario in [1,2,3]:
                    df = ContScen.F_systematic_C(df=df, 
                                                 Hazard_lvl= Inputz.Hazard_Lvl,
                                                 No_Cont_Clusters =Inputz.Cont_Cluster,
                                                 Cluster_Size= Inputz.Cluster_Size,
                                                 Partition_Weight = SCInputz.Partition_Weight )
                
            # Local Outputs: Initial Contamination     
            LV_Initial_CFU= sum(df.CFU) #Initial Contamination 
            Listz.List_Initial_CFU.append(LV_Initial_CFU) #Adding Initial Contamintion to List
            
            #Contprog Initial
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                         outputDF = df_Output_Contprog,
                                                         Step_Column = "Contam Event After PHS", 
                                                         i =Iteration_In )
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =  "PropCont_CE_A_PHS", 
                                                                   i = Iteration_In)
        
        
        #STEP 2 HARVEST ---------------------------------------------------------------------------------------------------------------------
        
        #Die-off
        LV_Die_Off_PHS_HS=LV_Die_off_Total - LV_Die_Off_CE_PHS
        df = Funz.Applying_dieoff(df=df, Dieoff =LV_Die_Off_PHS_HS ) #Updating Contmination to Show Total DieOff
        
        
        #Harvest Sampling
        if ScenCondz.H_Sampling == True:
            if ScenCondz.HS_Trad==True:
                df = Funz.F_Sampling_2(df =df,
                                             Test_Unit =SCInputz.test_unit_H, 
                                             NSamp_Unit = SCInputz.n_samples_slot_H, 
                                             Samp_Size =SCInputz.sample_size_H, 
                                             Partition_Weight =SCInputz.Partition_Weight, 
                                             NoGrab =SCInputz.No_Grabs_H )
            elif ScenCondz.HS_Agg==True:
               df = Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_H, 
                                               NSamp_Unit = SCInputz.n_samples_slot_H, 
                                               Samp_Size =SCInputz.sample_size_H, 
                                               Partition_Weight =SCInputz.Partition_Weight, 
                                               NoGrab =SCInputz.No_Grabs_H )
        
        

            
        #Before pre harvest sampling
        LO_Cont_B_H = sum(df.CFU) #Contamination before sampling
        LO_Weight_B_H = sum(df.Weight)
        
        ##Listz.List_BHS_CFU.append(LO_Cont_B_H) #List of contaminations before sampling
        
        #Contprog Before Harvest
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Bef Harvest Samp", 
                                                     i =Iteration_In )
        
        #PropoContaminated
        df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                               outputDF =df_Output_Propprog, 
                                                               Step_Column =   "PropCont_B_HS", 
                                                               i = Iteration_In)
        
        
        #Filtering out the Rejected lots, Harvest Sampling
        if ScenCondz.HS_Trad == True:
            df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_H_Trad ,limit = SCInputz.Limit_H) 
        elif ScenCondz.HS_Agg == True:
            df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_H_Agg ,limit = SCInputz.Limit_H) 
        
        
        #Contprog Before Harvest
        df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Aft Harvest Samp", 
                                                     i =Iteration_In )   
        #PropoContaminated
        df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                               outputDF =df_Output_Propprog, 
                                                               Step_Column =   "PropCont_A_PHS", 
                                                               i = Iteration_In)

        
        df_Output_H = Dictionariez.Output_Collection_Final(df = df, 
                                                    outputDF = df_Output_H, 
                                                    Step = "H", 
                                                    Cont_Before = LO_Cont_B_H, 
                                                    Weight_Before = LO_Weight_B_H, 
                                                    i = Iteration_In, 
                                                    Niterations = SCInputz.N_Iterations)
        
        
   
        if (ScenCondz.Field_Pack == False):
            #STEP 3 RECEIVING ---------------------------------------------------------------------------------------------------------------------
            #Time Between Harvest and Pre-Cooling KoseKi. 
            GrowthOutsH_PC = Funz.Growth_Function_Lag(DF =df, 
                                                    Temperature = Inputz.Temperature_H_PreCooling, 
                                                    Time = Inputz.Time_H_PreCooling, 
                                                    Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
            
            df = GrowthOutsH_PC[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsH_PC[1]

            
            LV_bef_precool = sum(df.CFU)
            
            #Pre_Cooling
            if SCInputz.Pre_CoolingYN == True:
                #Aggresive Temperature Change. Reset Lag Time
                Inputz.Lag_Consumed_Prev = 0 #Reseting Lag Consumed.
                #Pre-Cooling Process Cooling.
                GrowthOutsPC = Funz.Growth_Function_Lag(DF =df, 
                                                Temperature = Inputz.Temperature_PreCooling, 
                                                Time = Inputz.Time_PreCooling, 
                                                Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
                df = GrowthOutsPC[0]
                Inputz.Lag_Consumed_Prev = GrowthOutsPC[1]
            LV_af_precool = sum(df.CFU) 
            #print(LV_af_precool-LV_bef_precool, "pre-cool effect")
            
            
            #New Lettuce Temperature approximately 3C if precooling
            
            #Storage at Receiving
            GrowthOutsSto_R = Funz.Growth_Function_Lag(DF =df, 
                                            Temperature = Inputz.Temperature_Storage_R, 
                                            Time = Inputz.Time_Storage_R, 
                                            Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
            
            df = GrowthOutsSto_R[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsSto_R[1]
            
            #print("lag consumed after sto",Inputz.Lag_Consumed_Prev )
            
            #Paletization
            df = Funz.F_Palletization(df=df,
                                      Field_Weight=SCInputz.Field_Weight,
                                      Pallet_Weight=Inputz.Pallet_Weight,
                                      Partition_Weight = SCInputz.Partition_Weight,
                                      )
            
        
            #LV_Time_Agg = LV_Time_Agg + Inputz.Time_H_RS #Cummulative time so far in the process. 
            
            LO_Cont_B_R = sum(df.CFU)
            LO_Weight_B_R = sum(df.Weight)
            
            #Listz.List_BRS_CFU.append(LO_Cont_B_R) #Contamination before receiving sampling
            
            #Contprog Before Receiving
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Bef Receiving Samp", 
                                                     i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =   "PropCont_B_RS", 
                                                                   i = Iteration_In)
            
            if ScenCondz.R_Sampling == True:
                #Sampling at Reception
                df = Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_R, 
                                               NSamp_Unit = SCInputz.n_samples_pallet, 
                                               Samp_Size =SCInputz.sample_size_R, 
                                               Partition_Weight =SCInputz.Partition_Weight, 
                                               NoGrab =SCInputz.No_Grabs_R )
            
            #Rejecting Inidividual pallets if 1 positive
            df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_R_Trad, limit = SCInputz.Limit_R ) 
            
            
            #Contrprog after receiving
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                           outputDF = df_Output_Contprog,
                                           Step_Column =  "After Receiving Samp", 
                                           i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_A_RS", 
                                                                   i = Iteration_In)
            
            
            df_Output_R = Dictionariez.Output_Collection_Final(df = df, 
                                            outputDF = df_Output_R, 
                                            Step = "R", 
                                            Cont_Before = LO_Cont_B_R, 
                                            Weight_Before = LO_Weight_B_R, 
                                            i = Iteration_In, 
                                            Niterations = SCInputz.N_Iterations)

            #STEP 4 Value Addition ---------------------------------------------------------------------------------------------------------------------
            #Splitting pallets into processing lines. 
            gb2 = Funz.F_ProLineSplitting(df =df, Processing_Lines = Inputz.Processing_Lines)

            #Splitting Processing Lines into Mini Batches
            gb2 = Funz.F_Partitioning_ProcLines(gb3 = gb2 , NPartitions = int(Inputz.Pallet_Weight/SCInputz.Partition_Weight))

            #Value Addition Steps
        
            #Cross-Contamination Processing by processing line between 100 lb. batches
            
            #1 Preliminary Spray Wash ------------------------------------------
            
            df_gb2_bsw = (pd.concat(gb2))
            
            #Collecting outputs cont progression 
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bsw,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef SprayWash", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bsw, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_SprayWash", 
                                                                   i = Iteration_In)
            
            if SCInputz.Spray_WashYN == True:
                gb2 = Funz.F_Simple_Reduction_PLines(gb2, Inputz.Harvest_Cspray_red)
            
            #Contamination event, if it happens
            #Pending
                        
            #2 Shredder --------------------------------------------------------
            df_gb2_bs = (pd.concat(gb2))
            
            #Prewashreduction
            TotalRed_PreWash =df_gb2_bs["CFU"].sum()-df_gb2_bsw["CFU"].sum() #-(np.log10(df_gb2_bsw["CFU"].sum()/df_gb2_bs["CFU"].sum()))
            
            #Collecting outputs cont progression 
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bs,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Shredding", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bs, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_Shredding", 
                                                                   i = Iteration_In)
            #Contamination event, if it happens
            #Pending
                
            ShredderOuts = Funz.F_CrossContProLine(gb2 =gb2, 
                                                   Tr_P_S = Inputz.Tr_P_Sh, 
                                                   Tr_S_P= Inputz.Tr_Sh_P,
                                                   Sanitation_Freq_lb = Inputz.Sh_San_freq,
                                                   StepEff = Inputz.Sh_San_Eff , 
                                                   compliance = Inputz.Sh_Compliance)
            gb2 = ShredderOuts[0]
            ShredCont = ShredderOuts[1]  
            
            #3 Conveyor Belt ----------------------------------------------------
            #Contamination before conveyor belt
            df_gb2_bcb = (pd.concat(gb2))
            
            #Collecting outputs cont progression 
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bcb ,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Conveyor Belt", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bcb, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_CBelt", 
                                                                   i = Iteration_In)
            
            #Contamination event, if it happens
            #Pending
                
            CVOuts = Funz.F_CrossContProLine(gb2 = gb2, 
                                             Tr_P_S = Inputz.Tr_P_Cv, 
                                             Tr_S_P = Inputz.Tr_Cv_P,
                                             Sanitation_Freq_lb = Inputz.Cv_San_freq,
                                             StepEff = Inputz.Cv_San_Eff , 
                                             compliance = Inputz.Cv_Compliance)
            gb2 = CVOuts[0]
            CvCont = CVOuts[1]
            
            
            #4Washing:
            df_gb2_bw = (pd.concat(gb2))
            
            #Making it into a dataframe
            gb2 = Funz.F_ProLineSplitting(df =df_gb2_bw, Processing_Lines = Inputz.Processing_Lines)

            #Changing it into washing batches
            gb2 = Funz.F_Partitioning_ProcLines(gb3 = gb2 , NPartitions = int(Inputz.Pallet_Weight/Inputz.Wash_Rate))
            
            
            #Collecting outputs cont progression 
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bw,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Washing", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bw, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_Washing", 
                                                                   i = Iteration_In)
        
            #Contamination event, if it happens
            #Pending
                
            #Washing Step.   
            if SCInputz.Washing_YN == True: 
                gb2 = Funz.F_Washing_ProcLines3(List_GB3 =gb2, Wash_Rate = Inputz.Wash_Rate, Cdf =  Inputz.DF_Chlevels)
            
            #5Shaker Table ----------------------------------------------------
            df_gb2_bst = (pd.concat(gb2))
            
            #Total Reduction Wash
            TotalRed_Wash = df_gb2_bst["CFU"].sum()-df_gb2_bw["CFU"].sum()#np.log10(df_gb2_bw["CFU"].sum()/df_gb2_bst["CFU"].sum())
            
            #Spliting it back into 5 lb chunks. 
            gb2 = Funz.F_ProLineSplitting(df =df_gb2_bst, Processing_Lines = Inputz.Processing_Lines)

            #Splitting Processing Lines into Mini Batches
            gb2 = Funz.F_Partitioning_ProcLines(gb3 = gb2 , NPartitions = int(Inputz.Pallet_Weight/SCInputz.Partition_Weight))
            
            
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bst,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Shaker Table", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bst, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_ST", 
                                                                   i = Iteration_In)
            
            #Contamination event, if it happens
            #Pending
                
            StOuts = Funz.F_CrossContProLine(gb2 =gb2, 
                                             Tr_P_S = Inputz.Tr_P_St, 
                                             Tr_S_P= Inputz.Tr_St_P,
                                             Sanitation_Freq_lb = Inputz.St_San_freq,
                                             StepEff = Inputz.St_San_Eff , 
                                             compliance = Inputz.St_Compliance)
            gb2 = StOuts[0]
            StCont = StOuts[1]
            #print(sum(StCont))
            
            #6 Centrifuge-----------------------------------------------------
            

            
            df_gb2_bcf = (pd.concat(gb2))
            
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df_gb2_bcf,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Centrifuge", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df_gb2_bcf, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_Cent", 
                                                                   i = Iteration_In)
            
            
            
            
            #Contamination event, if it happens
            #Pending
                
            CentrifugeOuts = Funz.F_CrossContProLine(gb2 =gb2, 
                                                     Tr_P_S = Inputz.Tr_P_C, 
                                                     Tr_S_P= Inputz.Tr_C_P,
                                                     Sanitation_Freq_lb = Inputz.C_San_freq,
                                                     StepEff = Inputz.C_San_Eff , 
                                                     compliance = Inputz.C_Compliance)
            gb2 = CentrifugeOuts[0]
            CentrifugeCont = CentrifugeOuts[1]
                
            #Adding Contamination from Scenario to each lot
            
            
            #Joining Data Frames into one again, with contamination from lines. 
            df=(pd.concat(gb2))
            #Outputs after value addition.
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Aft Value Addition", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_A_VA", 
                                                                   i = Iteration_In)
            
            
            
            #LO_Cont_A_VA= sum(df.CFU)
            #Listz.List_AVA_CFU.append(LO_Cont_A_VA)
            
            df['Lot'] =1#Updating lot column to represent finished product lots
            df = Funz.Production_Hours(df = df,FP_Hour = Inputz.FP_Hour)

     
            #Environmental Monitoring Program
            
            
            #STEP 5 PACKING AND MIXING ---------------------------------------------------------------------------------------------------------------------
            
            #Mixing products into one batch
            df = df.reset_index(drop=True) 
            LV_NewPart_Weight  = df.at[0,"Weight"]
            LV_N_Partitions = int(LV_NewPart_Weight/Inputz.Pack_Weight_FP) 
            df = Funz.F_Partitioning(DF=df,
                                     NPartitions= LV_N_Partitions)
            
            if Inputz.N_Lots_FP==2:
                df =Funz.F_Lots_FP(df=df, 
                                   Nolots = 2)
                
            #Dividing the pallets dataframe into different processing lines.  
            gb2 = df.groupby('ProLine')#Creating Listby procesing line
            gb2 =[gb2.get_group(x) for x in gb2.groups] #Creating list of separate dataframe by processing lines
            
            #Contamination event, if it happens
            #Pending
                
            
            df=(pd.concat(gb2))
            #df["Accept"] = True
            #df['PositiveSamples'] = [list() for x in range(len(df.index))]

            df["PositiveSamples"] =""
            df["PositiveSamples"] = [list() for x in range(len(df.index))]
            
            LO_Cont_B_FP = sum(df.CFU) #Total CFU before FP Sampling
            LO_Weight_B_FP = sum(df.Weight)
            
            
            
            #Contamination before sampling
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef Final Prod S", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_FPS", 
                                                                   i = Iteration_In)
            
            
            
            #Sampling Step
            if ScenCondz.FP_Sampling == True:
                if ScenCondz.FPS_Trad ==True:
                    df["FP_Strat"]=np.concatenate([([i]*int(len(df["HourProd"])/SCInputz.N_Strat_FP)) for i in list(range(1,(SCInputz.N_Strat_FP+1)))], axis=0)
                    df =Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_FP, 
                                               NSamp_Unit = SCInputz.n_samples_FP, 
                                               Samp_Size =SCInputz.sample_size_FP, 
                                               Partition_Weight =Inputz.Pack_Weight_FP, 
                                               NoGrab = SCInputz.N_Packages_Samples)
                    
                elif ScenCondz.FPS_Agg ==True:
                    df =Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_FP, 
                                               NSamp_Unit = SCInputz.n_samples_FP, 
                                               Samp_Size =SCInputz.sample_size_FP, 
                                               Partition_Weight =Inputz.Pack_Weight_FP, 
                                               NoGrab = SCInputz.N_Packages_Samples)
            
        
            #Filtering out the Rejected lots, Final product
            #Rejecting Inidividual pallets if 1 positive
            if ScenCondz.FPS_Trad ==True:
                df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_FP_Trad, limit = SCInputz.Limit_FP) 
            elif ScenCondz.FPS_Agg ==True:
                df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_FP_Agg, limit = SCInputz.Limit_FP) 
            
            #collecting outputs
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Final Product Facility", 
                   i =Iteration_In )
            
            #Prop Contaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =   "PropCont_A_FP", 
                                                                   i = Iteration_In)
            
            df_Output_Propprog = Dictionariez.Prop_Collection_WholeField(df = df, 
                                                       outputDF =df_Output_Propprog, 
                                                       Step_Column =   "PropCont_A_FP_Whole", 
                                                       i = Iteration_In)
            #Sum of positive Packages
            df_Output_Propprog.at[i,"TotalCont_A_FP"] = len(df[df.CFU>0])
            
            #Adding final Contamination series
            Series_Final_Conts.append(df.CFU)
            
            
            

            df_Output_FP = Dictionariez.Output_Collection_Final(df = df, 
                                outputDF = df_Output_FP, 
                                Step = "FP", 
                                Cont_Before = LO_Cont_B_FP, 
                                Weight_Before = LO_Weight_B_FP, 
                                i = Iteration_In, 
                                Niterations = SCInputz.N_Iterations)
            
    
        #......End of Field Pack Indentation
        
 
    #STEP 3..Cont Field Pack Lettuce Packing   
        elif ScenCondz.Field_Pack == True:
            
            #Partitioning into Cases
            df = Funz.F_Field_Packing(DF =df, Case_Weight = Inputz.Case_Weight_FieldPack, PartWeight = SCInputz.Partition_Weight)
            
            #Palletizing those cases. 
            
            #Receiving:
            #Post Harvest, Same growth function as other process. 
            GrowthOutsH_PC = Funz.Growth_Function_Lag(DF =df, 
                                                    Temperature = Inputz.Temperature_H_PreCooling, 
                                                    Time = Inputz.Time_H_PreCooling, 
                                                    Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
            
            df = GrowthOutsH_PC[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsH_PC[1]
            
            
            #Pre_Cooling
            #Aggresive Temperature Change. Reset Lag Time
            Inputz.Lag_Consumed_Prev = 0 #Reseting Lag Consumed.
            #New Lettuce Temperature aPproximately 5C
            
            #Storage at Receiving
            GrowthOutsSto_R = Funz.Growth_Function_Lag(DF =df, 
                                            Temperature = Inputz.Temperature_Storage_R, 
                                            Time = Inputz.Time_Storage_R, 
                                            Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
            
            df = GrowthOutsSto_R[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsSto_R[1]
            
            
                        #Contprog Before Receiving
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                                     outputDF = df_Output_Contprog,
                                                     Step_Column = "Bef Receiving Samp", 
                                                     i =Iteration_In )
            
            #Contamination Before: 
            LO_Cont_B_R= sum(df.CFU)
            LO_Weight_B_R= sum(df.Weight)
            
            
            if ScenCondz.R_Sampling == True:
                #Sampling at Reception
                df = Funz.F_Sampling_2(df =df,Test_Unit =SCInputz.test_unit_R_FP, 
                                               NSamp_Unit = SCInputz.n_samples_R_FP, 
                                               Samp_Size =SCInputz.sample_size_R_FP, 
                                               Partition_Weight =Inputz.Case_Weight_FieldPack, 
                                               NoGrab =SCInputz.No_GRabs_R_FP )
            
            #Rejecting Inidividual pallets if 1 positive
            df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_R_FP_Trad, limit = SCInputz.Limit_R_FP) 
            
            
            #Contrprog after receiving
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                                           outputDF = df_Output_Contprog,
                                           Step_Column =  "After Receiving Samp", 
                                           i =Iteration_In )
            

            df_Output_FP = Dictionariez.Output_Collection_Final(df = df, 
                    outputDF = df_Output_FP, 
                    Step = "FP", 
                    Cont_Before = LO_Cont_B_R, 
                    Weight_Before = LO_Weight_B_R, 
                    i = Iteration_In, 
                    Niterations = SCInputz.N_Iterations)
        
        #STEP 6 POST PROCESS STEPS ---------------------------------------------------------------------------------------------------------------------
        #Steps after Final Product
        if ScenCondz.Customer_Added_Steps ==True:
            #Putting Packages into Cases
            if ScenCondz.Field_Pack == False:
                if df["Weight"].sum() != 50:
                    df = Funz.Case_Packaging(df =df,Case_Weight = Inputz.Case_Weight, Pack_Weight = Inputz.Pack_Weight_FP)
            
            #Growth or Die-off during storage post processing storage:
            GrowthOutsPPS = Funz.Growth_Function_Lag(DF =df, 
                                        Temperature = Inputz.Temperature_ColdStorage, 
                                        Time = Inputz.Time_PostPStorage, 
                                        Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
        
            df = GrowthOutsPPS[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsPPS[1]
            
            #Transportation of final product via trucks.            
            
            GrowthOutsPPT = Funz.Growth_Function_Lag(DF =df, 
                            Temperature = Inputz.Transportation_Temp, 
                            Time = Inputz.Trasnportation_Time, 
                            Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
        
            df = GrowthOutsPPT[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsPPT[1]
            
            #Here a step where the lot might get split into different Customers AKA restuarants
            LO_Cont_B_CS = sum(df.CFU) #Total CFU before FP Sampling
            LO_Weight_B_CS = sum(df.Weight)
            
            
            
            #Contamination before sampling
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "Bef CS Samp", 
                   i =Iteration_In )
            
            #PropoContaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =    "PropCont_B_CS", 
                                                                   i = Iteration_In)
            
            
            
            if ScenCondz.C_Sampling == True:
                df["FP_Strat"]=np.concatenate([([i]*int(len(df["HourProd"])/SCInputz.N_Strat_CS)) for i in list(range(1,(SCInputz.N_Strat_CS+1)))], axis=0)
                df =Funz.F_Sampling_2(df =df,
                                        Test_Unit =SCInputz.test_unit_CS, 
                                        NSamp_Unit = SCInputz.n_samples_CS, 
                                        Samp_Size =SCInputz.sample_size_CS, 
                                        Partition_Weight =Inputz.Pack_Weight_FP, 
                                        NoGrab = SCInputz.No_Grabs_CS)
            
            #Filtering out the Rejected lots, Final product
            
            #Rejecting Inidividual pallets if 1 positive
            df = Funz.F_Rejection_Rule3 (df =df, Test_Unit = SCInputz.RR_CS, limit = SCInputz.Limit_CS) 
            
            #collecting outputs
            df_Output_Contprog =  Dictionariez.Output_Collection_Prog(df = df,
                   outputDF = df_Output_Contprog,
                   Step_Column =  "After CS Samp", 
                   i =Iteration_In )
            
            #Prop Contaminated
            df_Output_Propprog = Dictionariez.Pop_Output_Colection(df = df, 
                                                                   outputDF =df_Output_Propprog, 
                                                                   Step_Column =   "PropCont_A_CS", 
                                                                   i = Iteration_In)
            
            df_Output_Propprog = Dictionariez.Prop_Collection_WholeField(df = df, 
                                                       outputDF =df_Output_Propprog, 
                                                       Step_Column =   "PropCont_A_CS_Whole", 
                                                       i = Iteration_In)
            #Sum of positive Packages
            df_Output_Propprog.at[i,"TotalCont_A_CS"] = len(df[df.CFU>0])
            
            #Adding final Contamination series
            #Series_Final_Conts.append(df.CFU)
            
            
            df_Output_CS = Dictionariez.Output_Collection_Final(df = df, 
                                outputDF = df_Output_CS, 
                                Step = "C", 
                                Cont_Before = LO_Cont_B_CS, 
                                Weight_Before = LO_Weight_B_CS, 
                                i = Iteration_In, 
                                Niterations = SCInputz.N_Iterations)
                    
            
            ''' 
            #Storage at customer: 
            GrowthOutsPPCS = Funz.Growth_Function_Lag(DF =df, 
                            Temperature = Inputz.Temperature_PostPCS, 
                            Time = Inputz.Time_PostPCS, 
                            Lag_Consumed_Prev  = Inputz.Lag_Consumed_Prev)
        
            df = GrowthOutsPPCS[0]
            Inputz.Lag_Consumed_Prev = GrowthOutsPPCS[1]
            
            #Washing at consumer: #wash every 2 packs  
            df = Funz.Washing_Batch(df = df, New_water_every_xpacks = 2)
            '''
            
        #SEnsitivity Analaysis Outsputs
        
        if SCInputz.Sensitivity_Analysis==True:
            #Adding Sens outputs
            df_Sensitivity= Dictionariez.Func_LoadInputs(df_Sensitivity,Iteration_In,df,LV_Die_off_Total,TotalRed_PreWash,TotalRed_Wash)
        
            
    #STEP 7: Outputs 
    
    
    if ScenCondz.Customer_Added_Steps ==True:
        df_outputs = pd.concat([df_Output_PH,df_Output_H,df_Output_R, df_Output_FP,df_Output_CS ], axis=1)
    else:
        df_outputs = pd.concat([df_Output_PH,df_Output_H,df_Output_R, df_Output_FP], axis=1)
    
    
    outputs = [df_Output_Contprog, df_outputs,df_Output_Propprog,gb2,df,df_Sensitivity,Series_Final_Conts]

        
    return outputs #Final 
#%%    
#Outs = F_MainLoop()
