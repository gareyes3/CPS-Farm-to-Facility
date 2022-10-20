import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import lognorm
import math
import random
import itertools
import ScenCondz
import matplotlib as plt
import SCInputz
from scipy.stats import beta
from numpy.random import Generator, PCG64
import Inputz
import ContScen
rng = Generator(PCG64())

#Function Source File
#%% Utility Functions
#Normal Truncated function
def Func_NormalTrunc(Min,Max, Mean, SD):
     X= stats.truncnorm((Min - Mean) / SD, (Max - Mean) / SD, loc=Mean, scale=SD)
     Y=(np.float(X.rvs(1)))
     return Y
 
def pert(a, b, c, *, size=1, lamb=4):
    r = c - a
    alpha = 1 + lamb * (b - a) / r
    beta = 1 + lamb * (c - b) / r
    return a + rng.beta(alpha, beta, size=size) * r

def betagen(a,b,mini,maxi):    
    return float(beta.rvs(a , b)* (maxi - mini) + mini)

def transfer_1(val):
    if val<=0.619:
        trans = 1.3*(val-0)/(0.619-0)
    elif (val>=0.619) and (val<=0.946):
        trans = 1.3+((340-1.3)*(val-0.619)/(0.946-0.619))
    else:
       trans = 340+((230000-340)*(val-0.946)/(1-0.946)) 
    return trans


def lognormal_max(mean,sigma):
    while True:
        age = np.random.lognormal(mean=mean,sigma=sigma,size=1)
        if age>=0.84 and age <= 8.4:
            break
    return age

#%% Contmination Scenarios
#Scenario 1: Uniform Contamination
def Irrigation_Water_Cont():
    #Ecoli_soil = np.random.normal(0.549,0.816) #log10 CFU/g
    #P_Ecoli_Water = 0.35
    #PosYN_Ecoli_Water = np.random.binomial(1,P_Ecoli_Water)
    C_Ecoli_Water = (10**np.random.normal(0.604,0.357))/10 #CFU/ml water
    
    Tr_Irr_Water = np.random.uniform(1.8,21.6) #ml/g produce
    #Pr_Irr_Splashing = pert(0.02,0.04,0.05)
    #Pr_Rain_Splashing = 1
    #Soil_Transfer = betagen(0.4,0.8,0.05,16.4)
    #P_Soil_Plant = np.random.uniform(0.35,0.9)
    
    Increase_Irrigation =C_Ecoli_Water*Tr_Irr_Water #CFU/g
    Increase_IrrigationCFU = Increase_Irrigation*(454*100_000)
    return Increase_IrrigationCFU

#Scenario 2: Uniform Contamination
def Irrigation_Soil_Splash():
    Ecoli_soil = np.random.normal(0.549,0.816) #log10 CFU/g
    #Pr_Irr_Splashing = pert(0.02,0.04,0.05)
    Soil_Transfer = betagen(0.4,0.8,0.05,16.4)
    P_Soil_Plant = np.random.uniform(0.35,0.9)
    Cont = (10**Ecoli_soil)*Soil_Transfer*P_Soil_Plant
    Increase_IrrigationCFU = Cont*(454*100_000)
    return int(Increase_IrrigationCFU) 

#Scenario 3
def Feces_Addition(Month): #Wrose case scenrio, High level of contamination.
    Month_Index = Month-1
    Cont_ferral_means =[-29.13,-29.13,-29.13,-29.13,-29.13,-29.13,-1.26,-2.1,-29.13,-2.75,-2.35,-2.13] 
    Cont_ferral_sds = [9.72,9.72,9.72,9.72,9.72,9.72,2.2,2.44,9.72,2.62,2.51,2.46]
    Ferral_Day = 4260 #g of feces per day
    #tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
    Cont_Ferral=lognormal_max(Cont_ferral_means[Month_Index],Cont_ferral_sds[Month_Index])
    Cont_Ferral_CFU=10**Cont_Ferral
    CFU_ferralloc_Soil_CFU = Ferral_Day*Cont_Ferral_CFU #CFU/g
    return float(CFU_ferralloc_Soil_CFU)
    #CFU_ferralloc_CFU*tr #CFU/g trafered to Sublot
    
def Feces_Splash(df,Soil_Slots):
    tr=transfer_1(np.random.uniform(0,1))/(1.29*10**8)
    Cont_Trans = Soil_Slots.apply(lambda x: rng.binomial(x,tr))
    df["CFU"] = df["CFU"]+Cont_Trans
    return df

#Scenario 4: 
def Feces_Runoff(Month): 
    amount_feralprec = 1_000 #g
    amount_cattleprec = 10_000 #g
    Month_Index = Month-1
    Irrigation_m = [0.861194671,0.921357972,0.972496777,0.484588038,0.135367426,0.03824667,0,0.012462398,0.03008165,0.280189085,0.496347228,0]
    Irrigation_m_choice = Irrigation_m[Month_Index]
    Cont_ferral_means =[-29.13,-29.13,-29.13,-29.13,-29.13,-29.13,-1.26,-2.1,-29.13,-2.75,-2.35,-2.13] 
    Cont_ferral_sds = [9.72,9.72,9.72,9.72,9.72,9.72,2.2,2.44,9.72,2.62,2.51,2.46]
    Cont_cattle_means = [-10.62,-7.07,-25.99,3.42,-6.15,-36.01,-1.15,-3.19,-3.62,-1.68,-6.75,-8.68]
    Cont_cattle_sds= [5.12,4.16,9.25,3.18,3.91,11.95,2.57,3.11,3.23,2.71,4.07,4.56]
    
    Cont_ferral = lognormal_max(Cont_ferral_means[Month_Index],Cont_ferral_sds[Month_Index])
    Cont_ferral_CFU = 10**Cont_ferral
    Cont_cattle =lognormal_max(Cont_cattle_means[Month_Index],Cont_cattle_sds[Month_Index])
    Cont_cattle_CFU = 10**Cont_cattle
    
    amount_cattleprec_Final=amount_cattleprec*Irrigation_m_choice
    amount_ferralprec_Final = amount_feralprec*Irrigation_m_choice
    
    CFU_ferralprec = amount_ferralprec_Final*Cont_ferral_CFU
    CFU_cattleprec = amount_cattleprec_Final*Cont_cattle_CFU
    
    Total_Soil_Add  =  CFU_ferralprec + CFU_cattleprec
    return float(Total_Soil_Add)

#%% Contamination Function
def Cont_Ini_PHS(df):
    TotalTime = 0 #Time Tracker
    Total_CFU_v = []
    #Contamination part 1: From Beginning-to-PHS
    for i in range(Inputz.Time_I_PHS_FullD): #Loop to repeat over the full days before PHS
        #Determining Contamination Scenario
        
        #Scenario #1
        if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
            #Calculating total CFU in Field for irrigation fields
            if TotalTime in Inputz.Final_Irrigation_Days:
                Cont = Irrigation_Water_Cont()
                #Adding Uniform Contmaination to the field
                df = ContScen.F_systematic_C(df =df, 
                              Hazard_lvl = Cont,
                              No_Cont_Clusters=1,
                              Cluster_Size=SCInputz.Field_Weight, 
                              Partition_Weight= SCInputz.Partition_Weight)
                

        #Scenario #2
        if Inputz.Scenario_no == 2:
            if TotalTime in Inputz.Final_Irrigation_Days:
                Cont = Irrigation_Soil_Splash()
                #Adding Uniform Contmaination to the field
                df = ContScen.F_systematic_C(df =df, 
                              Hazard_lvl = Cont,
                              No_Cont_Clusters=1,
                              Cluster_Size=SCInputz.Field_Weight, 
                              Partition_Weight= SCInputz.Partition_Weight)
        
        #Scenario 3
        if Inputz.Scenario_no == 3:
            Soil_Slots = df["CFU"]
            #Contaminating 10 random partitions
            Index_Cont=Soil_Slots.sample(n=10).index
            Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
            Soil_Slots[Index_Cont]=Updates_values[0]
            Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
            if TotalTime in Inputz.Final_Irrigation_Days:
                df = Feces_Splash(df,Soil_Slots)
        #Scenario 4
        if Inputz.Scenario_no == 4:
            Soil_Slots = df["CFU"]
            #Contaminating 10 random partitions
            Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
            Soil_Slots=Updates_values[0]
            Soil_Slots=pd.Series(Soil_Slots)
            Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
            if TotalTime in Inputz.Final_Irrigation_Days:
                df = Feces_Splash(df,Soil_Slots)        
            
        TotalCFU = df["CFU"].sum()
        Total_CFU_v.append(TotalCFU)
        
        #Applying dieoff to the pants every day.
        df = Applying_dieoff(df=df, Dieoff=F_Simple_DieOff(1))
        #Adding time to the time counter.
        TotalTime = TotalTime+1
    
    
        
    if Inputz.Time_I_PHS_PartD >0: #Only if the partial day is greater than 0
        #Part 2: Partial Day, deciding if contamination occurs in the first or second part of the day
        if rng.uniform(0,1)<Inputz.Time_I_PHS_PartD: 
            Happens_YN_1part=1 
        else: 
            Happens_YN_1part=0
        
        if Happens_YN_1part==1:
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)

    
            if Inputz.Scenario_no == 2: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
            
            if Inputz.Scenario_no == 3:  
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_I_PHS_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
            
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_I_PHS_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)  
                        
    
            
        TotalCFU = df["CFU"].sum()  
        Total_CFU_v.append(TotalCFU)
        #Dieoff for first part of the day:
        df = Applying_dieoff(df=df, Dieoff=F_Simple_DieOff(Inputz.Time_I_PHS_PartD))
        TotalTime = TotalTime+Inputz.Time_I_PHS_PartD
        
    if Inputz.Time_I_PHS_PartD ==0:
        Happens_YN_1part = 0
        
        
    return [df, Total_CFU_v, TotalTime, Happens_YN_1part]



def Cont_PHS_End(df, Happens_YN_1part, TotalTime, Total_CFU_v ):
    #Part 2 After Sampling
    if Inputz.Time_PHS_H_PartD > 0: #If Time between PHS-H is less than 1
        if Happens_YN_1part==0:
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                if (TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
                    
            if Inputz.Scenario_no == 2: #scenario 1: Irrigation water contamination
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
                    
            if Inputz.Scenario_no == 3:  
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_PHS_H_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
                
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**(-0.1744*Inputz.Time_PHS_H_PartD))) #applying daily die-off.            
                if math.ceil(TotalTime) in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots) 
            
        TotalCFU = df["CFU"].sum()  
        Total_CFU_v.append(TotalCFU)
        df = Applying_dieoff(df=df, Dieoff=F_Simple_DieOff(Inputz.Time_PHS_H_PartD))
        TotalTime = TotalTime+Inputz.Time_PHS_H_PartD
    
    
    for i in range(Inputz.Time_PHS_H_FullD):
            #Determining Contamination Scenario
            
            #Scenario #1
            if Inputz.Scenario_no == 1: #scenario 1: Irrigation water contamination
                #Calculating total CFU in Field for irrigation fields
                if TotalTime in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Water_Cont()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
                    
    
            #Scenario #2
            if Inputz.Scenario_no == 2:
                if TotalTime in Inputz.Final_Irrigation_Days:
                    Cont = Irrigation_Soil_Splash()
                    #Adding Uniform Contmaination to the field
                    df = ContScen.F_systematic_C(df =df, 
                                  Hazard_lvl = Cont,
                                  No_Cont_Clusters=1,
                                  Cluster_Size=SCInputz.Field_Weight, 
                                  Partition_Weight= SCInputz.Partition_Weight)
            
            #Scenario 3
            if Inputz.Scenario_no == 3:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Index_Cont=Soil_Slots.sample(n=10).index
                Updates_values=Soil_Slots[Index_Cont].values+ rng.multinomial(Feces_Addition(Inputz.Month_choice),[1/10]*10,1)
                Soil_Slots[Index_Cont]=Updates_values[0]
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
                if TotalTime in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)
            #Scenario 4
            if Inputz.Scenario_no == 4:
                Soil_Slots = df["CFU"]
                #Contaminating 10 random partitions
                Updates_values=Soil_Slots.values+ rng.multinomial(Feces_Runoff(Inputz.Month_choice),[1/2000]*2000,1)
                Soil_Slots=Updates_values[0]
                Soil_Slots=pd.Series(Soil_Slots)
                Soil_Slots = Soil_Slots.apply(lambda x: rng.binomial(x,10**-0.1744)) #applying daily die-off.            
                if TotalTime in Inputz.Final_Irrigation_Days:
                    df = Feces_Splash(df,Soil_Slots)        
                
            TotalCFU = df["CFU"].sum()
            Total_CFU_v.append(TotalCFU)
            #Applying dieoff to the pants every day.
            df = Applying_dieoff(df=df, Dieoff=F_Simple_DieOff(1))
            #Adding time to the time counter.
            TotalTime = TotalTime+1
    return [df, Total_CFU_v, TotalTime, Happens_YN_1part]

#%% Die off functions    

#Die-off model

def F_DieOff1():
    Die_off_rate1=Func_NormalTrunc(-16.52,-0.47, -7.07,3.41)
    return Die_off_rate1

def F_DieOff2():
    Die_off_rate2=Func_NormalTrunc(-1.94,3.04, -0.24,0.70)
    return Die_off_rate2

#Die -off from Irrigation to Pre-Harvest
def F_DieOff_IR_PH(Time, Break_Point,Dieoff1, Dieoff2):
    TimeD = Time
    if TimeD < Break_Point: 
        T_Die_off= Dieoff1*TimeD
    elif TimeD >= Break_Point:
        Seg1T = TimeD-Break_Point
        T_Die_off1=Dieoff1*Seg1T
        Seg2T = TimeD - Seg1T
        T_Die_off2= Dieoff2*Seg2T
        T_Die_off = T_Die_off1+T_Die_off2
    return T_Die_off

#Die-off from Pre-Harvest Sampling to Harvest sampling
def F_DieOff_PHS_HS(Time,Time_Agg,Break_Point ,Dieoff1, Dieoff2 ):
    if Time_Agg < Break_Point: 
        TimeLeft = Break_Point-Time_Agg
        if Time < TimeLeft:
            T_Die_off = Dieoff1*Time
        elif Time >=TimeLeft:
            Seg1T = Time-TimeLeft
            T_Die_off1=Dieoff1*Seg1T
            Seg2T = Time-Seg1T
            T_Die_off2=Dieoff2*Seg2T
            T_Die_off = T_Die_off1+T_Die_off2
    elif Time_Agg>=Break_Point:
        T_Die_off=Dieoff2*Time
    return T_Die_off


def F_Simple_DieOff (Time): 
    #Reduction = -((Time/(2.45/24))**0.3)
    dieoff= Func_NormalTrunc(-1.04,-0.33, -0.77,0.21) #Belias linear die-off. 
    Reduction  = dieoff*Time
    return Reduction

def Applying_dieoff (df,Dieoff):
    vector = np.array(df["CFU"])
    new_vector=[]
    for i in vector: 
        CFU_1 = i
        new_cont=np.random.poisson(10**Dieoff, CFU_1).sum()
        #new_cont = np.random.binomial(CFU_1,10**Dieoff)
        new_vector.append(new_cont)
    df["CFU"] = new_vector
    return df

def F_Simple_Reduction(df, Reduction):
    vector = np.array(df["CFU"])
    new_vector=[]
    for i in vector: 
        CFU_1 = i
        new_cont = np.random.binomial(CFU_1,10**-Reduction)
        #new_cont=rng.poisson(10**-Reduction, CFU_1).sum()
        new_vector.append(new_cont)
    df["CFU"] = new_vector
    return df

def F_Simple_Reduction_PLines(gb2, Reduction):
    for j in gb2: 
        vector = np.array(j["CFU"])
        new_vector=[]
        for i in vector: 
            CFU_1 = i
            new_cont = np.random.binomial(CFU_1,10**-Reduction)
            #new_cont=rng.poisson(10**-Reduction, CFU_1).sum()
            new_vector.append(new_cont)
        j["CFU"] = new_vector
    return gb2



#%% Growth  or Reduction Models
#Cold Storage growth Model


def F_Growth(DF,Temperature, TimeD ):
    Parition_Weight_g = 50*454
    CFUs = DF["CFU"]
    b = 0.023
    Tmin  = 1.335-5.766 *b
    New_CFUs=[]
    for i in CFUs:
        CFUs_g = i/Parition_Weight_g #CFU/g
        Total_CFU = i
        if CFUs_g < (10**7): #7 log max density
            DieoffRate = rng.triangular(0.0035, 0.013,0.040)/2.303 #log CFU /g h
            TotalGrowthRate = (b*(Temperature - Tmin))**(2)/(2.303) #log CFU /g h
            if Temperature >5:
                Growth  = 1
            else:
                Growth = 0
            if Growth == 1:
                TotalGrowth = (TotalGrowthRate*TimeD)
            else:
                TotalGrowth  = -DieoffRate * TimeD
            '''
            if TotalGrowth>=0:
                GrowthCeil = math.ceil(TotalGrowth)
                Difference = TotalGrowth-GrowthCeil
                MaxCont = Total_CFU*10**GrowthCeil
                Updated_CFUs = rng.binomial(MaxCont,10**Difference)
            else:
                Updated_CFUs= rng.binomial(Total_CFU,10**TotalGrowth)
            '''
            Updated_CFUs=rng.poisson(10**TotalGrowth, Total_CFU).sum()
        else:
            Updated_CFUs = Total_CFU
        New_CFUs.append(Updated_CFUs)
    DF["CFU"] =New_CFUs 
    return DF



#Calculate Lag time at given temperature
def Growth_Function_Lag(DF, Temperature,Time,Lag_Consumed_Prev):
    if Temperature > 5:
        Lag_Time = 7544*(Temperature**-3.11) #laf frmula
        Proportion_Lag_Consumed = Time/Lag_Time
        Cummulative_Lag_Consumed = Lag_Consumed_Prev + Proportion_Lag_Consumed
        if Cummulative_Lag_Consumed < 1: 
            df2 = DF
        if Cummulative_Lag_Consumed >1:
            #time not in lag phase
            if Lag_Consumed_Prev < 1:
                PropNotLag =(((Cummulative_Lag_Consumed - 1))/Cummulative_Lag_Consumed)
                Growth_Time = Time*PropNotLag
                df2 = F_Growth(DF = DF,Temperature = Temperature,TimeD = Growth_Time)
            elif Lag_Consumed_Prev >1:
                Growth_Time = Time
                df2=F_Growth(DF =DF,Temperature = Temperature,TimeD = Growth_Time)
        Lag_Consumed_Prev = Cummulative_Lag_Consumed
    elif Temperature <=5:
        df2= F_Growth(DF =DF, Temperature = Temperature, TimeD = Time)
    outputs = [df2,Lag_Consumed_Prev]
    return outputs


#%% Contamination Functions

#Calculation of E.coli in Water
def F_Ecoli_Water():   
    Cw = rng.uniform(1,235)
    Rw = 10**Func_NormalTrunc(Min= -100000, Max = 0, Mean=-1.9, SD = 0.6)
    W = Func_NormalTrunc(Min = 0, Max = 100000, Mean = 0.108, SD = 0.019)
    Ci = ((Cw/100)*Rw*W)*(454*SCInputz.Field_Weight)
    return Ci

def F_HarvestingCont ():
    Cs = 10**Func_NormalTrunc(Min=0, Max=3.67, Mean = 0.928, SD=1.11) #Soil, E.coli conc
    Rs = 10**Func_NormalTrunc(Min=-5, Max=0, Mean = -1.9, SD=0.6) #prob
    M = 10.22 #g soil in blades
    Nb = Cs *Rs*M #E coli cells in Blade CFU
    Rt1 = 0.0013 #Trasnfer Rate from soil to Lettuce
    Nh1 = Nb*Rt1 #Total E coli from harvesting blades to lettuce
    return Nh1

def F_InitialCont():
    #Using base pert distribution
    #Calculation of total CFUs
    Cont_CFU_g = float(pert(0,0,634))
    g_field = SCInputz.Field_Weight*454 #454 g in 1 lb. 
    Final_Cont = int(Cont_CFU_g*g_field)
    return Final_Cont

'''
Cs = 10**Func_NormalTrunc(Min= 0, Max = 3.67, Mean=0.928, SD =1.1)
Rs = 10**Func_NormalTrunc(Min= -100000, Max = 0, Mean=-1.9, SD = 0.6)
M = 10.22
Cs*Rs*M
 '''   
    
#%% Sampling Functions

#New Sampling Function
'''
#Outdates, Replaced with more efficient way to do the grabs iteration.
def F_Sampling_2 (df, Test_Unit, NSamp_Unit, Samp_Size, Partition_Weight, NoGrab):
    Unique_TestUnit = list(df[Test_Unit].unique())
    Grab_Weight = Samp_Size/NoGrab
    for i in (Unique_TestUnit): #From sublot 1 to sublot n (same for pallet,lot,case etc)
        for l in range (1, NSamp_Unit+1): #Number of samples per sublot or lot or pallet.
            for j in range(NoGrab):
                Sampled_Grab =df[df[Test_Unit] == i].sample(1, replace= True)
                Index = Sampled_Grab.index
                Index = Index[0]
                CFU = Sampled_Grab["CFU"]
                CFU_grab = CFU*(Grab_Weight/(Partition_Weight*454))
                P_Detection=1-math.exp(-CFU_grab)
                RandomUnif = random.uniform(0,1)
                if RandomUnif < P_Detection:
                    df.at[Index, 'PositiveSamples']. append(l)
    return (df)
'''

def F_Sampling_2 (df, Test_Unit, NSamp_Unit, Samp_Size, Partition_Weight, NoGrab):
    Unique_TestUnit = list(df[Test_Unit].unique())
    Grab_Weight = Samp_Size/NoGrab
    for i in (Unique_TestUnit): #From sublot 1 to sublot n (same for pallet,lot,case etc)
        for l in range (1, NSamp_Unit+1): #Number of samples per sublot or lot or pallet.
            for j in range(NoGrab):
                CFU_hh=np.array(df["CFU"])
                List_Random=random.choice(list(enumerate(CFU_hh)))
                CFU = List_Random[1]
                Index = List_Random[0]
                CFU_grab = CFU*(Grab_Weight/(Partition_Weight*454))
                P_Detection=1-math.exp(-CFU_grab)
                RandomUnif = random.uniform(0,1)
                if RandomUnif < P_Detection:
                    df.at[Index, 'PositiveSamples'].append(l)
    return (df)



#New Rejection Function
def F_Rejection_Rule2 (df, Test_Unit, limit):
    #Test_Unit = "Lot" or "Sublot"
    Listpositive = []
    for i, row in df.iterrows():
        Positives = len(set(df.at[i, "PositiveSamples"]))
        Listpositive.append(Positives)
    df.Positives =Listpositive
    Positives = df[df["Positives"]> limit]
    Unique_TestUnit=list(df[Test_Unit].unique())
    Unique_Positives = list(Positives[Test_Unit].unique())
    df.Positives = ""
    df.PositiveSamples = [list() for x in range(len(df.index))]
    if set(Unique_TestUnit)<= set(Unique_Positives):
        df_Blank = df.iloc[[0]]
        df_Blank.loc[:, ['CFU']] = 0
        df_Blank.loc[:, ['Weight']] = 1000
        df_Blank.loc[:, ['Accept']] = "All Rej"
        df = df_Blank
    else:
        df = df[~df[Test_Unit].isin(Unique_Positives)]
    return df


def F_Rejection_Rule3 (df, Test_Unit, limit):
    Unique_Test_Unit =list(df[Test_Unit].unique())
    Reject = []
    for  i in Unique_Test_Unit:
        df_Subset = df[df[Test_Unit] == i]
        List_of_grabs = df_Subset['PositiveSamples'].tolist()
        flat_list = [item for sublist in  List_of_grabs for item in sublist]
        Unique_Positives =list(np.unique(flat_list))
        if len(Unique_Positives)>limit:
            Reject.append(i)
    df.PositiveSamples = [list() for x in range(len(df.index))] #this is in case everything gets rejected
    if set(Unique_Test_Unit)<= set(Reject):
        df_Blank = df.iloc[[0]].copy()
        df_Blank.loc[:, ['CFU']] = 0
        df_Blank.loc[:, ['Weight']] = SCInputz.Partition_Weight
        df_Blank.loc[:, ['Accept']] = "All Rej"
        df = df_Blank
    else:
        df = df[~df[Test_Unit].isin(Reject)]
    return df
        


#Sampling Sublots
#df =dateframe,
#NoSampleLot = Number of samples per lot
#Sample_size = Composite sample weight
#Cluster_Unit_weight 
#Limit= Limit CFU
#Grabs = Total grabs per sublot
def F_Sampling(df,Test_Unit, NSamp_Unit, Samp_Size, Partition_Weight, Limit, NoGrab ):
    #Df where contaminations are located
    #Test_Unit: Testing lot, as a whole, or sublots
    #N Samples/ Units: How many samples are taken per Test_Unit
    #Samp_Size: Sample Size Composite sample size
    #Partition_Weight: Weight of the partition
    #Limit : m maximum level of positive samples
    #No Grab: N60, grabs. 
    Results=[]
    Unique_TestUnit=list(df[Test_Unit].unique())
    for i in (Unique_TestUnit):
        Reject_Lis=[]
        for l in range (NSamp_Unit):
            Sampled_Grabs =df[df[Test_Unit] == i].sample(NoGrab, replace= True)
            Sampled_Grabs =list(Sampled_Grabs.CFU)
            Grab_Weight = Samp_Size/NoGrab
            Detected = []
            for j in Sampled_Grabs: 
                CFU_grab = j*(Grab_Weight/(Partition_Weight*454))
                P_Detection=1-math.exp(-CFU_grab)
                if random.uniform(0,1)<P_Detection:
                    Reject_YN=1
                else:
                    Reject_YN=0
                Detected.append(Reject_YN)
                if sum(Detected)>0:
                    Detected_YN = 1
                elif sum(Detected) ==0:
                    Detected_YN =0
            Reject_Lis.append(Detected_YN)
        a=sum(Reject_Lis)
        if a > Limit:
            AR= False
        else:
            AR= True
        Results.append(AR)
    data1 =  {Test_Unit: Unique_TestUnit,
           'Accept_Reject': Results}
    dT = pd.DataFrame(data1)
    dT= dT.loc[dT['Accept_Reject'] == False]
    ListR= list(dT[Test_Unit])
    return(ListR)

def F_SamplingFProd (df, Test_Unit, N_SampPacks, Grab_Weight):
    Results=[]
    Clust_Weight = int(df.loc[1,"Weight"]*454)
    Sampled_Packs = list(df[Test_Unit].sample(N_SampPacks))
    for i in Sampled_Packs:
        CFU = df.loc[df['PackNo'] == i, 'CFU'].values[0]
        CFU_grab = CFU*(Grab_Weight/(Clust_Weight*454))
        P_Detection=1-math.exp(-CFU_grab)
        if random.uniform(0,1)<P_Detection:
            Reject_YN=1
        else:
            Reject_YN=0
        Results.append(Reject_YN)
    return Results


def F_Rejection_Rule (df, LL_Rej_Lots, Test_Unit):
    #Test_Unit = "Lot" or "Sublot"
    Unique_TestUnit=list(df[Test_Unit].unique())
    if set(Unique_TestUnit)<= set(LL_Rej_Lots):
        df_Blank = df.iloc[[0]]
        df_Blank.loc[:, ['CFU']] = 0
        df_Blank.loc[:, ['Weight']] = 0
        df = df_Blank
    else:
        df = df[~df[Test_Unit].isin(LL_Rej_Lots)]
    return df

#%% Partitioning and Mixing Functions
def F_Palletization (df, Field_Weight,Pallet_Weight, Partition_Weight):
    Partitions_Per_Pallet =  int(Pallet_Weight/Partition_Weight)
    Pallet_Field = int(Field_Weight/Pallet_Weight)
    Pallet_Pattern = [i for i in range(1, Pallet_Field+1) for _ in range(int(Partitions_Per_Pallet))]
    Crop_No = len(df.index)
    Pallet_Pattern=Pallet_Pattern[:Crop_No]
    df['PalletNo'] = Pallet_Pattern
    df = df[['Lot', 'Sublot','PalletNo','PartitionID','CFU','PositiveSamples','Accept', 'Weight']]
    return df



def F_ProLineSplitting(df, Processing_Lines,): #
    df2=df.groupby(['PalletNo'], as_index =False)[["CFU", "Weight"]].sum()
    #Splitting Pallets into processing lines. Faccept
    N_Pallets = len(df2.index)
    num, div = N_Pallets, Processing_Lines #Getting list of pallets per line
    N_Divs =  ([num // div + (1 if x < num % div else 0)  for x in range (div)])
    N_Lines = list(range(1,Processing_Lines+1))
    L_ProLine =list(itertools.chain(*(itertools.repeat(elem, n) for elem, n in zip(N_Lines, N_Divs))))
    df2["ProLine"] = L_ProLine
    #Dividing the pallets dataframe into different processing lines.  
    gb = df2.groupby('ProLine')#Creating Listby procesing line
    gb2 =[gb.get_group(x) for x in gb.groups] #Creating list of separate dataframe by processing lines
    return gb2



def F_CrossContProLine (gb2, Tr_P_S, Tr_S_P,Sanitation_Freq_lb = 0, StepEff = 0 , compliance = 0 ):
    ContS_L=[]
    for j in gb2:
        rateweight = j.iloc[0]['Weight']
        every_x_many = int(Sanitation_Freq_lb/rateweight)
        ContS=0
        if every_x_many > 0:
            Cleaning_steps = np.arange(0, len(j) , every_x_many )
        for row in j.itertuples():
            i  = row[0]
            if random.uniform(0,1)<compliance:
                if every_x_many > 0:
                    if i in Cleaning_steps:
                        if ContS>0:
                            ContS = ContS*(10**StepEff) 
            ContP = j.CFU[i] #Contamination product
            TotTr_P_S= rng.binomial(ContP,Tr_P_S) #Transfer from Product to Surfaces
            TotTr_S_P = rng.binomial(ContS,Tr_S_P) #Trasnfer from Surfaves to product
            ContPNew = ContP-TotTr_P_S+TotTr_S_P #New Contmination on Product
            ContS=ContS+TotTr_P_S-TotTr_S_P #Remiining Contamination in Surface for upcoming batches
            j.at[i,("CFU")]=ContPNew #Updating the Contamination in the Data Frame
        ContS_L.append(ContS)
    Outputs = [gb2,ContS_L]
    return Outputs


def F_CrossContProLine2 (gb2, Tr_P_S, Tr_S_P, Sanitation_Freq_lb = 0, StepEff = 0 , compliance = 0 ):
    ContS_L=[]
    for j in gb2:
        rateweight = j.iloc[0]['Weight']
        every_x_many = int(Sanitation_Freq_lb/rateweight)
        ContS=0
        vectorCFU = j.CFU
        newvector=[]
        if every_x_many > 0:
            Cleaning_steps = np.arange(0, len(vectorCFU) , every_x_many )
        for i in range(len(vectorCFU)):
            if random.uniform(0,1)<compliance:
                if every_x_many > 0:
                    if i in Cleaning_steps:
                        ContS = ContS*10**StepEff
                        print ("cleaned")
            ContP = vectorCFU[i] #Contamination product
            TotTr_P_S= rng.binomial(ContP,Tr_P_S) #Transfer from Product to Surfaces
            TotTr_S_P = rng.binomial(ContS,Tr_S_P) #Trasnfer from Surfaves to product
            ContPNew = ContP-TotTr_P_S+TotTr_S_P #New Contmination on Product
            ContS=ContS+TotTr_P_S-TotTr_S_P #Remiining Contamination in Surface for upcoming batches
            newvector.append(ContPNew)
        j["CFU"] = newvector
        ContS_L.append(ContS)
    Outputs = [gb2,ContS_L]
    return Outputs



#Paritioning Function
def F_Partitioning(DF,NPartitions):
    if ScenCondz.Field_Pack==False:
        AllParts_Cont = []
        for row in DF.itertuples():
            i = row[0]
            Cont = DF.at[i,'CFU']
            PartCont=rng.multinomial(Cont,[1/NPartitions]*NPartitions,size=1)
            PartCont = PartCont[0]
            AllParts_Cont.append(PartCont)
        b_flat = [j for i in AllParts_Cont for j in i]
        newdf = pd.concat([DF]*NPartitions,axis=0)
        newdf=newdf.sort_values(by=['PalletNo'])
        #Pallet_List=(list(range(1,NPartitions+1)))
        newdf["PackNo"] =list(range(1,len(newdf.index)+1))#np.tile(Pallet_List, len(newdf)//NPartitions)
        newdf = newdf.reset_index(drop=True)  
        newdf.Weight=newdf.Weight/NPartitions
        newdf.CFU = b_flat
        newdf["Sublot"] = 1
        newdf = newdf[['PalletNo','HourProd','PackNo','CFU', 'Weight', 'Sublot','ProLine','Lot']]
    elif ScenCondz.Field_Pack == True:
        AllParts_Cont = []
        for row in DF.itertuples():
            i = row[0]
            PartCont=rng.multinomial(Cont,[1/NPartitions]*NPartitions,size=1)
            PartCont = PartCont[0]
            AllParts_Cont.append(PartCont)
        b_flat = [j for i in AllParts_Cont for j in i]
        newdf = pd.concat([DF]*NPartitions,axis=0)
        #Pallet_List=(list(range(1,NPartitions+1)))
        newdf["CaseNo"] =list(range(1,len(newdf.index)+1))#np.tile(Pallet_List, len(newdf)//NPartitions)
        newdf = newdf.reset_index(drop=True)  
        newdf.Weight=newdf.Weight/NPartitions
        newdf.CFU = b_flat
        newdf["Sublot"] = 1
        newdf = newdf[['CaseNo','CFU', 'Weight', 'Sublot','Lot']]
    return newdf

def F_Field_Packing(DF, Case_Weight, PartWeight):
    NPartitions = int(PartWeight/Case_Weight) 
    AllParts_Cont = []
    for row in DF.itertuples():
        i = row[0]
        Cont = DF.at[i,'CFU']
        PartCont=rng.multinomial(Cont,[1/NPartitions]*NPartitions,size=1)
        PartCont = PartCont[0]
        AllParts_Cont.append(PartCont)
    b_flat = [j for i in AllParts_Cont for j in i]
    newdf = pd.concat([DF]*NPartitions,axis=0)
    newdf=newdf.sort_values(by=['Sublot'])
    #Pallet_List=(list(range(1,NPartitions+1)))
    newdf["CaseNo"] =list(range(1,len(newdf.index)+1))#np.tile(Pallet_List, len(newdf)//NPartitions)
    newdf = newdf.reset_index(drop=True)  
    newdf.Weight=newdf.Weight/NPartitions
    newdf.CFU = b_flat
    newdf = newdf[['CaseNo','CFU', 'Weight', 'Sublot','Lot']]
    return newdf
    

def F_Lots_FP(df, Nolots):
    l = len(df.index) // Nolots 
    df.loc[:l - 1, "Sublot"] = Nolots-1
    df.loc[l:, "Sublot"] = Nolots
    return df

def F_Mixing(DF):
    CFU_Summation = sum(DF.CFU)
    gram_Summation = sum(DF.Weight)*454
    Cont = CFU_Summation/gram_Summation
    ArrayUnique= pd.unique(DF['Sublot'])
    data1 = {'Sublot': [ArrayUnique],
       'Lot': 1,
       'Cont':Cont,
       'CFU':CFU_Summation,
       'Accept': True,
       'Weight': sum(DF.Weight)}
    df1 = pd.DataFrame(data1)  
    return df1


def parts(a, b):
    q, r = divmod(a, b)
    return [q + 1] * r + [q] * (b - r)


def F_Partitioning2(DF, Partition_Weight):
    LWeights = []
    LXX_2 = []
    for row in DF.itertuples():
        i = row[0]
        Weight= int(DF.at[i,'Weight'])
        xx_2=int(Weight//Partition_Weight)
        LXX_2.append(xx_2)
        Weight2 = parts(Weight,xx_2)
        LWeights.append(Weight2)  
    LWeightsFlat = [item for sublist in  LWeights for item in sublist]
    newDF= DF.loc[DF.index.repeat(LXX_2)]
    newDF["Weight"] = LWeightsFlat
    AllParts_Cont = []
    b_flat=[]
    DF['Parts'] =LXX_2
    for row in DF.itertuples():
        i = row[0]
        Cont = DF.at[i,'CFU']
        Parts = int(DF.at[i,'Parts'])
        PartCont=rng.multinomial(Cont,[1/Parts]*Parts, size =1)
        PartCont = PartCont[0]
        AllParts_Cont.append(PartCont)
    b_flat = [j for i in AllParts_Cont for j in i]
    newDF.CFU = b_flat
    return newDF


#%%
def F_Partitioning_W(DF,NPartitions):
    AllParts_Cont = []
    for row in DF.itertuples():
        i = row[0]
        Cont = DF.at[i,'CFU']
        PartCont=rng.multinomial(Cont,[1/NPartitions]*NPartitions,size=1)
        PartCont = PartCont[0]
        AllParts_Cont.append(PartCont)
    b_flat = [j for i in AllParts_Cont for j in i]
    newdf = pd.concat([DF]*NPartitions,axis=0)
    newdf=newdf.sort_values(by=['PalletNo'])
    #Pallet_List=(list(range(1,NPartitions+1)))
    newdf["Part"] =list(range(1,len(newdf.index)+1))#np.tile(Pallet_List, len(newdf)//NPartitions)
    newdf = newdf.reset_index(drop=True)  
    newdf.Weight=newdf.Weight/NPartitions
    newdf.CFU = b_flat
    return newdf



#Washing Chloride:
    #This function creates the CDF of chlorine levels (mg/L) ppm, per unit time, in this case every 0.1 minutes. 
    #Timewash is the time washed per minutes default is set to 300 min. 
    #Chlorine levels increase every 12 min
def F_Chloride_lvl (Time_Wash, Treatment):
    #Function Inputs. 
    #Changing times to 0.1 increments.
    Times = np.arange(0, Time_Wash+0.01, 0.01).tolist()
    Times = [round(num, 2) for num in Times]
    #Addition Rates
    if Treatment ==1:
        r1= 12.75 #(mg/(ml/min**2))
        r2 = 7.47 #(mg/(ml/min**2))
        r3 = 5.56 #(mg/(ml/min**2))
    elif Treatment == 0:
        r1= 0 #(mg/(ml/min**2))
        r2 =0 #(mg/(ml/min**2))
        r3 =0 #(mg/(ml/min**2))
    #Dose
    Ro = 12 #Chlorine dosing period, ever7 12 minutes
    Ro0 = 2 #Minutes duration of dose
    #Time
    Pre_runningT = 0 #Runing time variable
    K0 = 32.3 # free chrolirine demand per min 
    C= 0 # initial #(mg/L) #Concentration of Free Chrloride available
    O = 301 # Initial Oxygen demand, as per luos initial oxygen demand
    #Other parameters
    SigC = 1.70*(10**-3) #Natural decay of FC
    BC =5.38*(10**-4) #Depletion rate of FC in water. 
    A_Per =0
    List_Time_Ints = list(range(Ro,500,Ro))
    List_C=[]
    for i in Times: 
        Running_Time = i
        if(Running_Time in List_Time_Ints):
            A_Per=A_Per+1
        Time_Interval = Running_Time-(Pre_runningT)
        if 0<= Running_Time <= (0+Ro0) :
            Rate = r1
            X = 1
        elif Ro <= Running_Time <= (Ro+Ro0) :
            Rate = r2
            X = 1
        elif 2*Ro <= Running_Time <= (2*Ro+Ro0) : 
            Rate = r3
            X = 1
        elif (A_Per*Ro) <= Running_Time <= (A_Per*Ro+Ro0) : 
            Rate = r3
            X = 1
        else: 
            X = 0
        dO = K0*Time_Interval #Demand per time interval
        O = O+dO # Current oxygen demand
        decay = ((-SigC*Time_Interval)*C) - ((BC*Time_Interval)*O*C)  #Decay due to demand of chlorine
        Increase = (Rate*X*Time_Interval) #increase due to dosing period. 
        dC = decay + Increase #Total chanfe in Free Chlorine
        C = C+dC #FRee Chlorine after set time.
        if C < 0:
            C = 0 
        Pre_runningT = i #Running Time.
        if(i==10):
            print(O)
        List_C.append(C)
    Cdf = pd.DataFrame(
    {'Time': Times,
     'C': List_C,
    })
    return Cdf



'''
import seaborn as sns
df_c= F_Chloride_lvl(300)
sns.lineplot(df_c["Time"], df_c["C"])


Cdf = pd.DataFrame(
{'Time': list(range(0,301)),
 'C': 2,
})

sns.lineplot(Cdf["Time"], Cdf["C"])
'''

#CFU_Non = (TR(decimal)*(CFU non inoculatred + CFU wash Water))
def Washing_Batch(df, New_water_every_xpacks):
    Contamination_Vector = df['CFU']
    Rangeofiterations = list(range(0,len(Contamination_Vector)))
    if New_water_every_xpacks == 0:
         every_so = []
    else:
        every_so = Rangeofiterations[::New_water_every_xpacks]
    Log_Red_WashW = rng.uniform(1.87,2.23)
    TrRatetoNI = (1*10**rng.normal(0.0,0.3))/100 #check this fit
    Cont_Water =0
    for i in range(len(Contamination_Vector)):
        if i in every_so:
            Cont_Water = 0
        Cont =  Contamination_Vector[i]
        if Cont>0:
            New_Cont = Cont*10**-Log_Red_WashW
            Cont_Water = Cont - New_Cont
            Contamination_Vector[i] = New_Cont
        elif Cont ==0:
            Transfer_W_NI = Cont_Water*TrRatetoNI
            New_Cont = Transfer_W_NI
            Cont_Water = Cont_Water -Transfer_W_NI 
            Contamination_Vector[i] = New_Cont
    df["CFU"] = Contamination_Vector
    return df

#Washing

def F_Partitioning_ProcLines(gb3 , NPartitions):
    List_GB3 = []
    for j in gb3:
        if len(j) != 1:
            j = F_Partitioning_W(DF= j,NPartitions= NPartitions) 
            List_GB3.append(j)
        elif len(j) == 1:
            List_GB3.append(j)
    return List_GB3

def F_DF_Clvl(Time):
    List_Wash = []
    Times = list(range(0,Time+1))
    for i in Times:
     C = F_Chloride_lvl (i)
     List_Wash.append(C)
     
    dataTime = {
        "Time": Times,
        "Clvl":List_Wash }
    df_Clvl = pd.DataFrame(dataTime)
    return (df_Clvl)

#F_Chloride_lvl(200)


def F_Washing_ProcLines (List_GB3, Wash_Rate, Cdf):
    for j in List_GB3:
        WashT = len(j.index)
        #DF_Clvl = F_DF_Clvl(WashT)
        
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        
        Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
        alpha = 0.75#Inactivation rate of pathogen via FC L/mgmin
        V = (3200 *1000) #L #From Luo et al 2012. 
        Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
        Wash_Time = 2.3 #min 
        c1 = 1/Wash_Time #Reciprocal of average time. 
        L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
        Xl = 0
        Xw =0  #pathogen in process water MPN/ml
        
        L_Xw = []
        L_Xl = []
        for i in Times_W:
            #Defining Initial Contamination
            Time = i
            AvCont = j.at[i,"CFU"] /(j.at[i,"Weight"]*454)
            AvContAfter = AvCont*10**-0.8
            C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
            #C= F_Chloride_lvl(Time_Wash= Time)
            Bws = ((AvCont- AvContAfter)*Rate)/V
            CXw = Bws - (Blw*Xw*(L/V)) - (alpha*Xw*C)
            Xw = Xw+CXw
            if Xw<0:
                Xw = 0
            L_Xw.append(Xw)
            Xl = AvCont
            CXl = (Blw*Xw) - (alpha*Xl*C) - (c1*Xl)
            Xl =Xl +CXl
            if Xl < 0:
                Xl = 0
            L_Xl.append(Xl)
            AvCont = Xl
            CFU_2 = AvCont*((j.at[i,"Weight"]*454))
            j.at[i,"CFU"] =  CFU_2
    return (List_GB3) 


def F_Washing_ProcLines2 (List_GB3, Wash_Rate, Cdf):
    for j in List_GB3:
        WashT = len(j.index)
        #DF_Clvl = F_DF_Clvl(WashT)
        
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        
        Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
        alpha = 0.75#Inactivation rate of pathogen via FC L/mgmin
        V = (3200 *1000) #L #From Luo et al 2012. 
        Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
        Wash_Time = 2.3 #min 
        c1 = 1/Wash_Time #Reciprocal of average time. 
        L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
        Xl = 0
        Xw =0  #pathogen in process water MPN/ml
        
        L_Xw = []
        L_Xl = []
        for i in Times_W:
            Xs = np.random.triangular(0.003,0.055,0.149)
            index =i
            #Defining Initial Contamination
            Time = i
            AvCont = j.at[i,"CFU"] /(j.at[i,"Weight"]*454)
            #AvCont_CFU = df.at[i,"CFU"]
            #AvContAfter = AvCont*10**-0.8
            C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
            Bws = (((AvCont)-(AvCont*Xs))*Rate*1000)/V
            #Bws = ((AvCont- AvContAfter)*Rate)/V
            #print(Bws)
            CXWfirst = Bws - (Blw*Xw*(L/V))
            CXw =  CXWfirst - (alpha*Xw*C)
            Xw = Xw+CXw
            if Xw<0:
                Xw = 0
            L_Xw.append(Xw)
            Xl = (AvCont*Xs)
            #print(Xl)
            CXL23t = (alpha*Xl*C) - (c1*Xl)
            #print(CXL23t, "CXL23t")
            if CXL23t>Xl:
                Xl = 0
            CXl = (Blw*Xw)  #- (alpha*Xl*C) - (c1*Xl)
            #print(Blw*Xw, "fist section")
            #print(CXl, "CXL")
            #print(Xl, "XL")
            Xl =Xl +CXl
            if Xl < 0:
                Xl = 0
            L_Xl.append(Xl)
            AvCont = Xl
            CFU_2 = AvCont*((j.at[index,"Weight"]*454))
            if CFU_2<1:
                CFU_2= np.random.binomial(1,CFU_2)
            j.at[index,"CFU"] =  CFU_2
    return (List_GB3) 

def F_Washing_ProcLines3 (List_GB3, Wash_Rate, Cdf):
    for j in List_GB3:
        WashT = len(j.index)
        #DF_Clvl = F_DF_Clvl(WashT)
        
        Times_W = np.arange(0, WashT, 1).tolist()
        Times_W = [round(num, 1) for num in Times_W]
        
        Blw = 0.38 #ml/g min: is the pathogen binding rate to pieces of shredded lettuce heads
        alpha = 0.75#Inactivation rate of pathogen via FC L/mgmin
        V = (3200 *1000) #L #From Luo et al 2012. 
        Rate = Wash_Rate/2.2  #45.45 #kg/min #From Luo et al 2012. 
        Wash_Time = 2.3 #min 
        c1 = 1/Wash_Time #Reciprocal of average time. 
        L = (Rate*1000)/(c1) #g of lettuce in the tak at the same time
        Xl = 0
        Xw =0  #pathogen in process water MPN/ml
        
        L_Xw = []
        L_Xl = []
        for i in Times_W:
            Xs = np.random.triangular(0.003,0.055,0.149)
            index =i
            #Defining Initial Contamination
            Time = i
            AvCont = j.at[i,"CFU"] /(j.at[i,"Weight"]*454)
            #AvCont_CFU = df.at[i,"CFU"]
            #AvContAfter = AvCont*10**-0.8
            C =   float(Cdf.loc[Cdf['Time'] == Time, 'C'])
            Bws = (((AvCont)-(AvCont*Xs))*(Rate*1000))/V
            
            #Bws = ((AvCont- AvContAfter)*Rate)/V
            #print(Bws)
            CXWfirst = Bws - (Blw*Xw*(L/V))
            CXw =  CXWfirst - (alpha*Xw*C)
            Xw = Xw+CXw
            if Xw<0:
                Xw = 0
            L_Xw.append(Xw)
            Xl = (AvCont*Xs)
            #print(Xl)
            #CXL23t = (alpha*Xl*C) - (c1*Xl)
            #print(CXL23t, "CXL23t")
            if C>0.5:
                CXL23t = 0.214*np.log(C)+0.220
            if C<0.5:
                CXL23t = 0
            Xl = Xl*(10**-CXL23t)
            CXl = (Blw*Xw)  #- (alpha*Xl*C) - (c1*Xl)
            #print(Blw*Xw, "fist section")
            #print(CXl, "CXL")
            #print(Xl, "XL")
            Xl =Xl +CXl
            if Xl < 0:
                Xl = 0
            L_Xl.append(Xl)
            AvCont = Xl
            CFU_2 = AvCont*((j.at[index,"Weight"]*454))
            if CFU_2<1:
                CFU_2= np.random.binomial(1,CFU_2)
            elif CFU_2>1:
                partial =math.modf(CFU_2)
                part1= np.random.binomial(1,partial[0])
                part2= partial[1]
                CFU_2 = part1+part2
            j.at[index,"CFU"] =  CFU_2
    return (List_GB3) 


#Adding production hours to final product

def Production_Hours(df,FP_Hour):
    df.reset_index(drop =True, inplace= True)
    Part_per_hour = FP_Hour/df["Weight"][0]
    Total_Part = len(df.index)
    
    Hours_Prod = Total_Part/Part_per_hour
    
    #If the ones do not add to a whole number. 
    Splits =math.modf(Hours_Prod)
    Complete_H= int(Splits[1])
    Partial = Splits[0]
    
    Partial_Parts = Part_per_hour*Partial
    
    
    #Creating the sequence
    Seq_P1=np.repeat(range(1,Complete_H+1),Part_per_hour)
    Seq_P2 = np.repeat(Complete_H+1,Partial_Parts)
    
    Full_Seq =np.concatenate((Seq_P1,Seq_P2))
    df["HourProd"] = Full_Seq
    return df


#Final Product Case

def Case_Packaging(df,Case_Weight, Pack_Weight):

    Packages_Case = Case_Weight/Pack_Weight
    Total_Packages = len(df.index)
    Total_Cases = Total_Packages/Packages_Case
    Case_Pattern = [i for i in range(1, int(Total_Cases)+1) for _ in range(int(Packages_Case))]
    Crop_No = len(df.index)
    Case_Pattern=Case_Pattern[:Crop_No]
    df.insert(1,"Case",Case_Pattern)
    
    Total_Pallets = Total_Packages/320
    Pallet_Pattern = [i for i in range(1, int(Total_Pallets)+2) for _ in range(int(320))]
    Pallet_Pattern=Pallet_Pattern[:Crop_No]
    df.insert(1,"PalletNoFS",Pallet_Pattern)
    
    return df



#OTher Functions

def F_SummingGB2Cont(gb2):
    List_x_Sum=[]  
    for i in gb2:
        x_Sum=i.CFU.sum()
        List_x_Sum.append(x_Sum)
    Out =sum(List_x_Sum)
    return Out


