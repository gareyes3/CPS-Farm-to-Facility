# -*- coding: utf-8 -*-
import pandas as pd 
import numpy as np
import itertools
import scipy.stats as stats
import math
import os
import random
import seaborn as sns
from matplotlib import pyplot as plt
import Funz
import Inputz
from numpy.random import Generator, PCG64
rng = Generator(PCG64())

#%% Contamination Scenarios

def F_Background_C(df, Hazard_lvl, Partition_Units):
    Ci = Hazard_lvl/Partition_Units
    df.CFU =Ci
    return df

def F_Point_Source_C (df, Hazard_lvl, No_Cont_Clusters, Cluster_Size, Partition_Weight):
    No_Cont_PartitionUnits = int(Cluster_Size/Partition_Weight) #Cluster Units per Contaminatrion Cluster
    Hazard_lvl_C= Hazard_lvl/No_Cont_Clusters
    Ci = Hazard_lvl_C
    X_1= df.sample(No_Cont_Clusters)
    X_1 = list(X_1.PartitionID)
    df.loc[X_1,'CFU']= df['CFU'] + Ci
    return df
    
def F_systematic_C(df, Hazard_lvl,No_Cont_Clusters,Cluster_Size, Partition_Weight):
    #number of contminated partitions per cluster
    No_Cont_PartitionUnits = int(Cluster_Size/Partition_Weight) 
    #Determining the hazard level per cluster
    Hazard_lvl_percluster= Hazard_lvl / No_Cont_Clusters #(No_Cont_PartitionUnits*No_Cont_Clusters)
    for i in range(0,No_Cont_Clusters):
        n = random.randint(0,len(df.index)- No_Cont_PartitionUnits)
        x_random_consecutive_rows = df[n:n + No_Cont_PartitionUnits]
        x_random_consecutive_rows = list(x_random_consecutive_rows.index)
        #using multinomial creating array to contaminated the desired outputs
        #probaility defined by the number of partitions. 
        Contamination_Pattern = rng.multinomial(Hazard_lvl_percluster,[1/No_Cont_PartitionUnits]*No_Cont_PartitionUnits,1)
        df.loc[x_random_consecutive_rows,'CFU']= df.loc[x_random_consecutive_rows,'CFU'] + Contamination_Pattern[0]
    return df
        
def F_Crew_C(df, Hazard_lvl,No_Cont_Clusters, Cluster_Size, Partition_Weight):
        No_Cont_PartitionUnits = int(Cluster_Size/Partition_Weight)
        Hazard_lvl_C= Hazard_lvl/(No_Cont_PartitionUnits*No_Cont_Clusters)
        Ci = Hazard_lvl_C
        x_2 =[]
        PartitionIDs = list(df['PartitionID'])
        for i in range(No_Cont_Clusters):
            TrimPartitionIDS = PartitionIDs[:len(PartitionIDs)-No_Cont_PartitionUnits]
            n = random.sample(TrimPartitionIDS,1)
            n  = n[0]
            x_random_consecutive_rows = PartitionIDs[n:n + No_Cont_PartitionUnits]
            x_2.append(x_random_consecutive_rows)
            PartitionIDs = [x for x in PartitionIDs if x not in x_random_consecutive_rows]
        x_2 =[j for i in x_2 for j in i]
        df.loc[ x_2,'CFU']= df['CFU'] + Ci
        return df


def F_Harvester_C(df, Hazard_lvl,No_Cont_Clusters, Cluster_Size, Partition_Weight):
        No_Cont_PartitionUnits = int(Cluster_Size/Partition_Weight)
        Hazard_lvl_C= Hazard_lvl/(No_Cont_PartitionUnits*No_Cont_Clusters)
        Ci = Hazard_lvl_C
        x_2 =[]
        PartitionIDs = list(df['PartitionID'])
        for i in range(No_Cont_Clusters):
            TrimPartitionIDS = PartitionIDs[:len(PartitionIDs)-No_Cont_PartitionUnits]
            n = random.sample(TrimPartitionIDS,1)
            n  = n[0]
            x_random_consecutive_rows = PartitionIDs[n:n + No_Cont_PartitionUnits]
            x_2.append(x_random_consecutive_rows)
            PartitionIDs = [x for x in PartitionIDs if x not in x_random_consecutive_rows]
        x_2 =[j for i in x_2 for j in i]
        df.loc[ x_2,'CFU']= df['CFU'] + Ci
        return df
    
    
def F_PEC_C(gb2,Hazard_lvl,Processing_Lines, Lines_Cont):     
        L_Lines_Cont= random.sample(range(Processing_Lines),Lines_Cont)
        Hazard_PLine = Hazard_lvl/Lines_Cont
        Hazard_Pall=0
        for i in L_Lines_Cont:
            DfPick = gb2[i]
            Hazard_Pall=Hazard_PLine/ len(DfPick.index)
            DfPick['CFU'] =DfPick['CFU']+Hazard_Pall
        return gb2
    
