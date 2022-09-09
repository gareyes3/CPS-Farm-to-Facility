# -*- coding: utf-8 -*-
"""
Created on Tue May 10 10:06:02 2022

@author: Gustavo Reyes
"""
import pandas as pd
import numpy as np

Col_Names_NI = "BaselineNI PHS4D_NI PHS4H_NI PHSInt_NI HTrad_NI RSTrad_NI FPSTrad_NI CS_NI".split()
Col_Names_AI = "BaselineAI PHS4D_AI PHS4H_AI PHSInt_AI HTrad_AI RSTrad_AI FPSTrad_AI CS_AI".split()


AI_B = Outputs_Df_AI_1["Final_CFU_Acc_Portion_mean"]


AI_10 = Outputs_Df_AI_3["Final_CFU_Acc_Portion_mean"]
AI_1 = Outputs_Df_AI_2["Final_CFU_Acc_Portion_mean"]

df_AI_B = pd.DataFrame({"Expossure": AI_B,
               "Cont Scenario": "Random Contamination",
               "Sampling Plan": Col_Names_AI,
               "Baseline System": "All-Intervention"
    })


df_AI_10 = pd.DataFrame({"Expossure": AI_10,
               "Cont Scenario": "10% Cluster Contamination",
               "Sampling Plan": Col_Names_AI,
               "Baseline System": "All-Intervention"
    })

df_AI_1 = pd.DataFrame({"Expossure": AI_1,
               "Cont Scenario": "1% Cluster Contamination",
               "Sampling Plan": Col_Names_AI,
               "Baseline System": "All-Intervention"
    })


NI_B = Outputs_Df_NI_1["Final_CFU_Acc_Portion_mean"]
NI_10 = Outputs_Df_NI_3["Final_CFU_Acc_Portion_mean"]
NI_1 = Outputs_Df_NI_2["Final_CFU_Acc_Portion_mean"]

df_NI_B = pd.DataFrame({"Expossure": NI_B,
               "Cont Scenario": "Random Contamination",
               "Sampling Plan": Col_Names_NI,
               "Baseline System": "No-Intervention"
    })


df_NI_10 = pd.DataFrame({"Expossure": NI_10,
               "Cont Scenario": "10% Cluster Contamination",
               "Sampling Plan": Col_Names_NI,
               "Baseline System": "No-Intervention"
    })

df_NI_1 = pd.DataFrame({"Expossure": NI_1,
               "Cont Scenario": "1% Cluster Contamination",
               "Sampling Plan": Col_Names_NI,
               "Baseline System": "No-Intervention"
    })




final_exp = pd.concat([df_AI_B,df_AI_10,df_AI_1,df_NI_B,df_NI_10,df_NI_1])

final_exp.to_csv("C:\\Users\\Gustavo Reyes\\Box\\CPS Project- Farm to Facility\\Papers\\CSV Data\\final_exp.csv")
