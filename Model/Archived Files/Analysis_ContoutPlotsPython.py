# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:18:08 2021

@author: gareyes3
"""

#%%
import sys, os
sys.path
sys.path.append('C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model')
#sys.path.append('C:\\Users\gareyes3\Documents\GitHub\CPS-Farm-to-Facility\Model')

os.listdir()
#%%
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
Data1=pd.read_csv(r'C:\\Users\Gustavo Reyes\Documents\GitHubFiles\CPS-Farm-to-Facility\Model\3d-df-11-13.csv')

#%%
Summarized= Data1.groupby(['Clusters','ContLevel','Grabs','SampleMass']).mean()
Summarized["CFUperR"] =1-Summarized["PH_CFU_PerA"]


#%%
P1 = Summarized[0:25]
Matrix=P1.pivot_table( index = "Grabs",columns="SampleMass",values = 'CFUperR')

Matrix=Matrix.to_numpy()

#%%
fig = go.Figure(data =
    go.Contour(
        z=Matrix,
        x=[60,120,300,600,1200], # horizontal axis
        y=[1,30,60,120,300] # vertical axis
    ))

fig.show()



