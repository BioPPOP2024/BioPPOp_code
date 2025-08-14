#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

#from descriptors import compute_all_features
from stabilogram.stato import Stabilogram

# In[2]:


forceplate_file_selected = "test.csv"


# In[3]:


data_forceplatform = pd.read_csv(forceplate_file_selected,header=[31],sep=",",index_col=0)
#print(data_forceplatform.head(100))


# In[4]:

pd.Series(dtype='float64')
dft = data_forceplatform

X = dft.get(" Cx")
Y = dft.get(' Cy')

Cx1=X*1.5
Cy1=Y*1.5
dft.loc[:, "Cx1"]=Cx1
dft.loc[:, "Cy1"]=Cy1
X1 = dft.Cx1
Y1 = dft.Cy1

#print(dft.head(100))
X = X - np.mean(X)
Y = Y - np.mean(Y)
X1 = X1 - np.mean(X1)
Y1 = Y1 - np.mean(Y1)

X = 100*X
Y = 100*Y
X1 = 100*X1
Y1 = 100*Y1


RightSwayLs=list()
for i in range(len(X)-1):
    Xdis=X[i+1]-X[i]
    Ydis=Y[i+1]-Y[i]
    SwayL=np.sqrt(Xdis**2+Ydis**2)
    RightSwayLs.append(SwayL)
#print(SwayLs)
RightSwayLs=pd.Series(RightSwayLs)

LeftSwayLs=list()
for i in range(len(X)-1):
    Xdis=X1[i+1]-X1[i]
    Ydis=Y1[i+1]-Y1[i]
    SwayL=np.sqrt(Xdis**2+Ydis**2)
    LeftSwayLs.append(SwayL)
#print(SwayLs)
LeftSwayLs=pd.Series(LeftSwayLs)




#SwayLs=pd.Series(SwayLs)
#print(type(SwayLs))
#print(type(X))
#print(SwayLs)
time= pd.Series(range(0,3000))
time=time/1000
RightSwayLs=RightSwayLs.to_numpy()[4000:7000]
RMeanVel=RightSwayLs/3
plt.plot(time,RMeanVel)
#LeftSwayLs=LeftSwayLs.to_numpy()[4000:7000]
#LMeanVel=LeftSwayLs/3
#plt.plot(LMeanVel)
plt.xlabel('Time (s)', fontsize=15)
plt.show()

X = X.to_numpy()[4000:7000]
Y= Y.to_numpy()[4000:7000]
X1 = X1.to_numpy()[4000:7000]
Y1= Y1.to_numpy()[4000:7000]
plt.plot(X,Y,color='b')
plt.plot(X1,Y1,color='r')
plt.plot(0,0,marker='o',markersize=8, markerfacecolor='green')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')



#print(SwayLs)

# In[5]:

#fig, ax = plt.subplots(1)
#ax.plot([X])

#ax.plot([X1])
plt.legend(['Right Leg','Left Leg','Centre of Forceplate'],loc='upper left')
plt.axis([-3,3,-3,3])
plt.xlabel('Mediolateral (cm)', fontsize=15)
plt.ylabel('Anteroposterior (cm)', fontsize=15)
plt.savefig('CoP.png')

plt.show()

Timez=pd.Series(range(0,10))
Timez=Timez*10

plt.plot(Timez,Timez)
plt.xlabel('Normalised Movement Time (%)', fontsize= 15)
plt.show()

# In[6]:


data = np.array([X,Y]).T


# In[7]:

# Verif if NaN data
valid_index = (np.sum(np.isnan(data),axis=1) == 0)

if np.sum(valid_index) != len(data):
    raise ValueError("Clean NaN values first")


# In[8]:


stato = Stabilogram()
stato.from_array(array=data, original_frequency=100)


# In[9]:


fig, ax = plt.subplots(1)
ax.plot(stato.medio_lateral)
ax.plot(stato.antero_posterior)


# In[10]:

sway_density_radius = 0.3 # 3 mm

params_dic = {"sway_density_radius": sway_density_radius}

#features = compute_all_features(stato, params_dic=params_dic)


# In[11]:


#print(features)

