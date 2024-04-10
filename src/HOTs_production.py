'''
HOTs_production.py

HOTs cruise 347 



created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser/.....
'''

from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   
import seaborn as sns


###############################

#  Importing the data frame

###############################


df_all = pd.read_excel("../data/HOTs_main_assays-compiled.xlsx",sheet_name = 'Main-Assays', header = 0)


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time (hrs)':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_all.fillna(0)

#df_all['D_tech1'].mask(df_all['D_tech1']<1,inplace=True)
#creating log and stats for data 
df_all['Alog1'] = np.log(df_all['A_tech1'])
df_all['Alog2'] = np.log(df_all['A_tech2'])
df_all['Blog1'] = np.log(df_all['B_tech1'])
df_all['Blog2'] = np.log(df_all['B_tech2'])
df_all['Clog1'] = np.log(df_all['C_tech1'])
df_all['Clog2'] = np.log(df_all['C_tech2'])
df_all['Dlog1'] = np.log(df_all['D_tech1'])
df_all['Dlog2'] = np.log(df_all['D_tech2'])

df_all.fillna(100)

df_all['A_mean'] =  np.nanmean(np.r_[[df_all[i] for i in ['A_tech1', 'A_tech2']]],axis=0)
df_all['B_mean'] =  np.nanmean(np.r_[[df_all[i] for i in ['B_tech1', 'B_tech2']]],axis=0)
df_all['C_mean'] =  np.nanmean(np.r_[[df_all[i] for i in ['C_tech1', 'C_tech2']]],axis=0)
df_all['D_mean'] =  np.nanmean(np.r_[[df_all[i] for i in ['D_tech1', 'D_tech2']]],axis=0)


df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['A_tech1', 'A_tech2','B_tech1', 'B_tech2','C_tech1', 'C_tech2','D_tech1', 'D_tech2']]],axis=0)
df_all['sigma'] = np.std(np.r_[[df_all[i] for i in ['A_tech1', 'A_tech2','B_tech1', 'B_tech2','C_tech1', 'C_tech2','D_tech1', 'D_tech2']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['Alog1','Alog2','Blog1','Blog2','Clog1','Clog2','Dlog1','Dlog2']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['Alog1','Alog2','Blog1','Blog2','Clog1','Clog2','Dlog1','Dlog2']]],axis=0)


##############################

#   data arrays

###############################


sns.relplot(data=df_all, x="abundance", y="sigma", palette = 'cool',size ='ID', hue="organism", style = 'ID', markers =True,  kind="scatter").set(title='Raw Data' )

sns.relplot(data=df_all, x="log_abundance", y="log_sigma", palette = 'cool',size ='ID', hue="organism", style = 'ID', markers =True,  kind="scatter").set(title='Log Data' )



assays = df_all['ID'].unique()
nassays = assays.shape[0]
treats  = df_all['treatment'].unique()
ntreats = treats.shape[0]


colors = ('r','g')
markers = ('o','*')

##############################

#    graphing the data 

##############################
'''

for a,n in zip(assays,range(nassays)):
    fig1,(ax1)= plt.subplots(1,2, figsize = (12,7))
    fig2,(ax2) = plt.subplots(1,2,figsize = (12,7))
    for t,nt in zip(treats,range(ntreats)): 
        df = df_all[(df_all['ID'] == a)]
        count = nt
        df = df[(df['ID']==a) & (df['Treatment']==t)]
        if df.shape[0] == 0:
            print('empty df in'+str(a))
            continue
        if df.shape[0] > 0:
            pass
        pdf = df[df['organism'] == 'P']
        hdf = df[df['organism'] == 'H']
        ax1[0].plot(pdf['time'], pdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        ax1[1].plot(hdf['time'], hdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        fig1.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)
        fig1.suptitle('Raw Dynamics '+ str(a))
        ax1[1].set_ylabel('HOOH concentration (\u03BCM)')
        ax1[1].set_xlabel('Time (days)')
        ax1[0].set_ylabel('Pro abundance (cell/ml)')
        ax1[0].set_xlabel('Time (days)')
        #fig1.supxlabel('Time (days)')
        l1 = ax1[0].legend(loc = 'lower left', prop={"size":12}) 
        l1.draw_frame(False)#print(df)
        ax1[0].semilogy()
        ax1[1].semilogy()
        #log data
        ax2[0].plot(pdf['time'], pdf['log_abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        ax2[1].plot(hdf['time'], hdf['log_abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        fig2.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)
        fig2.suptitle('Log Dynamics '+ str(a))
        ax2[1].set_ylabel('HOOH concentration (\u03BCM)')
        ax2[1].set_xlabel('Time (days)')
        ax2[0].set_ylabel('Pro abundance (cell/ml)')
        ax2[0].set_xlabel('Time (days)')
        #fig2.supxlabel('Time (days)')
        l2 = ax1[0].legend(loc = 'center right', prop={"size":12}) 
        l2.draw_frame(False)#print(df)
        ax2[0].semilogy()
        ax2[1].semilogy()
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()
    fig1.savefig('../figures/HproductionRaw_'+str(a)+'_data.png')
    fig2.savefig('../figures/HproductionLog_'+str(a)+'_data.png')



df_H =df_all[df_all['Treatment']== 'HEPES']
df_T =df_all[df_all['Treatment']== 'TAPS']

df = df_H

fig, (ax1,ax2) = plt.subplots(1,2,figsize = (10,6))
fig.subplots_adjust(right=0.95, wspace = 0.45, left = 0.1, hspace = 0.30, bottom = 0.2)
fig.suptitle('Batch cultures in Hepes')
ax2.set_ylabel('HOOH concentration (\u03BCM)')
ax2.set_xlabel('Time')
ax1.set_ylabel("Pro abundance (cell/ml)")
ax1.set_xlabel('Time (days)')
ax2.semilogy()
ax1.semilogy()
ax2.plot( df[df['ID'] == 'UH18301+EZ55'][df['organism'] == 'H']['time'],  df[df['ID'] == 'UH18301+EZ55'][df['organism'] == 'H']['log_abundance'],marker = 's', color  = 'b', label = 'H in UH18301+EZ55')
ax2.plot( df[df['ID'] == 'UH18301'][df['organism'] == 'H']['time'],  df[df['ID'] == 'UH18301'][df['organism'] == 'H']['log_abundance'], marker = 's',color  = 'r', label = 'H in UH18301')
ax1.plot( df[df['ID'] == 'UH18301+EZ55'][df['organism'] == 'P']['time'],  df[df['ID'] == 'UH18301+EZ55'][df['organism'] == 'P']['log_abundance'], marker = 'd',color  = 'b', label = 'P in UH18301+EZ55')
ax1.plot( df[df['ID'] == 'UH18301'][df['organism'] == 'P']['time'],  df[df['ID'] == 'UH18301'][df['organism'] == 'P']['log_abundance'], marker = 'd',color  = 'r', label = 'P in UH18301')
plt.legend()
plt.show

'''


print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
