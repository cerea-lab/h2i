#-----------------------------------------------------------------------
#     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
#     The H2I is distributed under the GNU General Public License v3
#-----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import statistics
import math

shift = 60
delta = 100

name_list = ['O3','NO','NO2','HONO','HNO3ad','HONOad','NO2ad','NOad'] # list of species to visualize

# READS SIMULATION DATA
# ----- concentrations simulated in the shaded box
shade_file = 'dem_Exp_shade_demo-test.dat' 
sconc = {'O3':[],'NO':[],'NO2':[],'HONO':[],'NOx':[],'HNO3ad':[],'HONOad':[],'NO2ad':[],'NOad':[]}
stime = []
file = open(shade_file,"r")

get_time = 0
for line in file:
    lect = line.split()
    if lect :
        if lect[0] == '======':
            time = float(lect[-1])
            if (time - shift)%delta == 0:
                get_time = 1
                stime.append(time)
                if (len(sconc['NO']) > 0): # if this is not the first file read
                    sconc['NOx'].append(NOx)
                NOx = 0
            else:
                get_time = 0
        else:
            name = lect[0]
            if get_time:
                if name in name_list:
                    val = lect[1]
                    sconc[name].append(float(val))
                if (name == 'NO'):
                    NOx = NOx + float(val)
                if (name == 'NO2'):
                    NOx = NOx + float(val)                
sconc['NOx'].append(NOx)
NOx = 0
NOy = 0

# ----- concentrations simulated in the sunlit box
light_file = 'dem_Exp_light_demo-test.dat' 
lconc = {'O3':[],'NO':[],'NO2':[],'HONO':[],'NOx':[],'HNO3ad':[],'HONOad':[],'NO2ad':[],'NOad':[]}
file = open(light_file,"r")
get_time = 0
for line in file:
    lect = line.split()
    if lect :
        if lect[0] == '======':
            time = float(lect[-1])
            if (time - shift)%delta == 0:
                get_time = 1
                if (len(lconc['NO']) > 0): 
                    lconc['NOx'].append(NOx)
                NOx = 0
            else:
                get_time = 0
        else:
            name = lect[0]
            if get_time:
                if name in name_list:
                    val = lect[1]
                    lconc[name].append(float(val))
                if (name == 'NO'):
                    NOx = NOx + float(val)
                if (name == 'NO2'):
                    NOx = NOx + float(val)
lconc['NOx'].append(NOx)


htime = []
for i in range(len(stime)):
    tmp = (stime[i]%86400)/3600 # conversion from seconds to hours
    htime.append(tmp)
    

# READS EXPERIMENTAL DATA
file = open('INORG.txt',"r") # measured concentrations
INORG_list = ['NO','NO2','HONO','O3','NOx'] # species measured in file INORG.txt
mes_INORG = {'O3':[],'NO':[],'NO2':[],'HONO':[],'NOx':[]}
time_INORG = []
for line in file :
    lect = line.split()
    time_INORG.append(float(lect[0])) # 1st column = hours
    for i in range(len(lect)-1):
        mes_INORG[INORG_list[i]].append(float(lect[i+1]))
    mes_INORG['NOx'].append(float(lect[1])+float(lect[2]))

    
# VISUALIZATION
fig = plt.figure()
for i in range(len(INORG_list)):
    plt.subplot(3,3,i+1)
    plt.plot(htime,sconc[INORG_list[i]],'-', ms=2)
    plt.plot(htime,lconc[INORG_list[i]],'--', ms=2)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
    plt.plot(time_INORG[:len(mes_INORG[INORG_list[i]])],mes_INORG[INORG_list[i]],'.', ms=2)        
    plt.title(INORG_list[i])
    if (i%3==0) :
        plt.ylabel('C (µg.m$^{-3}$)')
# NO2ad
plt.subplot(3,3,len(INORG_list)+1)
plt.plot(htime,sconc['NO2ad'],'.', ms=2)
plt.plot(htime,lconc['NO2ad'],'.', ms=2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
plt.title('NO2ad')
# NOad
plt.subplot(3,3,len(INORG_list)+2)
plt.plot(htime,sconc['NOad'],'.', ms=2)
plt.plot(htime,lconc['NOad'],'.', ms=2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
plt.title('NOad')
plt.xlabel('Time (GMT hour)')
plt.ylabel('C (µg.m$^{-3}$)')
# HNO3ad
plt.subplot(3,3,len(INORG_list)+3)
plt.plot(htime,sconc['HNO3ad'],'.', ms=2)
plt.plot(htime,lconc['HNO3ad'],'.', ms=2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
plt.title('HNO3ad')
plt.xlabel('Time (GMT hour)')
# HONOad
plt.subplot(3,3,len(INORG_list)+4)
plt.plot(htime,sconc['HONOad'],'.', ms=2)
plt.plot(htime,lconc['HONOad'],'.', ms=2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.title('HONOad')
plt.xlabel('Time (GMT hour)')
fig.tight_layout()


plt.show()


# STATISTICS
# files containing simulated and measured concentrations 
# have identical start times and time steps;
for i in range(len(INORG_list)):
    Ndata = min(len(lconc[INORG_list[i]]), len(mes_INORG[INORG_list[i]])) 
    mL = statistics.mean(lconc[INORG_list[i]])
    mS = statistics.mean(sconc[INORG_list[i]])
    mSim = int((mL + mS)/2*100)/100
    mExp = int(statistics.mean(mes_INORG[INORG_list[i]])*100)/100
    Err = abs(int((mSim - mExp)/mExp * 10000)/100)
    count = 0
    MNBE = 0
    MNGE = 0
    RMSE = 0
    for j in range(Ndata):
        if mes_INORG[INORG_list[i]][j]==0 :
            count = count + 1
            continue
        model = (lconc[INORG_list[i]][j]+sconc[INORG_list[i]][j])/2
        MNBE = MNBE + (model-mes_INORG[INORG_list[i]][j])/mes_INORG[INORG_list[i]][j]
        MNGE = MNGE + abs(model-mes_INORG[INORG_list[i]][j])/mes_INORG[INORG_list[i]][j]
        RMSE = RMSE + (mes_INORG[INORG_list[i]][j] - model)**2
    MNBE = int( MNBE/(Ndata-count)*100 * 100)/100
    MNGE = int( MNGE/(Ndata-count)*100 * 100)/100
    RMSE = int( math.sqrt(RMSE/(Ndata-count)) * 100)/100
    print(INORG_list[i], ': ')
    print('mExp, mSim, Err  :  ', mExp, '  ', mSim, '  ', Err, '%')
    print('MNBE, MNGE, RMSE :', MNBE, '% ', MNGE, '% ', RMSE)
    print('')

plt.show()

