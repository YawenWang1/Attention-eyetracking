from __future__ import division
import pandas as pd
import numpy as np
import sndlib as sl

screenDistance = 57 # cm
DegDis         = 4
PixPerCm       = 27.7
PosC           = [4,7]
position       = np.round(np.tan(np.deg2rad(PosC)) * screenDistance * PixPerCm)



Esize, Cpos, Espace, Cchg, CueFile, Sori = [1,1.5], [4,7], [1.5,3],[0,+1, +3, +7, -1, -3,-7],['<','>'],[-15,15]
NumOfElements = [1,6]

SF = 3
Tcond       = []
for n in range(len(Esize)):
    if Esize[n] == np.max(Esize):
        SF = 3
    else:
        SF = 3 * (Esize[1]/Esize[0])
    for m in range(len(Espace)):
        RadiusOfCircle = [0,Espace[m]]
        p = [[RadiusOfCircle[g]*np.cos(h*(2*np.pi/NumOfElements[g])),RadiusOfCircle[g]*np.sin(h*(2*np.pi/NumOfElements[g]))] 
            for g in range(len(RadiusOfCircle)) for h in range(NumOfElements[g])]
        for i in range(len(CueFile)):
            for j in range(len(Cchg)):
                if Cchg[j] > 0:
                    corrans = 'right'
                elif Cchg[j] < 0:
                    corrans = 'left'
                else:
                    corrans = 'Vertical'
                for k in range(len(Sori)):
                    t = [Cchg[j]]
                    t.extend([Sori[k] for s in range(6)])
                    Tcond.append({'esize':Esize[n],'pos':p,'Cue':CueFile[i], 'cchg':Cchg[j],'allOris':t,'corrans':corrans,'sf':SF})
                    
np.random.shuffle(Tcond)
#asize,aoris,apos,acue,acchg = np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
asize,aoris,apos,acue,acchg,acorrans, asf = [],[],[],[],[],[],[]

for i in range(len(Tcond)):
    asize.extend([Tcond[i]['esize']])
    aoris.extend([Tcond[i]['allOris']])
    apos.extend([Tcond[i]['pos']])
    acchg.extend([Tcond[i]['cchg']])
    acue.extend([Tcond[i]['Cue']])
    acorrans.extend([Tcond[i]['corrans']])
    asf.extend([Tcond[i]['sf']])


newpath = 'D:\Yawen\VBehavior'
os.chdir(newpath)


#asize,aoris,apos,acue,acchg = np.array(asize),np.array(aoris),np.array(apos), np.array(acue),np.array(acchg)
    
d = {'asize':asize,'aoris':aoris,'apos':apos,'acchg':acchg,'acue':acue,'acorrans':acorrans,'asf':asf}

#create a pandas dataframe from the data
df = pd.DataFrame(data=d)

#create a pandas excel writer using xlsxwriter as the engine
writer = pd.ExcelWriter('VTiltIllusionContios.xlsx',engine='xlsxwriter')
#convert the dataframe to an xlsxwriter excel object
df.to_excel(writer,sheet_name='Sheet1',index = False)

writer.save()