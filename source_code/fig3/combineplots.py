#script used to generate the plots out of the data

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#mpl.rcParams['axes.spines.right'] = False
#mpl.rcParams['axes.spines.top'] = False

mpl.rcParams['axes.linewidth'] = 4.0
nspins=[50,200]

SMALL_SIZE = 30
MEDIUM_SIZE = 33
BIGGER_SIZE = 35
LW=7
CS=10
MS=15

primer='mol'
output='Analytical'+primer
imp=dict()
analytical100=np.loadtxt('n100stabel.csv',delimiter=',')
analytical50=np.loadtxt('n50stabel.csv',delimiter=',')
analyticalunstabel100=np.loadtxt('n100unstabel.csv',delimiter=',')
analyticalunstabel50=np.loadtxt('n50unstabel.csv',delimiter=',')
analyticalunstabel2_100=np.loadtxt('n100unstabel2.csv',delimiter=',')
analyticalunstabel2_50=np.loadtxt('n50unstabel2.csv',delimiter=',')

analytical200=np.loadtxt('n200stabel.csv',delimiter=',')
analyticalunstabel200=np.loadtxt('n200unstabel.csv',delimiter=',')
analyticalunstabel2_200=np.loadtxt('n200unstabel2.csv',delimiter=',')

#load for embeded subplot
spins=np.loadtxt('spins400.csv',delimiter=',')
analytical_spins=np.loadtxt('spins_analytic.csv',delimiter=',')
zeroline_spins=np.array(([analytical_spins[0,-1],0],[2,0]))

unStabel50=np.vstack([analyticalunstabel2_50,analyticalunstabel50])
unStabel100=np.vstack([analyticalunstabel2_100,analyticalunstabel100])
unStabel200=np.vstack([analyticalunstabel2_200,analyticalunstabel200])

zeroline50=np.array(([unStabel50[5,0],0],[0.5,0]))
zeroline100=np.array(([unStabel100[5,0],0],[0.5,0]))
zeroline200=np.array(([unStabel200[5,0],0],[0.5,0]))

unStabel200=np.delete(unStabel200,range(6),0)
analytical50=np.delete(analytical50,-1,0)
unStabel200=np.delete(unStabel200,range(5,7),0)

print(unStabel50)
for n in nspins:
    output+=str(n)
    file=primer+str(n)+'.csv'
    imp[n]=np.loadtxt(file,delimiter=',')

plt.figure()
plt.xlabel('Temperature ($k_B T$)',fontsize=BIGGER_SIZE)
plt.ylabel('Magnatisation ($m$)',fontsize=BIGGER_SIZE)
plt.xticks(fontsize=MEDIUM_SIZE)
plt.yticks(fontsize=MEDIUM_SIZE)
color=['#1f77b4','#ff7f0e','#2ca02c']
i=0
for n in nspins:
    x=imp[n][:,0]
    y=imp[n][:,1]
    sd=np.sqrt(imp[n][:,2])
    plt.errorbar(x,y,yerr=sd,fmt='o',label=('n='+str(n)),capsize=CS,color=color[i],markersize=MS,lw=LW-1)
    i+=1

plt.plot(analytical50[:,0],analytical50[:,1],label='analytical n=50',color=color[0],lw=LW)
#plt.plot(analytical100[:,0],analytical100[:,1],label='analytical n=100',color=color[1])
plt.plot(analytical200[:,0],analytical200[:,1],label='analytical n=200',color=color[1],lw=LW)


plt.plot(unStabel50[:,0],unStabel50[:,1],color=color[0],linestyle='--',lw=LW)
#plt.plot(unStabel100[:,0],unStabel100[:,1],color=color[1],linestyle='--')
plt.plot(unStabel200[:,0],unStabel200[:,1],color=color[1],linestyle='--',lw=LW)

plt.plot(zeroline200[:,0],zeroline200[:,1],color=color[1],lw=LW)
#plt.plot(zeroline100[:,0],zeroline100[:,1],color=color[1])
plt.plot(zeroline50[:,0],zeroline50[:,1],color=color[0],lw=LW)

plt.legend(frameon=False,bbox_to_anchor=(0.12,0.78),loc='lower left',fontsize=SMALL_SIZE,ncol=2)
plt.locator_params(axis='y', nbins=6)
plt.locator_params(axis='x', nbins=5)
plt.xlim([0,0.5])
plt.ylim([-0.2,1.49])


a = plt.axes([0.55, 0.4, .4, .36])
#plt.title('without molecules',fontsize=MEDIUM_SIZE)
x_spins=spins[::3,0]
y_spins=spins[::3,1]
sd_spins=np.sqrt(spins[::3,2])
plt.errorbar(x_spins,y_spins,yerr=sd_spins,fmt='o',label='n=400',capsize=CS,color=color[0],markersize=MS,lw=LW-1)
plt.plot(analytical_spins[:,0],analytical_spins[:,1],label='analytical',color=color[0],lw=LW)
plt.plot(zeroline_spins[:,0],zeroline_spins[:,1],color=color[0],lw=LW)
plt.locator_params(axis='y', nbins=2)
plt.locator_params(axis='x', nbins=3)
plt.xticks(fontsize=SMALL_SIZE)
plt.yticks(fontsize=SMALL_SIZE)
plt.xlim([0,2])
#plt.xticks([])

plt.legend(frameon=False,bbox_to_anchor=(0.31,0.53),loc='lower left',fontsize=SMALL_SIZE)
plt.savefig('output.png')
plt.show()
