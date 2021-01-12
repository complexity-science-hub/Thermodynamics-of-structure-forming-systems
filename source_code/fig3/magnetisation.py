#script that runs the mc-Simulation for all temperature steps on multiple threads

import numpy as np
import matplotlib.pyplot as plt
import mc
from multiprocessing import Pool

def TemperatureStep(Temp):
    m1=np.zeros(StepsMean)
    for j in range(StepsMean):
        mRun=mc.MonteCarlo(
        T=Temp,J=1,h=0.0,n=50,mcSteps=10000,sud=[0.0,0.5,0.5],mol=1)
        m1[j]=np.abs(mRun)
        print(Temp,mRun,j)
    mean=np.mean(m1)
    var=np.var(m1)

    print(f'------------------STEP T= {Temp} DONE!------------------')

    return(mean,var)


def main():
    global StepsMean
    output='mol50'
    steps=48
    StepsMean=1000
    Tmin=0.01
    Tmax=0.5

    Temperature=np.linspace(Tmin,Tmax,steps)
    p=Pool()
    m=p.map(TemperatureStep,Temperature)
    p.close()
    p.join()
    m=np.array(list(m))
    transT=Temperature.reshape(steps,1)
    c=np.hstack([transT,m])
    np.savetxt(output+'.csv',c, delimiter=",")
    plt.figure()
    plt.errorbar(Temperature,m[:,0],yerr=np.sqrt(m[:,1]),fmt='o')
    plt.xlabel('Temperature ($k_B T$)')
    plt.ylabel('Magnatisation ($m$)')
    plt.savefig(output+'.png')


if __name__=='__main__':
    main()
