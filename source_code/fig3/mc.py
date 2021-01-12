import numpy as np
import random
import matplotlib.pyplot as plt

#returns the fraction of n1! and n2! with the bigger one as the numerator
def dfac(n1,n2):
    if n1>n2:
        b=n1
        s=n2
    else:
        b=n2
        s=n1
    if b==s:
        return(1)
    elif b==s+1:
        return(b)
    elif b==s+2:
        return(b*(b-1))
    else:
        print('error')
        print(n1,n2)
        return(1)

#return True if value is not 0
def isNotZero(x):
    if x==0:
        return(False)
    else:
        return(True)

#calculates full hamiltonian
def hamiltonian(Spins,J,h,n):
    nud=len(Spins)
    H=0
    if h!=0:
        for i in range(nud):
            H-=h*Spins[i]
            for j in range(i):
                H-=J*Spins[i]*Spins[j]/(n)
    else:
        for i in range(nud):
            for j in range(i):
                H-=J*Spins[i]*Spins[j]/(n)
    return(H)

#more efficent hamiltonian calculation for spin flips
def hamiltonianChangeSpinFlip(OldH,i,Spins,J,h,n):
    nud=len(Spins)
    if h!=0:
        H=hamiltonian(Spins,J,h,n)
    else:
        tempH=0
        for j in range(nud):
            tempH+=(Spins[j]*Spins[i])
        H=OldH-2*(J/n)*(tempH-1)
    return(H)

#initializes the system with the distribution pSpin
def initMC(ParticleNumber,pSpins=[0,0.5,0.5]):
    if sum(pSpins)!=1:
        print('error')
        return('error')
    Spins=np.zeros([2,ParticleNumber])
    i=0
    while (i < ParticleNumber):
        rnd=random.random()
        if i>ParticleNumber-2:
#            print("Last Spin random")
            if rnd<pSpins[1]:
                Spins[0,i]=1
                i+=1
            else:
                Spins[0,i]=-1
                i+=1
        elif (rnd <pSpins[0]*0.5):
            Spins[1,i]=2
            i+=2
        elif (rnd>=pSpins[0]*0.75+pSpins[1]):
            Spins[0,i]=1
            i+=1
        elif (rnd>=pSpins[0]*0.5) and (rnd<pSpins[0]*0.75+pSpins[1]):
            Spins[0,i]=-1
            i+=1
    return(Spins)

#tries a Spin flip in the simulation
def SpinFlip(OldH,Spins,beta,J,h,n):
    NSpins=len(Spins)
    if (NSpins==0):
        return(OldH,Spins)
    i=np.random.choice(NSpins)
    temp=np.copy(Spins)
    temp[i]=-temp[i]
    rnd=random.random()
    NewH=hamiltonianChangeSpinFlip(OldH,i,temp,J,h,n)
    if rnd<min(1,np.exp(-beta*(NewH-OldH))):
        return(NewH,temp)
    else:
        return(OldH,Spins)

#does a Molecule move in the simulation
def MoleculeMove(OldH,Spins,Molecules,beta,J,h):
    Spins=np.array(Spins)
    NewSpins=np.copy(Spins)
    ns=len(Molecules)
    nud=len(Spins)
    n=2*ns+nud
    if nud!=0:
        Nu=(Spins==1).sum()
        Nd=(Spins==-1).sum()

        p=[Nu/nud,Nd/nud]
    else:
        p=[0.5,0.5]
        Nu=0
        Nd=0
    pdiv=0.5
    rnd1=random.random()
    rnd2=random.random()
    if rnd1<pdiv:
        if ns==0:
            return(OldH,Spins,Molecules)
        RndSpins=np.random.choice([1,-1], 2, p)
        NewSpins=np.append(NewSpins,RndSpins)
        NewH=hamiltonian(NewSpins,J,h,n)
        NewNu=(NewSpins==1).sum()
        NewNd=(NewSpins==-1).sum()
        if (NewNu+NewNd!=Nu+Nd+2):
            print(f'ERROR DIV{NewNu,NewNd,Nu,Nd}')

        pacc=((2*ns)/(dfac(NewNu,Nu)*dfac(NewNd,Nd)))*np.exp(-beta*(NewH-OldH))
#        pacc=((2*ns)/((nud+1)*(nud+2)))*np.exp(-beta*(NewH-OldH))
        if rnd2<min(1,pacc):
            NewMolecules=Molecules[1:]
            pick=np.random.choice(nud+2, 2, replace=False)
            NewSpins[pick[0]],NewSpins[-1]=NewSpins[-1],NewSpins[pick[0]]
            NewSpins[pick[1]],NewSpins[-2]=NewSpins[-2],NewSpins[pick[1]]
            return(NewH,NewSpins,NewMolecules)
        else:
            return(OldH,Spins,Molecules)
    else:
        if 2*ns>=n-1:
            return(OldH,Spins,Molecules)
        pick=np.random.choice(nud, 2, replace=False)
        NewSpins=np.delete(Spins,[pick[0],pick[1]])
        NewMolecules=np.append(Molecules,[2])
        NewH=hamiltonian(NewSpins,J,h,n)
        NewNu=(NewSpins==1).sum()
        NewNd=(NewSpins==-1).sum()
        if (NewNu+NewNd!=Nu+Nd-2):
            print(f'ERROR COM{NewNu,NewNd,Nu,Nd}')
            print(NewSpins)

        pacc=((dfac(NewNu,Nu)*dfac(NewNd,Nd))/(2*ns+2))*np.exp(-beta*(NewH-OldH))
#        pacc=((nud*(nud-1))/(2*ns+2))*np.exp(-beta*(NewH-OldH))
        if rnd2<min(1,pacc):
            NewMolecules=np.append(Molecules,[2])
            return(NewH,NewSpins,NewMolecules)
        else:
            return(OldH,Spins,Molecules)

#calculates the Distribution
def calcDistribution(t,Spins,Molecules,OldDist):
    OldDist[1,t]=np.count_nonzero(Spins==1)
    OldDist[2,t]=np.count_nonzero(Spins==-1)
    OldDist[0,t]=np.sum(Molecules)
    return(OldDist)

#plot the distribution of the 3 states
def plotDist(DistOfStates,n):
    D=DistOfStates/n
    plt.figure()
    plt.plot(D[0,1:],'r')
    plt.plot(D[1,1:],'b')
    plt.plot(D[2,1:],'y')
    plt.savefig('dist.png')

def MonteCarlo(T=100,J=2,h=0.0,n=50,mcSteps=2500,sud=[0,0.5,0.5],mol=1):
    mol=mol+1
    beta=1.0/T
    DistOfStates=np.zeros([3,mcSteps])
    FullSpins=initMC(n,sud)
    Spins=list(filter(isNotZero,FullSpins[0,:]))
    Molecules=list(filter(isNotZero,FullSpins[1,:]))
    H=hamiltonian(Spins,J=J,h=h,n=n)
    for t in range(mcSteps):
#        DistOfStates=calcDistribution(t,Spins,Molecules,DistOfStates)
        if t%mol==0:
            H,Spins=SpinFlip(H,Spins,beta,J,h,n)
        else:
            H,Spins,Molecules=MoleculeMove(H,Spins,Molecules,beta,J,h)

#    plotDist(DistOfStates,n)
    m=(1.0/n)*np.sum(Spins)
#    print(T,sum(Molecules),H,m)
    return(m)
