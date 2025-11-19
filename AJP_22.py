# AJP_22.py

# This Python file uses the Integrated Oscillator Model (IOM) to implement
# a negative feedback loop with a delay as a mechanism to synchronize islets. 

# Used to make all model simulation figures in:
#    Coordination of Pancreatic Islet Rhythmic Activity by Delayed Negative Feed
back,
#    N. Bruce, I.-A. Wei, W. Leng, Y. Oh, Y.-C. Chiu, M.G. Roper, R. Bertram
#    American Journal of Physiology, 323:E492-E502, 2022.

# Variables:
#    v -- voltage (mV)
#    n -- activation of delayed rectifier
#    c -- free cytosolic calcium concentration (muM)
#    adp -- cytosolic ADP concentration (muM)
#    cer -- concentration of free calcium in the endoplasmic reticulum (muM)
#    cam -- mitocondrial calcium concentration (muM)
#    f6p -- fructose 6-phosphate concentration (muM)
#    fbp -- fructose 1,6-bisphosphate concentration (muM)
#    I -- insulin, arbitrary units
#    Gi -- intracellular glucose (mM)
#    Ge -- extracellular glucose (mM)
#
# Units:
# time: ms
# volume: l
# conductance: pS
# V: mV
# I: fA
# J: muM/ms (mM/ms for Jgk and Jglut)
# alpha: mumol/fA/ms

import numpy as np
import matplotlib.pyplot as plt
from jitcdde import jitcdde, y, t
from symengine import exp, sqrt
import random
from jitcxde_common import conditional

# Number of cells
N=5

#PARAMETERS

# Electrical and calcium model
# Ca current ica
gca=1000
vca=25
vm=-20
sm=12

# Delayed-rectifying K+ current, ik
gk=2700
vk=-75
taun=20
nin=-16
sn=5

# Ca-activated K+ current, ikca
kd=0.5
gkca=50

# K-ATP channel current, ikatp
kdd=17
ktd=26
ktt=1
atot=3000
gkatpbar=27500

# Ca fluxes across the plasma membrane
alpha=5.18e-18
vcyt=1.15e-12
kpmca=0.2

# Ca fluxes into ER
kserca=0.4
pleak=2.0e-4

# Ca fluxes into the mitochondria
kuni=0.4
knaca=0.001

# The vector field
Cm=5300
fca=0.01
sigmam=100
sigmaer=31

# GLUT-2 reaction rate, Jglut
Vglut=8
Kglut=7

# Glucokinase reaction rate, Jgk
Kgk=11
ngk=15
Vgk=[0.0014,0.0015,0.0025,0.0018,0.0029]

# PFK reaction rate, Jpfk
k1=30
k2=1
k3=50000
k4=1000
famp=0.02
fatp=20
ffbp=0.2
fbt=20
fmt=20
kpfk=0.06
vpfk=0.01

# PDH reaction rate, Jpdh
kCaPDH=200
vpdh=0.009

# Extracellular glucose, Ge
tauG=50000
Gmin=7
Gmax=13
Ihat=15
SG=1

# Insulin equilibrium function, I
taui=10000
ki=0.1
Islope=1000
delta=8

# ADP equation
taua=300000

# Delay parameter
taud=0

# EQUATIONS

def f():
    Ge=y(10*(N-1)+10)
    for i in range(N):
        v=y(10*i)
        n=y(10*i+1)
        c=y(10*i+2)
        cer=y(10*i+3)
        cam=y(10*i+4)
        adp=y(10*i+5)
        f6p=y(10*i+6)
        fbp=y(10*i+7)
        Gi=y(10*i+8)
        I=y(10*i+9)
        
        # Ca current ica
        minf=1/(1+exp((vm-v)/sm))
        ica=gca*minf*(v-vca)
        
        # Delayed-rectifying K+ current, ik
        ninf=1/(1+exp((nin-v)/sn))
        ik=gk*n*(v-vk)
        
        # Ca-activated K+ current, ikca
        qinf=c**2/(kd**2+c**2)
        ikca=-gkca*qinf*(vk-v) 
        
        # K-ATP channel current, ikatp
        rad=sqrt(-4*adp**2+(atot-adp)**2)
        atp=(atot+rad-adp)/2
        mgadp=0.165*adp
        adp3m=0.135*adp
        atp4m=0.05*atp
        topo=0.08+0.89*mgadp**2/kdd**2+0.16*mgadp/kdd 
        bottomo=(1+mgadp/kdd)**2*(1+atp4m/ktt +adp3m/ktd) 
        katpo=topo/bottomo
        ikatp=gkatpbar*katpo*(v-vk)
        
        # Ca fluxes across the plasma membrane
        Jmem=-(alpha/vcyt*ica + kpmca*c)
        
        # Ca fluxes into ER
        Jer=kserca*c - pleak*(cer-c) 
        
        # Ca fluxes into the mitochondria
        Jm= kuni*c - knaca*(cam-c) 
        
        # GLUT-2 reaction rate, Jglut
        Jglut=Vglut*((Ge-Gi)*Kglut)/((Kglut+Ge)*(Kglut+Gi))
        
        # Glucokinase reaction rate, Jgk
        Jgk=Vgk[i]*(Gi)**ngk/(Kgk**ngk+Gi**ngk)
        
        # Insulin equilibrium function, I
        Iinf=conditional(c,ki,0,1)*Islope*(c-ki)
        
        # PFK reaction rate, Jpfk
        amp=adp**2/atp
        
        weight1=1
        topa1=0
        bottom1=1
        
        weight2=atp**2/k4
        topa2=topa1
        bottom2=bottom1+weight2
        
        # (0,0,1,0)
        weight3=f6p**2/k3
        topa3=topa2+weight3
        bottom3=bottom2+weight3
        
        # (0,0,1,1)
        weight4=(f6p*atp)**2/(fatp*k3*k4)
        topa4=topa3+weight4
        bottom4=bottom3+weight4
        
        # (0,1,0,0)
        weight5=fbp/k2
        topa5=topa4
        bottom5=bottom4+weight5
        
        # (0,1,0,1)
        weight6=(fbp*atp**2)/(k2*k4*fbt)
        topa6=topa5
        bottom6=bottom5+weight6
        
        # (0,1,1,0)
        weight7=(fbp*f6p**2)/(k2*k3*ffbp)
        topa7=topa6+weight7
        bottom7=bottom6+weight7
        
        # (0,1,1,1)
        weight8=(fbp*f6p**2*atp**2)/(k2*k3*k4*ffbp*fbt*fatp)
        topa8=topa7+weight8
        bottom8=bottom7+weight8
        
        # (1,0,0,0)
        weight9=amp/k1
        topa9=topa8
        bottom9=bottom8+weight9
        
        # (1,0,0,1)
        weight10=(amp*atp**2)/(k1*k4*fmt)
        topa10=topa9
        bottom10=bottom9+weight10
        
        # (1,0,1,0)
        weight11=(amp*f6p**2)/(k1*k3*famp)
        topa11=topa10+weight11
        bottom11=bottom10+weight11
        
        # (1,0,1,1)
        weight12=(amp*f6p**2*atp**2)/(k1*k3*k4*famp*fmt*fatp)
        topa12=topa11+weight12
        bottom12=bottom11+weight12
        
        # (1,1,0,0)
        weight13=(amp*fbp)/(k1*k2)
        topa13=topa12
        bottom13=bottom12+weight13
        
        # (1,1,0,1)
        weight14=(amp*fbp*atp**2)/(k1*k2*k4*fbt*fmt)
        topa14=topa13
        bottom14=bottom13+weight14
        
        # (1,1,1,0) -- the most active state of the enzyme
        weight15=(amp*fbp*f6p**2)/(k1*k2*k3*ffbp*famp)
        topa15=topa14
        topb=weight15
        bottom15=bottom14+weight15
        
        # (1,1,1,1)
        weight16=(amp*fbp*f6p**2*atp**2)/(k1*k2*k3*k4*ffbp*famp*fbt*fmt*fatp)
        topa16=topa15+weight16
        bottom16=bottom15+weight16
        
        Jpfk= vpfk*(topb + kpfk*topa16)/bottom16
        
        # PDH reaction rate, Jpdh
        sinfty= cam/(cam+kCaPDH)
        Jpdh=vpdh*sinfty*sqrt(fbp)
        
        # Differential Equations
        yield -(ica + ik + ikca + ikatp)/Cm
        yield -(n-ninf)/taun
        yield fca*(Jmem - Jm - Jer)
        yield fca*sigmaer*Jer
        yield fca*sigmam*Jm
        yield (atp-exp((1+2.2*Jpdh/(Jpdh+0.05))*(1-c/0.35))*adp)/taua
        yield 0.3*(Jgk-Jpfk)
        yield (Jpfk-Jpdh/2)
        yield Jglut-Jgk
        yield (Iinf-I)/taui
        
        # Delay implementation
    Iavg=0    
    for j in range(N):
        Iavg=Iavg+y(10*j+9,t-taud)
    Iavg=Iavg/N
    
    # Ge equilibrium function
    Ginf=Gmin+(Gmax-Gmin)/(1+exp((Iavg-Ihat)/SG))    
    yield conditional(t,2400000,0,1)*(Ginf-Ge)/tauG

DDE = jitcdde(f)

# Initial Conditions
size=10*N+1
x0=[0]*size
for i in range(N):
        x0[10*i]=random.uniform(0,35)-60
        x0[10*i+1]=0
        x0[10*i+2]=random.uniform(0,0.1)+0.1
        x0[10*i+3]=random.uniform(0,200)+200
        x0[10*i+4]=100
        x0[10*i+5]=random.uniform(0,60)+840
        x0[10*i+6]=random.uniform(0,30)+70
        x0[10*i+7]=random.uniform(0,40)
        x0[10*i+8]=10
        x0[10*i+9]=random.uniform(0,20)
x0[10*(N-1)+10]=10  
DDE.constant_past(x0)

# Numerical details 
DDE.step_on_discontinuities()
dt=int(7200000/20)
times=np.linspace(DDE.t,7200000,dt)
result=np.vstack([DDE.integrate(time) for time in times])

# Plotting 
ge=result[:,10*N]
tmin=times/60000
cavg=0   
c=[0]*N
for j in range(N):
    cavg=cavg+result[:,10*j+2]
cavg=(cavg/N)*1000
fig,ax = plt.subplots()
ax.plot(tmin, cavg)
plt.ylim([0, 500])
ax.set_xlabel("Time (min)",fontsize=14)
ax.set_ylabel("$[Ca^{2+}]$ (nM)",color="blue",fontsize=14)
ax2=ax.twinx()
ax2.plot(tmin, ge,color="red")
plt.ylim([0, 14])
ax2.set_ylabel("Glucose (nM)",color="red",fontsize=14)
plt.show()
