import numpy as np
import matplotlib.pyplot as plt
import glob
import seaborn as sns
import re
# FIG 2B, all with Ron=1e-7

markers = ['.','s','^', 'v', '*', 'X','P','D','o']
color=0

dx=7.0
D0=270000.0
xmax=100
Dhop = 2.0 * D0/(dx*dx)
flux_scale= Dhop/2
fig,ax=plt.subplots()


import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'

ax.tick_params(labelsize=14)


#(R=K_d,Ron) $R=Kd=koff/kon, Ron=kon/khop   # need strings, since str(Ron) switches to exponential at 1e-5
parameters = [(100.0,'0.0000001',.55),(10.0,'0.0000001',.55),(1.0,'0.0000001',.55),(.1,'0.0000001',.55),(.01,'0.0000001',.55)]

for (R, Son,factor) in parameters:
    Ron=float(Son)
    print(Son)
    kA = 1.0/R
    kD = R

    fractionFree = 1.0/(1.0+kA)
    khop = Dhop/fractionFree

    kOn = Ron*khop
    kOff= R*kOn
    pscale = np.sqrt(kOff/Dhop) * (1.0 + 1.0/kA)
    fluxes = []
    fluxes2 = []
    bound = []
    bound2 = []
    lengths = []

    fluxes.append(Dhop/2.0)
    fluxes2.append(Dhop/2.0)
    lengths.append(1)
    bound.append(1.0/(1.0 + R))
    bound2.append(1.0/(1.0 + R))
    for length in 2**np.arange(0,10,1): #to 64
    
        files = glob.glob('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale*R'+str(R)+'Ron'+Son)
        if len(files) > 0:
	        
            scale = np.max([ int(re.search('Scale(\d+)',f).group(1)) for f in files])
            if length<xmax:
                maxscale=scale
            
            print(R,Ron,length,scale)
              
            dens = np.loadtxt('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale'+str(scale)+'R'+str(R)+'Ron'+Son)
            fluxes.append(dens[int(length)+1,1])
            lengths.append(length)

            dens = np.loadtxt('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale'+str(scale-1)+'R'+str(R)+'Ron'+Son)
            fluxes2.append(dens[int(length)+1,1])

    lengths=np.array(lengths,dtype=float)
    print(lengths,fluxes2*(lengths)/flux_scale)
    plt.loglog(lengths,fluxes2*(lengths)/flux_scale,marker = markers[color],fillstyle='none',linewidth = 0,color = sns.color_palette()[color]) #scale-1
    plt.loglog(lengths,fluxes*(lengths)/flux_scale,marker = markers[color],linewidth=0,color = sns.color_palette()[color],label = '$K_A=$'+'{:}'.format(kA))
    # both empty and full, to show proper average 

    # approximate L-> infty flux as the slow component
    def tflux(phat):    # arctan part of approx flux equation
        return np.arctan((1+2*phat)/np.sqrt(3.0))
    D_0 = 0.5*khop*fractionFree
    phat_high=1/pscale
    phat_low=0.0/pscale
    flux_inf= 2*pscale/np.sqrt(3)*(tflux(phat_high)-tflux(phat_low))

#   two state model curve
    x=np.linspace(1,xmax,1000)
    tauf=1/(kOn*(x+1)/2)

    nbound=1+(x-1)/(2*(1+kD))*factor   # factor from parameters list above, smaller factor moves flux curve up
    taus=nbound/kOff    
    flux_slow=(x+1)/2/taus

    y=(1+kA)/(1+taus/tauf)  # no difference for smaller L
    plt.loglog(x,y,":",color=sns.color_palette()[color])

    color += 1
    if color>5: 
        color=0
        
# Horner 2018 Figure 8 data
lengths=np.array([6,12,12,12,20,22,26,30])
fluxes=np.array([216.88,52.56,28.91,18.23,37.35,23.55,6.01,1.70])   # length0=302.58
plt.loglog(lengths,0.0002*fluxes*lengths,color='k',marker='o',fillstyle='none',linewidth=0,label='Horner 2018',markeredgewidth=1.5,markersize=7)
# with arbitrary prefactor, and scaled by length




plt.xlim(.9,xmax)
plt.ylim(1e-3,1.2)
plt.ylabel(r'$\Phi/\Phi_0$',rotation='horizontal',size=18)
plt.xlabel(r'$L$',size=18)
plt.legend(loc = 'lower left',fontsize=11)

plt.text(.45,1.2,'B)',size=18)

# limit as pscale-> infty as black dotted horizontal 
x=np.linspace(1,1000,1000)
plt.loglog(x,1.0+0.0*x,":",color='k') # expected infinite system flux


# black dashed line
x=np.linspace(5,1000,1000)
plt.loglog(x,1/(x*x),"--",color='k') 

# save figs
plt.tight_layout()
plt.savefig('../Plots/fig2B.pdf', bbox_inches='tight',pad_inches=0)
