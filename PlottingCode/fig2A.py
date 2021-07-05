import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from numpy import trapz

dx=7.0
D0=270000.0

markers = ['.','s','^', 'v', '*', 'X','P','o']
color=0

#(R=K_d,Ron) $R=Kd=koff/kon, Ron=kon/khop
# need strings, since str(Ron) switches to exponential at 1e-5
parameters = [(.1,'1'),(.1,'0.1'), (.1, '0.01'), (.1, '0.001'), (.1, '0.0001'),(.1,'0.00001'),(.1,'0.000001'),(.1,'0.0000001')] #

xmax=100
Dhop = 2.0 * D0/(dx*dx)
flux_scale= Dhop/2
fig,ax=plt.subplots()

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'


ax.tick_params(labelsize=14)

for (R, Son) in parameters:
    if color > 6:
        thisc=sns.color_palette()[0]
    else:
        thisc=sns.color_palette()[color]

    Ron=float(Son)
    print(Son)
    kA = 1.0/R
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
    for length in 2**np.arange(0,11,1): #to 512
    
        files = glob.glob('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale*R'+str(R)+'Ron'+Son)     

        if len(files) > 0:
            
            scale = np.max([ int(re.search('Scale(\d+)',f).group(1)) for f in files])
            if length<xmax:
                maxscale=scale
            
            print(R,Ron,length,scale)

                
            dens = np.loadtxt('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale'+str(scale)+'R'+str(R)+'Ron'+Son)

            number_particles = trapz(x = dens[:length,0], y = dens[:length,1])
            number_bound = number_particles * 1.0/(1.0 + R)
            bound.append(number_bound)

            fluxes.append(dens[int(length)+1,1])
            lengths.append(length)

            if scale>0:
                smallscale=scale-1
            else:
                smallscale=scale   # so not to crash
            dens = np.loadtxt('../Sim/Data_profile/DensityProfile_Width'+str(length)+'Scale'+str(smallscale)+'R'+str(R)+'Ron'+Son)
            fluxes2.append(dens[int(length)+1,1])

            number_particles = trapz(x = dens[:length,0], y = dens[:length,1])
            number_bound = number_particles * 1.0/(1.0 + R)
            bound2.append(number_bound)

    bound = np.array(bound)
    bound2 = np.array(bound2)
    lengths=np.array(lengths,dtype=float)
    plt.loglog(lengths,fluxes2*(lengths)/flux_scale,marker = markers[color],fillstyle='none',linewidth = 0,color = thisc) #scale-1
    
    if np.log10(kOn/khop) > 0: # 1 
        plt.loglog(lengths,fluxes*(lengths)/flux_scale,marker = markers[color],linewidth=0,color = thisc,
                   label =r'$R_{{on}}=1$'.format(np.log10(kOn/khop)))
    else:
        plt.loglog(lengths,fluxes*(lengths)/flux_scale,marker = markers[color],linewidth=0,color = thisc,
                   label =r'$R_{{on}}=10^{{{:}}}$'.format(int(np.log10(kOn/khop))))

    # PLOT EXPECTED FLUXES AT LARGE L (times L)
    def tflux(phat):    # arctan part of approx flux equation
        return np.arctan((1+2*phat)/np.sqrt(3.0))
    
    D_0 = 0.5*khop*fractionFree
    phat_high=1/pscale
    phat_low=0.0/pscale
    Lflux_inf= 2* D_0*pscale/np.sqrt(3)*(tflux(phat_high)-tflux(phat_low))
    LX=1/Ron**0.333333
    x=np.linspace(LX,100,1000)
    plt.loglog(x,Lflux_inf/flux_scale+0.0*x,":",color=thisc) # expected infinite system flux
    if Ron==1e-06:
        plt.arrow(70,Lflux_inf/flux_scale,30,0,\
                  color=thisc,ls='--',width=.00002,head_width=.0003,\
                  length_includes_head=True,head_length=8)
 
    color += 1
#    if color>5: 
 #       color=0


plt.xlim(.9,xmax)
plt.ylim(1e-3,1.2)
plt.ylabel(r'$\Phi/\Phi_0$',rotation='horizontal',size=18)
plt.xlabel(r'$L$',size=18)
plt.legend(loc = 'lower left',fontsize=11)

plt.text(.45,1.2,'A)',size=18)

# limit as pscale-> infty as black dashed
pscale=1e10
D_0 = 0.5*khop*fractionFree
phat_high=1/pscale
phat_low=0.0/pscale
Lflux_inf= 2* D_0*pscale/np.sqrt(3)*(tflux(phat_high)-tflux(phat_low))
x=np.linspace(.1,1000,1000)
plt.loglog(x,Lflux_inf/flux_scale+0.0*x,":",color='k') # expected infinite system flux

x=np.linspace(10,1000,1000)

plt.tight_layout()
plt.savefig('../Plots/fig2A.pdf', bbox_inches='tight',pad_inches=0)
