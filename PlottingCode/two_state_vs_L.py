import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic, sem
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
colors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]
import sys

# compute fast/slow state index
def calc_timings(data):

    index_high = []
    index_low = []

    temp_index_high = []
    temp_index_low = []
    
    num_bound_high = []

    bound = False
    
    for i, (time, step, out, number_bound, number_tot)  in enumerate(data):
        
        num_bound_high.append(number_tot)
        if number_bound == 0:
            
            temp_index_high.append(i)
            if not bound and i != 0:
                temp_index_low.append(i)
                index_low.append(temp_index_low)
                temp_index_low = []
            bound=True

        else:
            
            temp_index_low.append(i)
            if bound:
                temp_index_high.append(i)
                index_high.append(temp_index_high)
                temp_index_high = []
            bound=False
        
    if bound:
        temp_index_high.append(i)
        index_high.append(temp_index_high)
        temp_index_high = []
    else:
        temp_index_low.append(i)
        index_low.append(temp_index_low)
        temp_index_low = []
    return index_high, index_low, num_bound_high


kOff=0.121224
R=0.1

kA = 1.0/R
kOn = kA * kOff
dx=7.0
D0=270000.0
Dhop = 2.0*D0/(dx*dx)
fractionFree = 1.0/(1.0 + kA)
khop = Dhop/fractionFree
Ron = kOn/khop
pscale = np.sqrt((kOn + kOff)/khop) * (1.0 + 1.0/kA)

tau_high_list = []
tau_low_list = []
flux_high_list = []
flux_low_list = []
prob_list = []
n_high_list = []
n_low_list = []

flux_eqn1 = []
flux_eqn2 = []
phi_slow = []
phi_free = []
tau_slow = []
tau_free = []


width_list = [2,4,8,16,32,64,128,256,512,1024]

for id,width in enumerate(width_list):

    data = np.loadtxt('../Sim/Data/TimeSeries_Width'+str(width)+'R'+str(R)+'koff'+str(kOff)+'ID'+str(id))
    
    time, step, out, number_bound, number_total = list(data.T)
    index_high,index_low, num_bound_high = calc_timings(data)
    
    total_time = time[-1] - time[0]
    
    nu_s = 0.55
    n_plug = (width-1)/2.0
    n_bound = 1 + nu_s * n_plug/(1 + 1.0/kA)
    phi_slow.append( (1.0 + n_plug) / (n_bound/kOff) )
    phi_free.append(khop/2.0/width)
    
    n_free = (width+1)/2.0
    tau_slow.append(n_bound/kOff)
    tau_free.append(1.0/(kOn * n_free))
    
    flux_eqn1.append(2.0*D0/(dx**2)*pscale/width/np.sqrt(3) * (np.arctan(1 + 2.0/pscale)  - np.arctan(1) )  )
    
    tau_high = np.array(list(map(lambda x: time[x][-1] - time[x][0], index_high)))
    tau_low = np.array(list(map(lambda x: time[x][-1] - time[x][0], index_low)))

    n_high = np.array(list(map(lambda x: out[x].sum(), index_high)))
    n_low = np.array(list(map(lambda x: out[x].sum(), index_low)))

    n_high = n_high[tau_high > 0]
    n_low = n_low[tau_low > 0]

    tau_high = tau_high[tau_high > 0]
    tau_low = tau_low[tau_low > 0]

    flux_high = (n_high / tau_high)
    flux_low = (n_low / tau_low)

    # duration weighted fluxes
    flux_high_av = np.nansum(flux_high * tau_high) / np.nansum(tau_high)
    flux_low_av = np.nansum(flux_low * tau_low) / np.nansum(tau_low)
    
    flux_high_list.append(flux_high_av)
    flux_low_list.append(flux_low_av)

    tau_high_av = np.mean(tau_high)
    tau_low_av = np.mean(tau_low)
    
    tau_high_list.append(tau_high_av)
    tau_low_list.append(tau_low_av)
    
    total_time_high = np.sum(tau_high/total_time)
    total_time_low = np.sum(tau_low/total_time)
    prob_list.append(total_time_high/(total_time_low + total_time_high))
    
    n_high_av = np.mean(n_high)
    n_low_av = np.mean(n_low)
    n_high_list.append(n_high_av)
    n_low_list.append(n_low_av)

    
fig,ax = plt.subplots(1,2,figsize=(9.75,4))
ax = ax.flatten()

width_list = np.array(width_list)
tau_high_list = np.array(tau_high_list)
tau_low_list = np.array(tau_low_list)
flux_high_list = np.array(flux_high_list)
flux_low_list = np.array(flux_low_list)

Lx = Ron**(-1.0/3.0)

ax[0].loglog(width_list[width_list < Lx], tau_low_list[width_list < Lx], marker='s',color = colors[1], label = '',markersize = 8, linestyle = '')
ax[0].loglog(width_list[width_list < Lx], tau_high_list[width_list < Lx], marker='o',color = colors[0],label = '',markersize = 8, linestyle = '')

ax[0].loglog(width_list[width_list >= Lx], tau_low_list[width_list >= Lx], marker='s',color = colors[1], label = 'slow/plugged',markersize = 8, linestyle = '')
ax[0].loglog(width_list[width_list >= Lx], tau_high_list[width_list >= Lx], marker='o',color = colors[0],label = 'fast/flowing',markersize = 8,linestyle = '')


ax[0].set_ylabel(r'$\langle\tau\rangle$',fontsize=20,rotation=0,labelpad=7)
ax[0].set_xlabel(r'$L$',fontsize=20)
ax[0].legend(loc = 'upper left',fontsize = 16)

ax[1].loglog(width_list[width_list < Lx], flux_high_list[width_list < Lx],marker='o',color = colors[0],label = 'high',markersize = 8, linestyle = '')
ax[1].loglog(width_list[width_list < Lx], flux_low_list[width_list < Lx],marker = 's',color = colors[1], label = 'low',markersize = 8, linestyle = '')

ax[1].loglog(width_list[width_list >= Lx], flux_high_list[width_list >= Lx],marker='o',color = colors[0],label = 'high',markersize = 8,linestyle = '')
ax[1].loglog(width_list[width_list >= Lx], flux_low_list[width_list >= Lx],marker = 's',color = colors[1], label = 'low',markersize = 8,linestyle = '')

ax[1].set_ylabel(r'$\langle \Phi\rangle$',fontsize=20,rotation=0,labelpad=7)
ax[1].set_xlabel(r'$L$',fontsize=20)

flux_eqn1 = np.array(flux_eqn1)
phi_slow = np.array(phi_slow)
phi_free = np.array(phi_free)
tau_slow = np.array(tau_slow)
tau_free = np.array(tau_free)

ax[1].loglog(width_list[width_list>100], flux_eqn1[width_list>100],marker = '',color = 'k',linestyle = '--')

ax[1].loglog(width_list[width_list< Lx], phi_slow[width_list < Lx],marker = '',color = colors[1],linestyle = '--')
ax[1].loglog(width_list[width_list< Lx], phi_free[width_list < Lx],marker = '',color = colors[0],linestyle = '--')
ax[0].loglog(width_list[width_list< Lx], tau_slow[width_list < Lx],marker = '',color = colors[1],linestyle = '--')
ax[0].loglog(width_list[width_list< Lx], tau_free[width_list < Lx],marker = '',color = colors[0],linestyle = '--')

ax[1].text(0.75, 0.066, r'$\Phi_{Fick}$', fontsize=17, horizontalalignment='left', verticalalignment='center',transform=ax[1].transAxes, color = 'k')

ax[0].text(0.95, 0.97, 'A)', fontsize=20, horizontalalignment='right', verticalalignment='top',transform=ax[0].transAxes, color = 'k')
ax[1].text(0.95, 0.97, 'B)', fontsize=20, horizontalalignment='right', verticalalignment='top',transform=ax[1].transAxes, color = 'k')


ax[0].tick_params(labelsize=15)
ax[1].tick_params(labelsize=15)

ax[0].set_ylim(3E-3,5E5)
ax[1].set_ylim(5E-3,6E4)
ax[0].set_xlim(1,1150)
ax[1].set_xlim(1,1150)

plt.tight_layout()
plt.subplots_adjust(wspace=0.24)
plt.savefig('../Plots/fig4AB_kA'+str(int(kA))+'_Ron'+str(round(Ron,9))+'.pdf')
