import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
colors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]
cm = plt.get_cmap('Set1')
red = '#FE2600'

# calculate fast/slow state index
def calc_timings(data):

    state = []
    
    index_high = []
    index_low = []

    temp_index_high = []
    temp_index_low = []

    bound = False
    start = data[0,0]
    for i, (time, step, out, number_bound, _)  in enumerate(data):
        
        if number_bound == 0:
            
            temp_index_high.append(i)
            
            if not bound and i != 0:
                temp_index_low.append(i)
                index_low.append(temp_index_low)
                temp_index_low = []
                state.append([start, time, 1, len(index_low) - 1])
                start = time
            bound=True


        else:
            temp_index_low.append(i)
            
            if bound:
                temp_index_high.append(i)
                index_high.append(temp_index_high)
                temp_index_high = []
                state.append([start, time, 0, len(index_high) - 1])
                start = time
            bound=False
            

    return index_high, index_low, state



kOff_list = [0.012122, 0.121224, 1.212245]
id_list = [1,2,3]
exp = [-6,-5,-4]


R=0.1#10.0
width=8

fig,ax = plt.subplots(3,1,sharex=True,sharey=True)
ax_twin = np.array([ax[id].twinx() for id in range(3)])

for id,kOff in enumerate(kOff_list):
    
    
    kA = 1.0/R
    kOn = kA * kOff
    dx=7.0
    D0=270000.0
    Dhop = 2.0*D0/(dx*dx)
    fractionFree = 1.0/(1.0 + kA)
    khop = Dhop/fractionFree
    Ron = kOn/khop
    
    data = np.loadtxt('../Sim/Data_Old/Longer/TimeSeries_Width'+str(width)+'R'+str(R)+'koff'+str(kOff)+'ID'+str(id_list[id]))
    time, step, out, number_bound, number_total = list(data.T)

    index_high, index_low, state = calc_timings(data)
    state = np.array(state)

    
    low_bound = np.array(list(map(lambda x: number_bound[x][1:-1].mean() , index_low)))
    high_bound = np.array(list(map(lambda x: number_bound[x][1:-1].mean() , index_high)))
    
    low_tau = np.array(list(map(lambda x: time[x][-1] - time[x][0], index_low)))
    high_tau = np.array(list(map(lambda x: time[x][-1] - time[x][0], index_high)))
    
    low_flux = np.array(list(map(lambda x: out[x][1:-1].sum() if out[x][1:-1].sum() > 0 else 0.001, index_low))) / low_tau
    high_flux = np.array(list(map(lambda x: out[x][1:-1].sum() if out[x][1:-1].sum() > 0 else 0.001, index_high))) / high_tau
    
    t = 0
    times = []
    flux = []
    bound = []
    
    for i in range(len(state)):
        if state[i,2] == 0:
            times.append(state[i,0])
            times.append(state[i,1])
            bound.append(0)
            bound.append(0)
            flux.append(high_flux[state[i,3].astype(int)])
            flux.append(high_flux[state[i,3].astype(int)])
        else:
            times.append(state[i,0])
            times.append(state[i,1])
            bound.append(low_bound[state[i,3].astype(int)])
            bound.append(low_bound[state[i,3].astype(int)])
            flux.append(low_flux[state[i,3].astype(int)])
            flux.append(low_flux[state[i,3].astype(int)])       

    ax[id].set_ylabel(r'$\Phi$',color='k',fontsize = 14.5,rotation=0)
    ax_twin[id].set_ylabel(r'$n_{\mathrm{bound}}$',color='k', fontsize = 14.5)

    ax[id].set_zorder(ax_twin[id].get_zorder()+1) # put ax in front of ax_twin[id] 
    ax[id].patch.set_visible(False) # hide the 'canvas' 

    ax[id].set_yscale('log')
    ax[id].set_ylim(0.1,10000)
    ax_twin[id].set_ylim(0,8)
    
    time_scale = (times - times[0]) * Ron
    time_scale_b = (time - times[0]) * Ron
    
    if id == 2:
        ax[id].set_xlim(0,0.01)
        ax_twin[id].set_xlim(0,0.010)
        ax[id].set_xlabel(r'$R_{\mathrm{on}}t$', fontsize = 14)
        ax[id].plot(time_scale, flux,color='k',zorder=100,linewidth = 1,linestyle = '-',marker = '')
        ax_twin[id].fill_between(time_scale_b, np.zeros(number_bound.shape), number_bound,color=red,zorder=-1, linewidth = 1)
      
    if id == 1:
        ax[id].set_xlim(0,0.010)
        ax_twin[id].set_xlim(0,0.010)
        ax[id].plot(time_scale, flux,color='k',zorder=100,linewidth = 1,linestyle = '-',marker = '')
        ax_twin[id].fill_between(time_scale_b, np.zeros(number_bound.shape), number_bound,color=red,zorder=-1, linewidth = 1)
        
    
    if id == 0:
        ax[id].set_xlim(0,0.010)
        ax_twin[id].set_xlim(0,0.010)
        ax[id].plot(time_scale, flux,color='k',zorder=100,linewidth = 1,linestyle = '-',marker = '')
        ax_twin[id].fill_between(time_scale_b, np.zeros(number_bound.shape), number_bound,color=red,zorder=-1, linewidth = 1)
        
    
    ax[id].text(0.02, 1.0005, r'$R_{\mathrm{on}} = $' + r'$10^{{{}}}$'.format(exp[id]), fontsize=14.5, horizontalalignment='left', verticalalignment='bottom',transform=ax[id].transAxes, color = 'k')

    ax[id].set_xticks([0,0.002,0.004,0.006,0.008,0.01])
    ax[id].set_xticklabels(["0","0.002","0.004","0.006","0.008","0.01"])
    
    ax_twin[id].set_xticks([0, 0.002, 0.004, 0.006,0.008,0.01])
    ax_twin[id].set_xticklabels(["0","0.002","0.004","0.006","0.008","0.01"])
    
    if id == 0:
        ax_twin[id].annotate("",
            xy=(0.0098,4.75), xycoords='data',
			xytext=(0.00949, 4.75), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = red, linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
        ax_twin[id].annotate("",
            xy=(-0.00005,4.5), xycoords='data',
            xytext=(0.0005, 4.5), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = 'k', linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
    
    if id == 1:
        ax_twin[id].annotate("",
            xy=(0.00965,2.5), xycoords='data',
            xytext=(0.0091, 2.5), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = red, linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
        ax_twin[id].annotate("",
            xy=(-0.00005,4.5), xycoords='data',
            xytext=(0.00035, 4.5), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = 'k', linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
    
    if id == 2:
        ax_twin[id].annotate("",
            xy=(0.0098,4), xycoords='data',
            xytext=(0.0093, 4), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = red, linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
        ax_twin[id].annotate("",
            xy=(-0.00005,4.5), xycoords='data',
            xytext=(0.0003, 4.5), textcoords='data',
            arrowprops=dict(arrowstyle="->", color = 'k', linewidth = 2,
                            connectionstyle="arc3"), transform=ax_twin[id].transAxes, zorder=-1,annotation_clip=False
            )
    
plt.tight_layout()
plt.subplots_adjust(hspace=0.32)
fig.subplots_adjust(top=0.94,bottom=0.09)
plt.savefig('../Plots/fig3_binding_same_figure_kA'+str(int(kA))+'_L'+str(width)+'.pdf')
