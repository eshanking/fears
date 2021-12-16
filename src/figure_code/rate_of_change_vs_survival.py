from fears.utils import results_manager, plotter
import os
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.ticker import ScalarFormatter

fig,ax = plt.subplots(figsize=(3.5,7.75))
labelsize=12

data_folder = 'results_07212021_0001'
exp_info_file = 'experiment_info_07212021_0001.p'
exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
k_abs = exp_info.slopes


# make a dummy barchart
# k_abs = k_abs[2:]
x = np.arange(len(k_abs))
barchart_data = np.ones(len(k_abs))*100
rects = ax.barh(x,barchart_data,color='slategrey',facecolor='w')

tc_axes = []
drug_axes = []

exp_folders.reverse()
k_abs = np.flip(k_abs)
for exp in exp_folders:
    k_abs_t = exp[exp.find('=')+1:]
    k_abs_t = float(k_abs_t)
    num = np.argwhere(k_abs == k_abs_t)
    num = num[0,0]
    
    # generate timecourse axes
    
    width = 110
    height = rects.patches[num].get_height()
    ypos = rects.patches[num].get_y()
    xpos = 120
    
    tcax = ax.inset_axes([xpos,ypos,width,height],transform=ax.transData)
    # da = tcax.twinx()
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    survive_count = 0
    counts_total = None
    
    k=0
    while k < len(sim_files):
    # for sim in sim_files:
        sim = sim_files[k]
        sim = exp + os.sep + sim
        data = results_manager.get_data(sim)
        dc = data[:,-1]
        data = data[:,0:-1]
        # data = data/np.max(data)
        data_t = data[-1,:]
        
        # check to see if any genotypes are at least 10% of the max cell count
        if any(data_t >= 1):
            survive_count += 1
            if counts_total is None:
                counts_total = data
            else:
                counts_total += data
        # data = data/np.max(data)
        
        # exp_info.populations[num].counts_log_scale = True
        data = data/max_cells
        if k==0:
            drug_kwargs = {'alpha':0.7,
                           'color':'black',
                           'linewidth':2,
                           'label':'Drug Concentration ($\u03BC$M)'
                           }
            tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tcax,
                                                        drug_curve=dc,
                                                        drug_ax_sci_notation=True,
                                                        drug_kwargs=drug_kwargs,
                                                        legend_labels=False,
                                                        grayscale=True,
                                                        color='gray', 
                                                        linewidth=1,
                                                        labelsize=12,
                                                        alpha=0.7
                                                        )
            drug_ax.set_ylabel('')
            drug_axes.append( drug_ax )
        else:
            tcax,da = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tcax,
                                                        grayscale=True,
                                                        color='gray',
                                                        legend_labels=False,
                                                        linewidth=2,
                                                        labelsize=12,
                                                        alpha=0.2
                                                        )            
        # drug_ax.set_ylim(0,10**4)
        k+=1
        
    if survive_count > 0:
        counts_avg = counts_total/survive_count
        # counts_avg = counts_avg/np.max(counts_avg) 
        # counts_avg = counts_total
        counts_avg = counts_avg/np.max(counts_avg)
        tcax,temp = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                               counts_avg,
                                               tcax,
                                               labelsize=12)
    
    # t = np.arange(len(dc))
    # t = t*exp_info.populations[0].timestep_scale/24
    # da.plot(t,dc)
    
    tc_axes.append( tcax )
    barchart_data[num] = survive_count   

# for a in tc_axes:
    # a.set_yscale('log',base=2)
    # a.set_ylim(10,max_cells)

for da in drug_axes:
    da.ticklabel_format(style='sci',axis='y',scilimits=(0,4))
    da.yaxis.tick_right()
    # da.set_yticks(da.get_yticks())
    # yt = da.get_yticks
    # yt = yt/10**3

drug_axes[1].set_ylabel('Drug Concentration (uM)', color='gray',fontsize=labelsize)
tc_axes[1].set_ylabel('Proportion of \nmax cell count',fontsize=labelsize)
tc_axes[0].set_xlabel('Days',fontsize=labelsize)

rects = ax.barh(x,barchart_data,color='slategrey')
ax.set_yticks(x)
ax.set_yticklabels(k_abs*10**3)
# ax.yaxis.set_major_formatter(ScalarFormatter())
# ax.ticklabel_format(style='sci',axis='y')
ax.set_xlabel('% Resistant',fontsize=12)
ax.set_ylabel('$k_{abs} (x10^{-3})$',fontsize=12)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

# compute error bars
# rule of succession explanation: https://en.wikipedia.org/wiki/Rule_of_succession
p = (np.array(barchart_data) + 1)/(n_sims + 2) # uniform prior (rule of succession) 
n = n_sims
q = 1-p

sd = 100*(p*q/n)**0.5 # standard deviation of the estimator of the parameter of a bernoulli distribution

# rects = bar_ax.barh(x, barchart_data,color='slategrey')
# errorbar_pos = x + rects[0].get_height()/2
# bar_ax.errorbar(x=barchart_data, y=errorbar_pos, xerr=sd,linewidth=0,elinewidth=2,capsize=5,color='tomato')
ax.errorbar(x=barchart_data, y=x, xerr=sd,linewidth=0,elinewidth=2,capsize=5,color='black')

tc_axes[0].legend(frameon=False,fontsize=11,loc='lower left',
                  bbox_to_anchor=(-0.8,-0.95),ncol=4)
drug_axes[0].legend(frameon=False,fontsize=11,loc='lower left',
                    bbox_to_anchor=(-0.8,-1.15),ncol=1)

ax.annotate('$C(t) = C_{max}(1-e^{-k_{abs}*t})$',xy=(30,0)) # xy in data points

results_manager.save_fig(fig,'rate_of_change_v_survival.pdf',bbox_inches='tight')
