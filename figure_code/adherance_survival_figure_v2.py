# https://matplotlib.org/api/_as_gen/matplotlib.gridspec.GridSpec.html#matplotlib.gridspec.GridSpec
import os
import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter

data_folder = 'results_07202021_0000'
exp_info_file = 'experiment_info_07202021_0000.p'

fig,ax = plt.subplots(figsize=(3,7.75))
labelsize=12

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
max_cells = exp_info.populations[0].max_cells
n_sims = exp_info.n_sims
p_drop = exp_info.prob_drops


# make a dummy barchart
# p_drop = p_drop[2:]
x = np.arange(len(p_drop))
barchart_data = np.ones(len(p_drop))*100
rects = ax.barh(x,barchart_data,color='slategrey',facecolor='w')

tc_axes = []
drug_axes = []
pop_axes = []

exp_folders.reverse()
p_drop = np.flip(p_drop)

thresh = 1

# data_extinct = np.zeros((999,1))
for exp in exp_folders:
    p_drop_t = exp[exp.find('=')+1:]
    p_drop_t = p_drop_t.replace(',','.')
    p_drop_t = float(p_drop_t)
    num = np.argwhere(p_drop == p_drop_t)
    num = num[0,0]
    
    # generate timecourse axes
    
    width = 100
    height = rects.patches[num].get_height()
    ypos = rects.patches[num].get_y()
    xpos = 125
    popax = ax.inset_axes([xpos,ypos,width,height],transform=ax.transData)

    xpos = xpos + width + 30
    tcax = ax.inset_axes([xpos,ypos,width,height],transform=ax.transData)
    # da = tcax.twinx()
    
    sim_files = os.listdir(path=exp)
    sim_files = sorted(sim_files)
    survive_count = 0
    counts_total = None
    
    k=0
    while k < len(sim_files):
    # while k < 10:
        sim = sim_files[k]
        sim = exp + os.sep + sim
        data = results_manager.get_data(sim)
        dc = data[:,-2]
        data = data[:,0:-2]
        # data = data/np.max(data)
        data_t = data[-1,:]
        
        # check to see if any genotypes are at least 10% of the max cell count
        if any(data_t >= thresh):
            survive_count += 1
            if counts_total is None:
                counts_total = data
            else:
                counts_total += data
        
        # exp_info.populations[num].counts_log_scale = True
        d = np.sum(data,axis=1)
        popax = plotter.plot_population_count(exp_info.populations[num],
                                              d,
                                              popax,
                                              thresh=thresh,
                                              normalize=False,
                                              logscale=True)
        # popax.yaxis.tick_right()
        
        data = data/max_cells
        if k==0:
            drug_kwargs = {'alpha':1,
                            'color':'black',
                            'linewidth':1,
                            'label':'Drug Concentration ($\u03BC$M)'
                            }
            tcax,drug_ax = plotter.plot_timecourse_to_axes(exp_info.populations[num],
                                                        data,
                                                        tcax,
                                                        drug_curve=dc,
                                                        drug_curve_linestyle='-',
                                                        legend_size=3,
                                                        drug_ax_sci_notation=True,
                                                        drug_kwargs=drug_kwargs,
                                                        legend_labels=False,
                                                        grayscale=True,
                                                        color='gray', 
                                                        linewidth=1,
                                                        labelsize=10,
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
                                                        legend_size=3,
                                                        legend_labels=False,
                                                        linewidth=2,
                                                        labelsize=10,
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
                                                legend_size=3,
                                                labelsize=10)
    
    # t = np.arange(len(dc))
    # t = t*exp_info.populations[0].timestep_scale/24
    # da.plot(t,dc)
    
    tc_axes.append( tcax )
    pop_axes.append( popax )
    barchart_data[num] = survive_count   

# for a in tc_axes:
    # a.set_yscale('log',base=2)
    # a.set_ylim(10,max_cells)

# for da in drug_axes:
    # da.ticklabel_format(style='sci',axis='y',scilimits=(0,4))
    # da.set_yticks(da.get_yticks())
    # yt = da.get_yticks
    # yt = yt/10**3

# drug_axes[1].set_ylabel('Drug Concentration (uM)', color='gray',fontsize=labelsize)
# tc_axes[3].set_ylabel('Proportion of \nmax cell count',fontsize=labelsize)
# tc_axes[0].set_xlabel('Days',fontsize=labelsize)
pop_axes[0].set_xlabel('Days',fontsize=labelsize)

rects = ax.barh(x,barchart_data,color='slategrey')
ax.set_yticks(x)
ax.set_yticklabels(p_drop)
# ax.yaxis.set_major_formatter(ScalarFormatter())
# ax.ticklabel_format(style='sci',axis='y')
ax.set_xlabel('% Resistant',fontsize=12)
ax.set_ylabel('$p_{forget}$',fontsize=12)
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
                   bbox_to_anchor=(-0.25,-0.9),ncol=4)
drug_axes[0].legend(frameon=False,fontsize=11,loc='lower left',
                    bbox_to_anchor=(-0.8,-1.4),ncol=1)

# ax.annotate('Proportion of carrying capacity',rotation=90,
#             fontsize=labelsize,xy=(97,0.65),
#             ha='center') # xy in data points

# ax.annotate('Drug Concentration ($\u03BC$M)',rotation=90,
#             fontsize=labelsize,xy=(242,0.75),
#             ha='center',annotation_clip=False) # xy in data points

# ax.annotate('Cell count',rotation=90,
#             fontsize=labelsize,xy=(380,1.15),
#             ha='center',annotation_clip=False)

# ax.annotate('Cell count',rotation=90,
#             fontsize=labelsize,xy=(97,1.15),
#             ha='center',annotation_clip=False)

# for da in drug_axes:
#     da.set_yticks([0,1e3,2e3])
#     da.ticklabel_format(style='sci',axis='y',scilimits=(0,3))

for pa in pop_axes:
    xl = pa.get_xlim()
    xl = [0,xl[1]]
    pa.set_xlim(xl)
    


# # xl = tc_axes[0].get_xlim()
# # xt = tc_axes[0].get_xticks()
# for pa in pop_axes:
#     pa.set_xlim(xl)
#     pa.set_xticks(xt)
#     pa.set_ylim([10**0,10**12])

handles, labels = pop_axes[-1].get_legend_handles_labels()
r_index = labels.index('resistant')
e_index = labels.index('extinct')
pop_axes[-1].legend([handles[r_index],handles[e_index]],
                    ['resistant','extinct'],frameon=False)

results_manager.save_fig(fig,'adherance_v_survival.pdf',bbox_inches='tight')