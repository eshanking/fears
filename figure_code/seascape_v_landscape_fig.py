from fears.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

data_folder = 'results_08262021_0000'
exp_info_file = 'experiment_info_08262021_0000.p'
exp_folder,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
# fitness axes
fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(7,5))
linewidth = 2

labelsize=8

f,ax[0,0] = plotter.plot_fitness_curves(exp_info.p_landscape,
                                    ax=ax[0,0],
                                    show_legend=False,
                                    show_axes_labels=False,
                                    labelsize=labelsize,
                                    linewidth=linewidth)

ax[0,0].set_xticks([10**-3,10**-1,10**1,10**3,10**5])
ax[0,0].xaxis.tick_top()

f,ax[0,1] = plotter.plot_fitness_curves(exp_info.p_seascape,
                                    ax=ax[0,1], 
                                    show_legend=False,
                                    show_axes_labels=False,
                                    labelsize=labelsize,
                                    linewidth=linewidth)
ax[0,1].set_xticks([10**-3,10**-1,10**1,10**3,10**5])
ax[0,1].xaxis.tick_top()

# timecourse axes
landscape_exp = exp_folder[0]
data = results_manager.get_data(landscape_exp)
counts = data[:,0:4]
dc = exp_info.p_landscape.drug_curve
drug_kwargs = {'color':'black',
               'alpha':0.5}

ax[1,0],drug_ax = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,0],
                                    labelsize=labelsize,
                                    linewidth=linewidth,
                                    drug_curve=dc,
                                    drug_curve_linestyle='--',
                                    drug_curve_label=False,
                                    drug_kwargs=drug_kwargs)

drug_ax.set_ylim([10**-5,10**7])
drug_ax.set_yticks([10**-3,10**1,10**5])

seascape_exp = exp_folder[1]
data = results_manager.get_data(seascape_exp)
counts = data[:,0:4]
ax[1,1],drug_ax = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,1],
                                    labelsize=labelsize,
                                    linewidth=linewidth,
                                    drug_curve_linestyle='--',
                                    drug_curve=dc,
                                    drug_kwargs=drug_kwargs)

drug_ax.set_ylim([10**-5,10**7])
drug_ax.set_yticks([10**-3,10**1,10**5])
# landscape axes

null_ax = ax[0,0]
conc = [exp_info.first_dose,exp_info.second_dose,exp_info.third_dose]
cmap = 'Blues'
edgecolor='black'
textcolor='goldenrod'

for c in conc:
    plotter.add_landscape_to_fitness_curve(c,null_ax,exp_info.p_landscape,
                                           textcolor=textcolor,
                                           cmap=cmap,
                                           edgecolor=edgecolor,
                                           linewidths=0.5,
                                           textsize=9,
                                           position='bottom',
                                           pad=0.1)
    
sea_ax = ax[0,1]

for i in range(len(conc)-1):
    c = conc[i]
    plotter.add_landscape_to_fitness_curve(c,sea_ax,exp_info.p_seascape,
                                           textcolor=textcolor,
                                           cmap=cmap,
                                           edgecolor=edgecolor,
                                           linewidths=0.5,
                                           textsize=9)

c = conc[-1]
# cbax = fig.add_subplot()
l1 = plotter.add_landscape_to_fitness_curve(c,sea_ax,exp_info.p_seascape,
                                            textcolor=textcolor,
                                            cmap=cmap,
                                            edgecolor=edgecolor,
                                            linewidths=0.5,
                                            textsize=9,
                                            colorbar=True)

# reposition axes
w = 0.3
h = 0.25

wspace = (1-2*w)/3
hspace = (1-2*h)/2.7

bottom = np.array([[1-hspace-h,1-hspace-h],[hspace,hspace]])
left = np.array([[wspace,1-wspace-w],[wspace,1-wspace-w]])

for a in ax[0,:]:
    # a.set_ylabel('Growth rate',fontsize=labelsize)
    a.set_xlabel('Drug concentration ($\u03BC$M)',fontsize=labelsize)
    
for a in ax[1,:]:
    # a.set_ylabel('Cell count',labelpad=0,fontsize=labelsize)
    a.set_xlabel('Days',fontsize=labelsize)
    
ax[1,0].set_ylabel('Cell count',labelpad=0,fontsize=labelsize)
ax[0,0].set_ylabel('Growth rate',fontsize=labelsize)
    
ax[1,1].legend(frameon=False,fontsize=7,
               bbox_to_anchor=(-0.7, -0.40, 1., .102), loc='lower left',
               ncol=4, mode="expand", borderaxespad=0.)
    
for row in range(2):
    for col in range(2):
        a = ax[row,col]
        pos = [left[row,col],bottom[row,col],w,h]
        a.set_position(pos)
        
ax[0,0].annotate('a.', xy=(-0.15,1.05),  xycoords='axes fraction')
ax[1,0].annotate('c.', xy=(-0.15,1.05),  xycoords='axes fraction')
ax[0,1].annotate('b.', xy=(-0.15,1.05),  xycoords='axes fraction')
ax[1,1].annotate('d.', xy=(-0.15,1.05),  xycoords='axes fraction')
        
results_manager.save_fig(fig,'seascape_v_landscape.pdf',bbox_inches='tight')