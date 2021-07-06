from fears.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

data_folder = 'results_07022021_0000'
exp_info_file = 'experiment_info_07022021_0000.p'
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

f,ax[0,1] = plotter.plot_fitness_curves(exp_info.p_seascape,
                                    ax=ax[0,1], 
                                    show_legend=False,
                                    show_axes_labels=False,
                                    labelsize=labelsize,
                                    linewidth=linewidth)

# timecourse axes
landscape_exp = exp_folder[0]
data = results_manager.get_data(landscape_exp)
counts = data[:,0:4]

ax[1,0],d = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,0],
                                    labelsize=labelsize,
                                    linewidth=linewidth)

seascape_exp = exp_folder[1]
data = results_manager.get_data(seascape_exp)
counts = data[:,0:4]
ax[1,1],d = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,1],
                                    labelsize=labelsize,
                                    linewidth=linewidth)

# plot the location of drug conc transition in the timecourse axes

# ydata = np.logspace(1,11)
# xdata = np.ones(len(ydata))*exp_info.transition_times[0]
# ax[1,0].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)
# ax[1,1].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)

# ydata = np.logspace(1,11)
# xdata = np.ones(len(ydata))*exp_info.transition_times[1]
# ax[1,0].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)
# ax[1,1].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)

# landscape axes

null_ax = ax[0,0]
conc = [exp_info.first_dose,exp_info.second_dose,exp_info.third_dose]
cmap = 'Blues'
edgecolor='black'
textcolor='gray'

for c in conc:
    plotter.add_landscape_to_fitness_curve(c,null_ax,exp_info.p_landscape,
                                           textcolor=textcolor,
                                           cmap=cmap,
                                           edgecolor=edgecolor,
                                           linewidths=0.5,
                                           textsize=9)
    
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
    a.set_ylabel('Growth rate')
    a.set_xlabel('Drug concentration ($\u03BC$M)')
    
for a in ax[1,:]:
    a.set_ylabel('Cell count')
    a.set_xlabel('Days')
    
ax[1,1].legend(frameon=False,fontsize=10,loc=(1,0.2))
    
for row in range(2):
    for col in range(2):
        a = ax[row,col]
        pos = [left[row,col],bottom[row,col],w,h]
        a.set_position(pos)
        
results_manager.save_fig(fig,'test.png',bbox_inches=None)