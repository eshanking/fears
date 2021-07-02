from fears.utils import plotter, results_manager
import matplotlib.pyplot as plt
import numpy as np

data_folder = 'results_07012021_0002'
exp_info_file = 'experiment_info_07012021_0002.p'
exp_folder,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
# fitness axes
fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(7,5))
linewidth = 2

f,ax[0,0] = plotter.plot_fitness_curves(exp_info.p_landscape,
                                    ax=ax[0,0],
                                    show_legend=False,
                                    show_axes_labels=False,
                                    labelsize=10,
                                    linewidth=linewidth)

f,ax[0,1] = plotter.plot_fitness_curves(exp_info.p_seascape,
                                    ax=ax[0,1],
                                    show_legend=False,
                                    show_axes_labels=False,
                                    labelsize=10,
                                    linewidth=linewidth)

# timecourse axes
landscape_exp = exp_folder[0]
data = results_manager.get_data(landscape_exp)
counts = data[:,0:4]

ax[1,0],d = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,0],
                                    labelsize=10,
                                    linewidth=linewidth)

seascape_exp = exp_folder[1]
data = results_manager.get_data(seascape_exp)
counts = data[:,0:4]
ax[1,1],d = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1,1],
                                    labelsize=10,
                                    linewidth=linewidth)

# plot the location of drug conc transition in the timecourse axes

ydata = np.logspace(1,11)
xdata = np.ones(len(ydata))*exp_info.transition_times[0]
ax[1,0].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)
ax[1,1].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)

ydata = np.logspace(1,11)
xdata = np.ones(len(ydata))*exp_info.transition_times[1]
ax[1,0].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)
ax[1,1].plot(xdata,ydata,color='gray',alpha=0.5,linewidth=5)

# landscape axes

def get_pos_in_log_space(center,width):
    """
    Returns x or y coordinates to pass to inset_axes when the axis is a log 
    scale in order to make the inset axes uniform sizes.

    Parameters
    ----------
    center : float
        Center point of the axes in data points (not log-ed).
    width : float
        Visual width of the axes (in log units).

    Returns
    -------
    x : list of floats
        x[0] = left point, x[1] = right point.
        log10(x[1]-x[0]) = width

    """
    
    center = np.log10(center)
    x0 = center - width/2
    x1 = center + width/2
    
    x = [10**x0,10**x1]
    
    return x

null_ax = ax[0,0]
conc = [exp_info.first_dose,exp_info.second_dose,exp_info.third_dose]

for c in conc:
    x = get_pos_in_log_space(c,3)
    l1 = null_ax.inset_axes([x[0],1,x[1]-x[0],0.5],transform=null_ax.transData)
    l1 = plotter.plot_landscape(exp_info.p_landscape,c,ax=l1,node_size=200,
                            colorbar=False,
                            textcolor='gray',
                            square=True)
    
    yl = null_ax.get_ylim()
    ydata = np.arange(yl[0],yl[1],0.1)
    xdata = np.ones(len(ydata))*c
    null_ax.plot(xdata,ydata,'--',color='black',alpha=0.5)
    
    # xl = null_ax.get_xlim()
    # null_ax.set_xlim([xl[0]/10,xl[1]*10])
    
sea_ax = ax[0,1]
for c in conc:
    x = get_pos_in_log_space(c,3)
    l1 = sea_ax.inset_axes([x[0],1,x[1]-x[0],0.5],transform=sea_ax.transData)
    l1 = plotter.plot_landscape(exp_info.p_seascape,c,ax=l1,node_size=200,
                            colorbar=False,
                            textcolor='gray',
                            square=True)
    
    yl = sea_ax.get_ylim()
    ydata = np.arange(yl[0],yl[1],0.1)
    xdata = np.ones(len(ydata))*c
    sea_ax.plot(xdata,ydata,'--',color='black',alpha=0.5)
    
# cb.show()
    

# pos00 = ax[0,0].get_position()
# pos01 = ax[0,1].get_position()
# pos10 = ax[1,0].get_position()
# pos11 = ax[1,1].get_position()

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