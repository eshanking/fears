import matplotlib.pyplot as plt
import numpy as np
from fears.utils import results_manager, plotter

shift1 = 0.02
shift2 = 0.065
data_folder = 'results_06282021_0001'
info_file = 'experiment_info_06282021_0001.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                              info_file)
p = exp_info.p_landscape
plt.box(False)
fig,ax = plt.subplots(3,figsize=(3.5,6))

landscape_exp = exp_folders[0]
data = results_manager.get_data(landscape_exp)
data = data[:,0:16]
p.counts_log_scale = True

ax[0],da = plotter.gen_timecourse_axes(p,data,ax[0],labelsize=10,linewidth=2)

seascape_exp = exp_folders[1]
data = results_manager.get_data(seascape_exp)
data = data[:,0:16]
p.counts_log_scale = True
ax[1],da = plotter.gen_timecourse_axes(p,data,ax[1],labelsize=10,linewidth=2)

trans1 = exp_info.transition_times[0] - exp_info.ramp/2
trans2 = exp_info.transition_times[1] + exp_info.ramp/2

for i in range(2):
    # a.set_ylim(0,1.1*10**6)
    a = ax[i]
    c = [0.8,0.8,0.8]
    a.axvspan(0,trans1,facecolor=c)
    c = [0.1,0.1,0.1]
    a.axvspan(trans1,trans2,facecolor=c)
    c = [0.5,0.5,0.5]
    end_time = a.get_xlim()[1]
    a.axvspan(trans2,end_time,facecolor=c)
    a.set_ylabel('Cells',fontsize=12)
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.spines["left"].set_visible(False)
    a.spines["bottom"].set_visible(False)

# ax[0].set_title('Landscape',fontsize=15)
# ax[1].set_title('Seascape',fontsize=15)
x = np.log10(exp_info.first_dose*np.ones(20))
gr = np.arange(20)/19
ax[2].plot(x,gr,linewidth=4,color=[0.8,0.8,0.8],label='First dose')

x = np.log10(exp_info.second_dose*np.ones(20))
gr = np.arange(20)/19
ax[2].plot(x,gr,linewidth=4,color=[0.1,0.1,0.1],label='Second dose')

x = np.log10(exp_info.third_dose*np.ones(20))
gr = np.arange(20)/19
ax[2].plot(x,gr,linewidth=4,color=[0.5,0.5,0.5],label='Third dose')

ax[2] = plotter.plot_fitness_curves(exp_info.p_seascape,ax=ax[2],linewidth=2,labelsize=10)

ax[2].set_xlabel('Drug concentration ($\mathrm{\mu}$M)',fontsize=12)
ax[2].set_ylabel('Growth Rate',fontsize=12)

handles, labels = ax[2].get_legend_handles_labels()
# ax[2].legend(handles[-3:],labels[-3:],ncol=3,loc=(-0.3,-0.65),frameon=False)

pos_0 = ax[0].get_position()
pos_1 = ax[1].get_position()
pos_2 = ax[2].get_position()

pos_1.y1 = pos_1.y1-shift1
pos_1.y0 = pos_1.y0-shift1

pos_2.y1 = pos_2.y1-shift2
pos_2.y0 = pos_2.y0-shift2

ax[1].set_position(pos_1)
ax[2].set_position(pos_2)

ax[1].set_xlabel('Days',fontsize=12)

# ax[1].legend(frameon=False,fontsize=11,loc=(0,-1),bbox_to_anchor=(0,-0.75),ncol=4)

# plt.gca().set_axis_off()
# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
#             hspace = 0, wspace = 0)

results_manager.save_fig(fig,'ramp_up_down.png')