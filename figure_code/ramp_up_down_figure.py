import matplotlib.pyplot as plt
from fears.utils import results_manager, plotter

data_folder = 'results_06262021_0001'
info_file = 'experiment_info_06262021_0001.p'

exp_folders,exp_info = results_manager.get_experiment_results(data_folder,
                                                              info_file)
p = exp_info.p_landscape
plt.box(False)
fig,ax = plt.subplots(2,figsize=(5,4))

landscape_exp = exp_folders[0]
data = results_manager.get_data(landscape_exp)
data = data[:,0:16]
p.counts_log_scale = True
ax[0],da = plotter.gen_timecourse_axes(p,data,ax[0])

seascape_exp = exp_folders[1]
data = results_manager.get_data(seascape_exp)
data = data[:,0:16]
p.counts_log_scale = True
ax[1],da = plotter.gen_timecourse_axes(p,data,ax[1])

for a in ax:
    # a.set_ylim(0,1.1*10**6)
    a.set_ylabel('Cells',fontsize=15)
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

# ax[0].set_title('Landscape',fontsize=15)
# ax[1].set_title('Seascape',fontsize=15)
pos_0 = ax[0].get_position()
pos_1 = ax[1].get_position()
pos_1.y1 = pos_1.y1-0.05
pos_1.y0 = pos_1.y0-0.05
ax[1].set_position(pos_1)
ax[1].set_xlabel('Days',fontsize=15)
ax[1].legend(frameon=False,fontsize=11,loc=(0,-1),bbox_to_anchor=(0,-0.75),ncol=4)

# plt.gca().set_axis_off()
# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
#             hspace = 0, wspace = 0)

results_manager.save_fig(fig,'ramp_up_down.png')




