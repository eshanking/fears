from fears.utils import plotter, results_manager
import matplotlib.pyplot as plt

data_folder = 'results_06302021_0001'
exp_info_file = 'experiment_info_06302021_0001.p'
exp_folder,exp_info = results_manager.get_experiment_results(data_folder,
                                                             exp_info_file)
# fitness axes
fig,ax = plt.subplots(3,figsize=(3,6))

d,ax[0] = plotter.plot_fitness_curves(exp_info.p_landscape,
                                    ax=ax[0],
                                    show_legend=False,
                                    show_axes_labels=False)

# timecourse axes
landscape_exp = exp_folder[0]
data = results_manager.get_data(landscape_exp)
counts = data[:,0:4]

ax[1] = plotter.plot_timecourse_to_axes(exp_info.p_landscape,
                                    counts,
                                    ax[1])

# landscape axes

ep = plotter.plot_landscape(exp_info.p_landscape,ax=ax[2],colorbar=False,resize_param=0.2,square=False)

ax[0].legend()