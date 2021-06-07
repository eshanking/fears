# https://matplotlib.org/api/_as_gen/matplotlib.gridspec.GridSpec.html#matplotlib.gridspec.GridSpec
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import seaborn as sns
from cycler import cycler
# from population_class import int_to_binary

# method for computing binary allele label from allele number
def int_to_binary(num, pad=4):
    return bin(num)[2:].zfill(pad)

data_folder = 'results_09172020_0001'
max_cells = 0.6*10**6 # carrying capacity for determining if a population survived or not

###############################################################################
# generate figure and axes

results_dir = os.getcwd() + os.sep + data_folder # folder containing all of the results
experiment_folders = os.listdir(path=results_dir) # each of these folders corresponds to a different k_abs

n_params = len(experiment_folders)

# Full page figure with room for caption
fig = plt.figure(figsize=(6.25,7.75))

# Generate a 3x10 grid for axes
gs = GridSpec(n_params,2,figure=fig,wspace=0.1,width_ratios=[1,1])

bar_ax = fig.add_subplot(gs[:,0]) # axis for the bar chart

drug_ax = [] # axes for the drug concentration timecourse
all_ax = [] # axes where I will plot all time traces
for i in range(n_params):
    ax = fig.add_subplot(gs[i,1])
    ax.yaxis.tick_right()
    all_ax.append(ax)
    drug_ax.append(ax.twinx()) # make the drug concentration timecourse share axes with the population timecourse

###############################################################################
# data analysis

barchart_labels = [] # will eventually contain the k_abs values
barchart_data = [] # will eventually contain the survival data for the barchart

experiment_folders = sorted(experiment_folders)

row = 0 # counts the current row of the figure
for exp in experiment_folders:
    k_abs = exp[exp.find('=')+1:] # determine the k_abs from the folder name
    
    if len(k_abs)>7: # truncate the absorption parameter for labeling purposes
        k_abs = k_abs[0:7]
        
    exp_path = results_dir + os.sep + exp
    sim_files = os.listdir(path=exp_path) # each file corresponds to a single simulation
    
    n_survive = 0
    n_sim = len(sim_files)
    
    k=0
    for sim in sim_files:

        sim_path = exp_path + os.sep + sim
        
        data_df = pd.read_csv(sim_path)
        data = data_df.to_numpy()
        
        drug_conc = data[:,-1] # drug concentration is the 17th column of the data
        counts = data[:,0:-1]
        
        if any(counts[-1,:]>0.1*max_cells):
            n_survive+=1
            if k==0:
                counts_t = counts
            else:
                counts_t += counts
            k+=1
        # else:
        #     counts_t = np.zeros(counts.shape)
        
        counts = counts/max_cells # normalize cell counts to the carrying capacity     
        # plot each simulation trace in grey
        for allele in range(counts.shape[1]):
            all_ax[row].plot(counts[:,allele],linewidth=1.0,color='grey',alpha=0.5)
            
    counts_t = counts_t/n_survive
    counts_total = np.sum(counts,axis=0)
    
    # counts_total = counts_total/np.max(counts_total)    
    sorted_index = counts_total.argsort()
    sorted_index_big = sorted_index[8:] # sort the top 8 alleles with the largest cell count
    
    colors = sns.color_palette('bright')
    colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
    
    # shuffle colors
    
    colors[[14,15]] = colors[[15,14]]
    
    # cycle through dashed and solid lines
    cc = (cycler(color=colors) + 
          cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                            '--','--','--','--','--','--','--']))
    
    all_ax[row].set_prop_cycle(cc)
    
    counts_t = counts_t/max_cells
    for allele in range(counts.shape[1]):
        if allele in sorted_index_big:
            # only label the trace for the legend if the cell count is large enough (don't want to bother putting tiny populations in the legend)
            all_ax[row].plot(counts_t[:,allele],linewidth=2.0,label=str(int_to_binary(allele)))
#            print(str(allele))
        else:
            all_ax[row].plot(counts_t[:,allele],linewidth=2.0,label=None)
        # if k==0:
        #     counts_avg
        
    # plot the drug concentration
    drug_ax[row].plot(drug_conc,color='black',label='Drug Concentration ($\u03BC$M)')
    
                
    barchart_data.append(100*n_survive/n_sim)
    barchart_labels.append(k_abs)
    
    row+=1

##############################################################################
# plotting

absorption_rate = [float(k_abs) for k_abs in barchart_labels]
# x = np.arange(10,0,-1)-.5
# barchart_data = [10,9,8,7,6,5,4,3,2,1]
# x = [10,9,8,7,6,5,4,3,2,1]
# x = 10*(np.arange(n_params,0,-1))/n_params
x = np.arange(n_params,0,-1)

# compute error bars
p = np.array(barchart_data)/100
n = 100 # number of simulations
q = 1-p

sd = 100*(p*q/n)**0.5 # standard deviation of the estimator of the parameter of a bernoulli distribution

rects = bar_ax.barh(x, barchart_data,color='slategrey')
errorbar_pos = x + rects[0].get_height()/2
# bar_ax.errorbar(x=barchart_data, y=errorbar_pos, xerr=sd,linewidth=0,elinewidth=2,capsize=5,color='tomato')
bar_ax.errorbar(x=barchart_data, y=errorbar_pos, xerr=sd,linewidth=0,elinewidth=2,capsize=5,color='black')

# position the barchart rectangles exactly where I want them
i=0
for rect in rects:
    # rect.set_y(rect.get_y()-.1)
    # print(str(rect.get_width()))
    rect.set_y(x[i])
    i+=1

labels = ['0','25','50','75','100']
bar_ax.set_ylim(min(x)-0.05,max(x)+0.85)
bar_ax.set_xticks([0,25,50,75,100])
bar_ax.set_xticklabels(labels)

bar_ax.set_yticks(x+0.4)
bar_ax.set_yticklabels(barchart_labels)
# bar_ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

bar_ax.set_xlabel('% Resistant',fontsize=12)
bar_ax.set_ylabel('$k_{abs}$',fontsize=15)
bar_ax.spines["right"].set_visible(False)
bar_ax.spines["top"].set_visible(False)
    
for i in range(len(all_ax)-1):
    ax=all_ax[i]
    ax.tick_params(labelbottom=False)
    ax.yaxis.tick_right()
    # ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    
# all_ax[-1].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
all_ax[-1].set_xlabel('Time Step',fontsize=12)
all_ax[-1].yaxis.tick_right()
all_ax[-1].legend(frameon=False,fontsize=11,loc='lower left',bbox_to_anchor=(-1.1,-1.05),ncol=4)
drug_ax[-1].legend(frameon=False,fontsize=11,loc='lower left',bbox_to_anchor=(-1.1,-1.3),ncol=1)
# all_ax[-1].legend(loc=(1.25,-.12),frameon=False,fontsize=15)

for ax in drug_ax:
    ax.yaxis.tick_left()
    ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    ax.set_ylim(0,110)
    
for ax in all_ax:
    ax.set_ylim(-0.05,1.05)


all_ax[int(np.floor(len(all_ax)/2))].set_ylabel('Proportion of Max Cell Count',fontsize=12)
all_ax[int(np.floor(len(all_ax)/2))].yaxis.set_label_position("right")

drug_ax[int(np.floor(len(all_ax)/2))].set_ylabel('Drug Concentration ($\u03BC$M)',fontsize=12)
drug_ax[int(np.floor(len(all_ax)/2))].yaxis.set_label_position("left")


plt.savefig('pk_rate_barchart.svg',bbox_inches="tight")