import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import numpy as np
import os
from fears.src import utils

def plot_timecourse(pop,counts_t=None,title_t=None):
    
    if (pop.counts == 0).all() and counts_t is None:
        print('No data to plot!')
        return
    elif counts_t is None:
        counts = pop.counts
    else:
        counts = counts_t # an input other than pop overrides pop
    if title_t is not None:
        title = title_t
    else:
        title = pop.fig_title    
        
    left = 0.1
    width = 0.8
    
    if pop.plot_entropy == True:
        fig,(ax1,ax3) = plt.subplots(2,1,figsize=(6,4),sharex=True) 
        ax3.set_position([left, 0.2, width, 0.2]) # ax3 is entropy
    else:
        fig,ax1 = plt.subplots(1,1,figsize=(6,4),sharex=True)
    
    ax1.set_position([left, 0.5, width, 0.6]) # ax1 is the timecourse
            
    counts_total = np.sum(counts,axis=0)
    
    sorted_index = counts_total.argsort()
    sorted_index_big = sorted_index[-8:]
    
    colors = sns.color_palette('bright')
    colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
    
    # shuffle colors
    colors[[14,15]] = colors[[15,14]]
    
    cc = (cycler(color=colors) + 
          cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                            '--','--','--','--','--','--','--']))
    
    ax1.set_prop_cycle(cc)

    color = [0.5,0.5,0.5]
    
    if pop.fitness_data == 'generate':
        ax2 = ax1.twinx() # ax2 is the drug timecourse
        ax2.set_position([left, 0.5, width, 0.6])
        ax2.set_ylabel('Drug Concentration (uM)', color=color,fontsize=20) # we already handled the x-label with ax1
        
        if pop.drug_log_scale:
            if all(pop.drug_curve>0):
                drug_curve = np.log10(pop.drug_curve)
            yticks = np.log10([10**-4,10**-3,10**-2,10**-1,10**0,10**1,10**2,10**3])    
            ax2.set_yticks(yticks)
            ax2.set_yticklabels(['0','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                             '$10^1$','$10^2$','$10^3$'])
            ax2.set_ylim(-4,3)
        else:
            drug_curve = pop.drug_curve
            ax2.set_ylim(0,1.1*max(drug_curve))
    
        ax2.plot(drug_curve, color=color, linewidth=2.0)
        ax2.tick_params(axis='y', labelcolor=color)
            
        ax2.legend(['Drug Conc.'],loc=(1.25,0.93),frameon=False,fontsize=15)
        
        ax2.tick_params(labelsize=15)
        ax2.set_title(title,fontsize=20)
    
    # if pop.normalize:
    #     counts = counts/np.max(counts)
        
    for allele in range(counts.shape[1]):
        if allele in sorted_index_big:
            ax1.plot(counts[:,allele],linewidth=3.0,label=str(pop.int_to_binary(allele)))
        else:
            ax1.plot(counts[:,allele],linewidth=3.0,label=None)
            
    ax1.legend(loc=(1.25,-.12),frameon=False,fontsize=15)
        
    ax1.set_xlim(0,pop.x_lim)
    ax1.set_facecolor(color='w')
    ax1.grid(False)

    ax1.set_ylabel('Cells',fontsize=20)
    ax1.tick_params(labelsize=15)
    
    if pop.plot_entropy == True:
        e = pop.entropy(counts)
        
        ax3.plot(e,color='black')
        ax3.set_xlabel('Time',fontsize=20)
        ax3.set_ylabel('Entropy',fontsize=20)
        if pop.entropy_lim is not None:
            ax3.set_ylim(0,pop.entropy_lim)
        ax3.tick_params(labelsize=15)
    
    if pop.y_lim is not None:
        y_lim = pop.y_lim
    else:
        y_lim = np.max(counts) + 0.05*np.max(counts)
    
    if pop.counts_log_scale:
        ax1.set_yscale('log')
        ax1.set_ylim(1,5*10**5)
    else:
        ax1.set_ylim(0,y_lim)
    
    xlabels = ax1.get_xticks()
    xlabels = xlabels*pop.timestep_scale
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    ax1.set_xticklabels(xlabels)
    ax1.set_xlabel('Days',fontsize=20)

    plt.show()
    return fig

def plot_fitness_curves(pop,fig_title='',plot_r0 = False,save=False):
    
    drugless_rates = pop.drugless_rates
    ic50 = pop.ic50
    
    fig, ax = plt.subplots(figsize = (10,6))
    
    powers = np.linspace(-3,5,20)
    conc = np.power(10*np.ones(powers.shape[0]),powers)
    fit = np.zeros(conc.shape[0])
    
    colors = sns.color_palette('bright')
    colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
    colors[[14,15]] = colors[[15,14]]
    
    cc = (cycler(color=colors) + 
           cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                            '--','--','--','--','--','--','--']))
    ax.set_prop_cycle(cc) 
    
    for allele in range(16):
        for j in range(conc.shape[0]):
            fit[j] = pop.gen_fitness(allele,conc[j],drugless_rates,ic50)
            
        if plot_r0:
            fit = fit - pop.death_rate
            ylabel = '$R_{0}$'
            thresh = np.ones(powers.shape)
            ax.plot(powers,thresh,linestyle='dashdot',color='black',linewidth=3)
        else:
            ylabel = 'Growth Rate'
            
        ax.plot(powers,fit,linewidth=3,label=str(pop.int_to_binary(allele)))
    
    ax.legend(fontsize=15,frameon=False,loc=(1,-.10))
    ax.set_xticks([-3,-2,-1,0,1,2,3,4,5])
    ax.set_xticklabels(['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$',
                         '$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
    
    plt.title(fig_title,fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    plt.xlabel('Drug concentration ($\mathrm{\mu}$M)',fontsize=20)
    plt.ylabel(ylabel,fontsize=20)
    ax.set_frame_on(False)
    
    if save:
        plt.savefig('fitness_curve.pdf',bbox_inches="tight")
    
    return fig

def plot_msw(pop,fitness_curves,conc,genotypes,save=False):
    """
    plot_msw: method for plotting mutant selection window figures.

    Parameters
    ----------
    pop : population_class object
        
    fitness_curves : numpy array
        Columns 1-N represents a genotype that is a neighbor of column 0 
        (ancestor). Rows represent drug concentration.
    conc : numpy array
        Drug concentration used to calculate fitness_curves
    genotypes : list of ints
        Genotypes that were used to calculate the fitness_curves.
    save : bool
    
    Returns
    -------
    fig : figure object
        MSW figures

    """
    n_genotype = fitness_curves.shape[1]
    rows = int((n_genotype-1)/2)
    fig, ax = plt.subplots(rows,2)
    g = 1
    wt_fitness_curve = fitness_curves[:,0]
    for r in range(rows):
        for col in range(2):
           
            ax[r,col].plot(conc,wt_fitness_curve,label='wt',linewidth=3)
            
            cur_fitness_curve = fitness_curves[:,g]
            gt = genotypes[g]
            bitstring = pop.int_to_binary(gt)    
            ax[r,col].plot(conc,cur_fitness_curve,label=bitstring,linewidth=3)
            
            msw_left_assigned = False
            msw_right_assigned = False
            if wt_fitness_curve[0] > cur_fitness_curve[0] \
                and any(cur_fitness_curve>wt_fitness_curve):
                for c in range(len(conc)):
                    if wt_fitness_curve[c] < cur_fitness_curve[c] \
                        and msw_left_assigned is False:
                        msw_left = conc[c]
                        msw_left_assigned = True
                    if (cur_fitness_curve[c] < 1 
                        and msw_right_assigned is False):
                        msw_right = conc[c]
                        msw_right_assigned = True
                if msw_left < msw_right:
                    ax[r,col].axvspan(msw_left, msw_right, 
                                      facecolor='#2ca02c',alpha=0.5,
                                      label='MSW')
            
            ax[r,col].set_xscale('log')
            ax[r,col].legend(fontsize=10,frameon=False)

            g+=1
            
    for r in range(rows):
        ax[r,0].set_ylabel('$R_{0}$',fontsize=10)
    for c in range(2):
        ax[rows-1,c].set_xlabel('Drug concentration ($\mathrm{\mu}$M)',
                              fontsize=10)
    if save:
        r = utils.get_project_root()
        savename = str(r) + os.sep + 'figures' + os.sep + 'msw.pdf'
        plt.savefig(savename,bbox_inches="tight")
    
    return fig












