# import sys
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import numpy as np
# import os
import math
import scipy.stats
from fears.utils import fitness
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import networkx as nx
from labellines import labelLine


def gen_color_cycler(style=None,palette='bright',n_colors=16):
    """Generates a custom matplotlib color cycler

    Args:
        style (str, optional): Determines linestyle. If solid, all lines are solid.
        Else, half are dashed and half are solid. Defaults to None.
        palette (str, optional): Matplotlib palette. Defaults to 'bright'.
        n_colors (int, optional): Number of colors to cycle through. Defaults to 16.

    Returns:
        cycler: Matplotlib cycle object
    """
    
    if style is None:
        colors = sns.color_palette(palette)
        colors = np.concatenate((colors[0:9],colors[0:7]),axis=0)
        # colors[[14,15]] = colors[[15,14]]
        
        colors[[7,8]] = colors[[8,7]]
        
        cc = (cycler(color=colors) + 
            cycler(linestyle=['-', '-','-','-','-','-','-','-','-',
                                '--','--','--','--','--','--','--']))
    elif style == 'solid':
        colors = sns.color_palette(palette,n_colors+1)
        # reshuffle
        colors_t = colors
        colors[3] = colors_t[4]
        colors = colors[:4]
        cc = cycler(color=colors)
    return cc

def plot_timecourse(pop,counts_t=None,title_t=None,**kwargs):
    """Plot timecourse of evolution experiment

    Args:
        pop (population): Population class object
        counts_t (array-like, optional): Experiment counts. Columns are different 
        genotypes and rows are timesteps. If None, gets data from population object.
        Defaults to None.
        title_t (str, optional): Figure title. Defaults to None.

    Returns:
        figure: Matplotlib figure object
    """
    
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
    
    cc = gen_color_cycler()
    
    ax1.set_prop_cycle(cc)

    color = [0.5,0.5,0.5]
    
    if pop.plot_drug_curve:
        ax2 = ax1.twinx() # ax2 is the drug timecourse
        ax2.set_position([left, 0.5, width, 0.6])
        ax2.set_ylabel('Drug Concentration \n($\u03BC$M)', color=color,fontsize=14) # we already handled the x-label with ax1
        
        drug_curve = pop.drug_curve
        
        ax2.plot(drug_curve, color='black', linewidth=2.0)
        ax2.tick_params(axis='y', labelcolor=color)
        
        if pop.drug_log_scale:
            ax2.set_yscale('log')
            if min(drug_curve) <= 0:
                axmin = 10**-3
            else:
                axmin = min(drug_curve)
            ax2.set_ylim(axmin,2*max(drug_curve))
            ax2.legend(['Drug Conc.'],loc=(1.3,0.93),frameon=False,fontsize=12)
            
        else:
            ax2.set_ylim(0,1.1*max(drug_curve))
            ax2.legend(['Drug Conc.'],loc=(1.25,0.93),frameon=False,fontsize=12)

            
        ax2.tick_params(labelsize=12)
        ax2.axes.ticklabel_format(axis='y',style='sci',scilimits=(1,10))
        ax2.set_title(title,fontsize=14)
        
    for allele in range(counts.shape[1]):
        if allele in sorted_index_big:
            ax1.plot(counts[:,allele],linewidth=3.0,label=str(pop.int_to_binary(allele)))
        else:
            ax1.plot(counts[:,allele],linewidth=3.0,label=None)

    if pop.plot_pop_size:
        c = np.sum(counts,axis=1)
        ax1.plot(c,'--',color='black',linewidth=3)

        ext_time = np.argwhere(c<=1)
        if len(ext_time) == 0:
            ext_time = 'None'
        else:
            ext_time = str(ext_time[0][0])

        title = 'Min population size = ' + str(np.min(c)) +\
            ' Ext. time = ' + ext_time

        ax1.set_title(title)

    ax1.legend(loc=(1.25,-.12),frameon=False,fontsize=12)
        
    ax1.set_xlim(0,pop.x_lim)
    ax1.set_facecolor(color='w')
    ax1.grid(False)

    ax1.set_ylabel('Cells',fontsize=14)
    ax1.tick_params(labelsize=12)    
    
    # if pop.y_lim is not None:
    #     y_lim = pop.y_lim
    # else:
    #     y_lim = np.max(counts) + 0.05*np.max(counts)
    
    xlabels = ax1.get_xticks()
    xlabels = xlabels*pop.timestep_scale
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    xt = ax1.get_xticks()
    ax1.set_xticks(xt)
    ax1.set_xticklabels(xlabels)
    ax1.set_xlabel('Days',fontsize=14)

    # get the length of the simulation
    n_timestep = len(counts[:,0])
    ax1.set_xlim(0,n_timestep)
    ax1.set_ylim(0,10*np.max(counts))

    ax1.set_yscale('symlog',linthresh=1)

    return fig

def plot_fitness_curves(pop,
                        fig_title='',
                        plot_r0 = False,
                        fig=None,
                        ax=None,
                        labelsize=15,
                        linewidth=3,
                        show_legend=True,
                        legend_cols=1,
                        show_axes_labels=True,
                        raw_data = False,
                        color_kwargs={},
                        legend_loc = (1,-0.05),
                        xdata=None,
                        linthresh=10**-3):
    """Plots genotype-specific dose reponse curves (fitness seascape)

    Args:
        pop (population): Population class object
        fig_title (str, optional): Figure title. Defaults to ''.
        plot_r0 (bool, optional): If true, subtracts death rate from fitness. 
        Defaults to False.
        fig (figure, optional): Optional matplotlib figure object to plot to. 
        Defaults to None.
        ax (axes, optional): Optional matplotlib axes object to plot to. 
        Defaults to None.
        labelsize (int, optional): Fontsize of tick labels. Defaults to 15.
        linewidth (int, optional): Plot linewidth. Defaults to 3.
        show_legend (bool, optional): If true, plots the legend. Defaults to True.
        show_axes_labels (bool, optional): If true, adds x- and y-axis labels. Defaults to 
        True.
        raw_data (bool, optional): If true, plots growth rate point estimates rather than 
        estimated dose-response curves. Defaults to False.
        color_kwargs (dict, optional): Options for generating color cycler. Defaults to {}.

    Returns:
        tuple: fig,ax
    """
    if pop.seascape_lib is not None:
        
        if xdata is None:
            xdata = np.logspace(pop.drug_conc_range[0],pop.drug_conc_range[1],
                            num=1000)

        if ax is None:
            fig, ax = plt.subplots(figsize = (10,6))

        cc = gen_color_cycler(**color_kwargs)
        ax.set_prop_cycle(cc)

        if raw_data:
            gl = pop.growth_rate_lib
            for g in range(pop.n_genotype):

                f = gl[str(g)]
                f = [x*60**2 for x in f]
                ax.scatter(xdata,f,label = str(pop.int_to_binary(g)),linewidth=linewidth) 

        
        else:

            if not xdata[0] == 0:

                xdata = np.concatenate(([0],xdata))
            # sl = pop.seascape_lib

            for g in range(pop.n_genotype):

                f = []
                for c in xdata:

                    f.append(fitness.gen_fitness(pop,g,c))

                ax.plot(xdata,f,label = str(pop.int_to_binary(g)),linewidth=linewidth) 

        ax.set_xscale('symlog',linthresh=linthresh)
        
        ax.tick_params(labelsize=labelsize)
        ax.set_ylabel('Growth rate (hr$^-1$)',fontsize=labelsize)
        ax.set_xlabel('Drug concentration (ug/ml)',fontsize=labelsize)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if show_legend:
            ax.legend(fontsize=labelsize,frameon=False,loc=(1,0),ncol=legend_cols)

    else:
        if ax is None:
            fig, ax = plt.subplots(figsize = (10,6))
        
        min_conc = pop.drug_conc_range[0]
        max_conc = pop.drug_conc_range[1]
        if xdata is None:

            conc = np.logspace(min_conc,max_conc,1000)
        else:
            conc = xdata

        cc = gen_color_cycler(**color_kwargs)
        
        ax.set_prop_cycle(cc) 
        
        fit = np.zeros((pop.n_genotype,conc.shape[0]))
        
        for j in range(conc.shape[0]):
            fit[:,j] = fitness.gen_fit_land(pop,conc[j])
        
        if plot_r0:
            fit = fit-pop.death_rate
            ylabel = '$R_{0}$'
            thresh = np.ones(conc.shape)
            ax.plot(conc,thresh,linestyle='dashdot',color='black',linewidth=linewidth)
        
        else:
            ylabel = 'Growth rate (hr$^{-1})$'
        
        for gen in range(pop.n_genotype):
            ax.plot(conc,fit[gen,:],linewidth=linewidth,label=str(pop.int_to_binary(gen)))
        
        if show_legend:
            ax.legend(fontsize=labelsize,frameon=False,loc=legend_loc,ncol=legend_cols)
        
        ax.set_xscale('log')
        
        ax.set_title(fig_title,fontsize=labelsize)
        
        ax.tick_params(labelsize=labelsize)
        
        if show_axes_labels:
            unit = pop.drug_units
            xl = 'Drug concentration (' + unit + ')'
            ax.set_xlabel(xl,fontsize=labelsize)
            ax.set_ylabel(ylabel,fontsize=labelsize)
        # ax.set_frame_on(False)
    
    return fig,ax

def plot_msw(pop,wt,conc=None,fc=None,ncols=2,figsize=(2.5,8)):
    """Plot mutant selection window figures

    Args:
        pop (population): Population class object
        wt (int): wild-tyoe (reference) genotype
        conc (array-like, optional): Drug concentrations. Defaults to None.
        fc (dict, optional): Dict of dose reponse curves. If none, gets data from pop. 
        Defaults to None.
        ncols (int, optional): _description_. Defaults to 2.
        figsize (tuple, optional): _description_. Defaults to (2.5,8).

    Returns:
        _type_: _description_
    """

    if conc is None:
        min_conc = pop.drug_conc_range[0]
        max_conc = pop.drug_conc_range[1]
        conc = np.logspace(min_conc,max_conc,1000)
    if fc is None:
        fc =  fitness.gen_fitness_curves(pop,conc=conc)

    rows = int((pop.n_allele)/ncols)
    fig, ax = plt.subplots(rows,ncols,figsize=figsize)
    
    neighbors = pop.gen_neighbors(wt)
    wt_fitness_curve = fc[wt]

    i = 0

    if ax.ndim == 1:
        if rows > 1:
            for r in range(rows):
            # for col in range(ncols):
                n = neighbors[i]
                wtlabel = pop.int_to_binary(wt) + ' (ref)'
                mutlabel = pop.int_to_binary(n)
                
                ax[r] = plot_msw_to_ax(pop,ax[r],conc,wt_fitness_curve,fc[n],wtlabel,mutlabel)

                ax[r].set_xscale('log')
                ax[r].set_xlim([10**pop.drug_conc_range[0],10**pop.drug_conc_range[1]])
                ax[r].legend(fontsize=10,frameon=False,loc='best')
                shifty(ax[r],i*-0.05)
                i+=1
   

                ax[r].set_ylabel('Replication rate',fontsize=14)
                ax[r].tick_params(axis='both',labelsize=14)

            ax[-1].set_xlabel('Drug concentration \n($\mathrm{\mu}$g/mL)',
                                    fontsize=14)
            
        else:
            for c in range(ncols):
                n = neighbors[i]
                wtlabel = pop.int_to_binary(wt)
                mutlabel = pop.int_to_binary(n)

                ax[c] = plot_msw_to_ax(pop,ax[c],conc,wt_fitness_curve,fc[n],wtlabel,mutlabel)

                ax[c].set_xscale('log')
                ax[c].set_xlim([10**pop.drug_conc_range[0],10**pop.drug_conc_range[1]])
                ax[c].legend(fontsize=12,frameon=False,loc='best')
                shiftx(ax[c],i*0.05)
                i+=1

                ax[c].set_xlabel('Drug concentration \n($\mathrm{\mu}$g/mL)',fontsize=14)
                ax[c].tick_params(axis='both',labelsize=14)

            ax[-1].set_ylabel('Replication rate ($hr^{-1}$)',
                                        fontsize=14)
            ax[-1].yaxis.set_label_position("right")
    else:
        for r in range(rows):
            for col in range(ncols):
                n = neighbors[i]
                wtlabel = pop.int_to_binary(wt) + ' (ref)'
                mutlabel = pop.int_to_binary(n)   

                ax[r,col] = plot_msw_to_ax(pop,ax[r,col],conc,wt_fitness_curve,fc[n],wtlabel,mutlabel)

                ax[r,col].set_xscale('log')
                ax[r,col].set_xlim([10**pop.drug_conc_range[0],10**pop.drug_conc_range[1]])
                ax[r,col].legend(fontsize=10,frameon=False,loc='best')
            i+=1
        
        for r in range(rows):
            # ax[r,0].set_ylabel('$R_{0}$',fontsize=10)
            ax[r,0].set_ylabel('Replication \nrate',fontsize=17)
        # for c in range(rows):
        #     ax[c].set_xlabel('Drug concentration \n($\mathrm{\mu}$M)',
        #                           fontsize=15)
        for col in range(ncols):
            ax[-1,col].set_xlabel('Drug concentration \n($\mathrm{\mu}$g/mL)',
                                    fontsize=20)

    return fig, ax

def plot_timecourse_to_axes(pop,
                        counts,
                        counts_ax,
                        drug_curve=None,
                        drug_curve_label='Drug Concentration \n($\u03BC$M)',
                        drug_curve_legend_label = None,
                        # drug_curve_linestyle='--',
                        drug_ax_sci_notation=False,
                        drug_ax=None,
                        labelsize=15,
                        linewidth=3,
                        legend_labels = True,
                        label_lines = False,
                        select_labels=None,
                        label_xpos = None,
                        grayscale=False,
                        legend_size=8,
                        color_kwargs = {},
                        drug_kwargs = {},
                        label_kwargs={},
                        **kwargs):
    """
    Plots simulation timecourse to user defined axes (counts_ax).

    Parameters
    ----------
    pop : Population class object
        Population class object containing population visualization options.
    counts : numpy array
        Simulation data to be plotted.
    counts_ax : matplotlib axes object
        Axes on which data is plotted.
    drug_curve : numpy array, optional
        Optional drug concentration curve to plot. Requires additional drug
        axes. The default is None.
    drug_ax : matplotlib axes, optional
        Axes on which drug curve is plotted. The default is None.
    labelsize : float, optional
        Font size of the labels. The default is 15.
    linewidth : float, optional
        Width parameter passed to matplotlib plot function. The default is 3.

    Raises
    ------
    Exception
        Error given if no drug axes are provided but the drug curve is not 
        None (drug data needs drug axes to plot to).

    Returns
    -------
    counts_ax : matplotlib axes
        Axes with counts data plotted.
    drug_ax : matplotlib axes
        Axes with drug curve data plotted.

    """
    
    counts_total = np.sum(counts,axis=0)
    sorted_index = counts_total.argsort()
    sorted_index_big = sorted_index[-legend_size:]
    
    if grayscale is False:     
        cc = gen_color_cycler(**color_kwargs)
        counts_ax.set_prop_cycle(cc)
    
    if drug_curve is not None:
        if 'color' in drug_kwargs:
            color = drug_kwargs['color']
        else:
            color='black'
        if drug_ax is None:
            drug_ax = counts_ax.twinx() # ax2 is the drug timecourse
            if drug_curve_label is None:
                yax_label = ''
            else:
                yax_label = drug_curve_label

            drug_ax.set_ylabel(yax_label, 
                            color=color,fontsize=labelsize)
            
        drug_ax.plot(drug_curve,zorder=0,**drug_kwargs)
        
        if pop.drug_log_scale:
            drug_ax.set_yscale('log')
            if min(drug_curve) <= 0:
                axmin = 10**-3
            else:
                axmin = min(drug_curve)
            drug_ax.set_ylim(axmin,2*max(drug_curve))
        else:
            drug_ax.set_ylim(0,1.1*max(drug_curve))
            
        # drug_ax.yaxis.label.set_color('gray')
        drug_ax.tick_params(labelsize=labelsize,color=color)
        plt.setp(drug_ax.get_yticklabels(), color=color)
        if drug_ax_sci_notation:
            drug_ax.ticklabel_format(style='scientific',axis='y',
                                    scilimits=(0,3))
    
    for genotype in range(counts.shape[1]):
        if genotype in sorted_index_big:
            if legend_labels:
                counts_ax.plot(counts[:,genotype],linewidth=linewidth,
                            zorder=10,
                            label=str(pop.int_to_binary(genotype)),
                            **kwargs)
            else:
                counts_ax.plot(counts[:,genotype],linewidth=linewidth,
                            zorder=10,
                            **kwargs)
        else:
            counts_ax.plot(counts[:,genotype],linewidth=linewidth,
                        zorder=10,
                        label=None)
    
    if pop.counts_log_scale:
        counts_ax.set_yscale('log')
        yl = counts_ax.get_ylim()
        yl = [10**1,yl[1]]
        counts_ax.set_ylim(yl)
    
    counts_ax.set_xlim(0,pop.x_lim)
    counts_ax.set_facecolor(color='w')
    counts_ax.grid(False)
    # counts_ax.set_ylabel('Cells',fontsize=20)
    counts_ax.tick_params(labelsize=labelsize)

    xticks = counts_ax.get_xticks()
    xlabels = xticks
    xlabels = xlabels*pop.timestep_scale
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    counts_ax.set_xticks(xticks)
    
    counts_ax.set_xticklabels(xlabels)
    
    xl = [0,len(counts[:,0])]
    counts_ax.set_xlim(xl)
    counts_ax.spines["top"].set_visible(False)
    counts_ax.spines["right"].set_visible(False)
    
    counts_ax.patch.set_alpha(0)
    if drug_ax is not None:
        drug_ax.zorder = 0
    # counts_ax.zorder = 10
    
    if label_lines:
        lines = counts_ax.get_lines()
        
        for i in range(len(label_xpos)):
            sl = select_labels[i]
            labelLine(lines[sl],label_xpos[i],
                    fontsize=5,
                    zorder=100,
                    outline_color='white',
                    outline_width=6,
                    **label_kwargs)
    
    return counts_ax, drug_ax

def plot_landscape(p,conc=10**0,
                arrowprops=None,
                fit_land=None,
                relative=True,
                rank=True,
                ax=None,
                ignore_zero=False,
                colorbar_lim=None,
                colorbar=True,
                node_size = 800,
                node_label = 'base2',
                textsize=11,
                resize_param=0.2,
                square=False,
                arrows=False,
                trajectory_list=None,
                textcolor='black',
                cbax=None,
                cblabel=None,
                cbloc = [0.1,0.8,0.3,0.5],
                network_only=False, # plots just the network without any fit_land data
                edge_color='gray',
                edge_alpha=1,
                plot_sub_network=False,
                sub_network=None,
                sub_network_color='white',
                weight_list=None,
                arrow_alpha_list=None,
                **kwargs):
    """
    Plots a graph representation of this landscape on the current matplotlib figure.
    If p is set to a vector of occupation probabilities, the edges in the graph will
    have thickness proportional to the transition probability between nodes.
    """

    if fit_land is None:
        fit_land = fitness.gen_fit_land(p,conc)
    
    if relative:
        fit_land = fit_land-min(fit_land)
        fit_land = fit_land/max(fit_land)
        
    if ax is None:
        fig,ax=plt.subplots()
        
    if rank:
        fit_land = scipy.stats.rankdata(fit_land)
        if cblabel is None:
            cblabel = 'Rank'
    
    if ignore_zero:
        fit_land_t = [f==0 for f in fit_land]
        fit_land[fit_land==0] = 'NaN'

    if network_only:
        colorbar = False
    
    # Figure out the length of the bit sequences we're working with
    N = int(np.log2(len(fit_land)))

    # Generate all possible N-bit sequences
    genotypes = np.arange(2**N)
    genotypes = [p.int_to_binary(g) for g in genotypes]

    # Turn the unique bit sequences array into a list of tuples with the bit sequence and its corresponding fitness
    # The tuples can still be used as nodes because they are hashable objects
    genotypes = [(genotypes[i], fit_land[i]) for i in range(len(genotypes))]

    # Build hierarchical structure for N-bit sequences that differ by 1 bit at each level
    hierarchy = [[] for i in range(N+1)]
    for g in genotypes: hierarchy[g[0].count("1")].append(g)

    # Add all unique bit sequences as nodes to the graph
    G = nx.DiGraph()
    G.add_nodes_from(genotypes)

    # Add edges with appropriate weights depending on the TM
    TM = p.random_mutations(len(genotypes))
    for i in range(len(TM)):
        for j in range(len(TM[i])):
            if TM[i][j] != 0 and i != j:
                G.add_edge(genotypes[i], genotypes[j], weight=1)
    

    # just using spring layout to generate an initial dummy pos dict
    pos = nx.spring_layout(G)

    # # calculate how many entires in the longest row, it will be N choose N/2
    # # because the longest row will have every possible way of putting N/2 1s (or 0s) into N bits
    maxLen = math.factorial(N) / math.factorial(N//2)**2

    # Position the nodes in a layered hierarchical structure by modifying pos dict
    y = 1
    for row in hierarchy:
        if len(row) > maxLen: maxLen = len(row)
    for i in range(len(hierarchy)):
        levelLen = len(hierarchy[i])
        # algorithm for horizontal spacing.. may not be 100% correct?
        offset = (maxLen - levelLen + 1) / maxLen
        xs = np.linspace(0 + offset / 2, 1 - offset / 2, levelLen)
        for j in range(len(hierarchy[i])):
            pos[hierarchy[i][j]] = (xs[j], y)
        y -= 1 / N
    
    labels = dict(pos)
    for k in labels.keys():
        labels[k] = k[0]
    
    xy = np.asarray([pos[v] for v in list(G)])
    
    # draw edges
    
    edgelist = list(G.edges())
    edge_pos = np.asarray([(pos[e[0]], pos[e[1]]) for e in edgelist])

    edge_colors = []
    edge_widths = []
    if plot_sub_network:
        sub_network_str = [p.int_to_binary(b) for b in sub_network]

    for e in edgelist:

        if plot_sub_network:
            if e[0][0] in sub_network_str and e[1][0] in sub_network_str:
                edge_colors.append('black')
                edge_widths.append(2)
        else:
            edge_colors.append(edge_color)
            edge_widths.append(1)


    edge_collection = LineCollection(
        edge_pos,
        linewidths=edge_widths,
        antialiaseds=(1,),
        linestyle='solid',
        zorder=1,
        color=edge_colors,
        alpha=edge_alpha)
    edge_collection.set_zorder(1)
    ax.add_collection(edge_collection)
 
 
    # Filter edge list and draw arrows: 

    filtered_edges = []
    unique_edges = set()

    for edge in edgelist:
        # Sort the edge tuple to ensure order independence
        sorted_edge = tuple(sorted(edge))
        reversed_edge = tuple(reversed(edge))

        # Check if either the sorted or reversed edge is already in the set
        if sorted_edge not in unique_edges and reversed_edge not in unique_edges:
            unique_edges.add(sorted_edge)
            filtered_edges.append(edge)

    if(arrows):

        for edge in filtered_edges:
            if((edge[0][1]-edge[1][1])>0):
                ax.annotate('', xy=pos[edge[0]], xytext=pos[edge[1]],
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5, shrinkA=14, shrinkB=14))
            
            if((edge[0][1]-edge[1][1])<0):
                ax.annotate('', xy=pos[edge[1]], xytext=pos[edge[0]],
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.5, shrinkA=14, shrinkB=14))


    # #Path trajectory with arrows
                
    # if trajectory and not isinstance(trajectory, list):
    #     raise ValueError("If trajectory parameter is set to True, trajectory data must be provided as a list.")
    
    # draw nodes
    
    if colorbar_lim is not None:
        vmin = colorbar_lim[0]
        vmax = colorbar_lim[1]
    else:
        vmin=min(fit_land)
        vmax=max(fit_land)

    if not network_only:
        ax.scatter(xy[:,0],xy[:,1],
                s=node_size,
                c=fit_land,
                vmin=vmin,
                vmax=vmax,
                clip_on=False,
                **kwargs)
        
    else:
        ax.scatter(xy[:,0],xy[:,1],
                s=node_size,
                c='white',
                edgecolors='black',
                vmin=vmin,
                vmax=vmax,
                clip_on=False,
                **kwargs)
    if plot_sub_network:
        if sub_network is None:
            raise ValueError('sub_network should be a list of ints')
        ax.scatter(xy[sub_network,0],xy[sub_network,1],
                s=node_size,
                c=sub_network_color,
                vmin=vmin,
                vmax=vmax,
                clip_on=False,
                **kwargs)
    
    # if you don't want to include nodes with fitness = 0
    if ignore_zero:
        fit_land_t = np.array(fit_land_t)
        indx = np.argwhere(fit_land_t==True)
        for i in indx:
            ax.scatter(xy[i,0],xy[i,1],
            s=node_size,
            c='gray',
            clip_on=False,
            **kwargs)

    if textcolor is not None:
        for n, label in labels.items():
            (x, y) = pos[n]
            if not isinstance(label, str):
                label = str(label)  # this makes "1" and 1 labeled the same
            
            if node_label == 'base10':
                l = 0
                for i in range(len(label)):
                    l += int(label[i]) * 2**(len(label) - i - 1)

                label = str(l)
            
            if plot_sub_network and n[0] in sub_network_str:
                color = 'white'
            else:
                color = textcolor
            ax.text(
                x,
                y,
                label,
                size=textsize,
                color=color,
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transData,
                clip_on=True,
            )
        
    # display colorbar    
    if colorbar:
        if 'cmap' in kwargs:
            cmap=kwargs['cmap']
            sm = plt.cm.ScalarMappable(cmap=cmap, 
                                        norm=plt.Normalize(vmin = vmin, vmax=vmax))
        else:
            sm = plt.cm.ScalarMappable(norm=plt.Normalize(
                vmin = vmin, vmax=vmax))
        sm._A = []

        cbax = ax.inset_axes(cbloc)
        cbax.set_frame_on(False)
        cbax.set_xticks([])
        cbax.set_yticks([])
        
        cb = plt.colorbar(sm,
                        drawedges=False,
                        ax=cbax,
                        location='right',
                        aspect=10)
        cb.outline.set_visible(False)
        cb.set_label(cblabel,fontsize=12)
        
        if rank:
            ticks = [min(fit_land),max(fit_land)]
            cb.set_ticks(ticks)
            ticks = [max(fit_land),min(fit_land)]
            ticks = np.array(ticks).astype('int')
            ticks = [str(t) for t in ticks]
            cb.set_ticklabels(ticks,fontsize=10)
    
    if square:
        ydata_range = max(xy[:,1])-min(xy[:,1])
        xdata_range = max(xy[:,0])-min(xy[:,0])
        ax.set_aspect(xdata_range/ydata_range)
        
    xl = ax.get_xlim()
    xrange = xl[1]-xl[0]
    xl = [xl[0]-resize_param*xrange,xl[1]+xrange*resize_param]
    ax.set_xlim(xl)
    
    yl = ax.get_ylim()
    yrange = yl[1]-yl[0]
    yl = [yl[0]-resize_param*yrange,yl[1]+yrange*resize_param]
    ax.set_ylim(yl)


    if trajectory_list is not None:
        for i,trajectory in enumerate(trajectory_list):
                
                if weight_list is not None:
                    weight = weight_list[i]
                    arrowprops['lw'] = weight
                
                if arrow_alpha_list is not None:
                    alpha = arrow_alpha_list[i]
                    arrowprops['alpha'] = alpha

                trajectory_pairs = []

                if arrowprops is None:
                    arrowprops = dict(arrowstyle='->', color='black', lw=2.5, shrinkA=14.5, shrinkB=14.5, mutation_scale=15)

                for i in range(len(trajectory) - 1):
                    trajectory_pairs.append([trajectory[i], trajectory[i + 1]])

                binary_trajectory_pairs = []

                for pair in trajectory_pairs:
                    binary_pair = []
                    for node in pair:
                        binary_representation = format(node, '04b')  # Convert node number to 4-bit binary representation
                        binary_pair.append("" + binary_representation + "")
                    binary_trajectory_pairs.append(binary_pair)

                for pair in binary_trajectory_pairs:
                    start_binary = pair[0]
                    end_binary = pair[1]
            
                    # Find the corresponding edge in filtered_edges based on the binary values
                    for edge in filtered_edges:
                        if edge[0][0] == start_binary:
                            start_pos = pos[edge[0]]
                            break
            
                    for edge in filtered_edges:
                        if edge[0][0] == end_binary:
                            end_pos = pos[edge[0]]
                            break

                    # For genotypes only in the end position in filtered_edges (genotype 15)
                    for edge in filtered_edges:
                        if edge[1][0] == start_binary:
                            start_pos = pos[edge[1]]
                            break

                    for edge in filtered_edges:
                        if edge[1][0] == end_binary:
                            end_pos = pos[edge[1]]
                            break
        
                # Draw trajectory arrows
                    
                    ax.annotate('', xy=end_pos, xytext=start_pos,
                            arrowprops=arrowprops)
    
    ax.set_axis_off()
    return ax

def add_landscape_to_fitness_curve(c,ax,pop,
                                vert_lines=True,
                                position = 'top',
                                pad = 0,
                                ypos=1,
                                vert_lines_ydata = None,
                                width=3,
                                height=None,
                                vert_lines_kwargs={},
                                **kwargs):

    if position == 'top':
        ypos = ypos+pad
    elif position == 'bottom':
        ypos = -ypos-pad
    else:
        raise Exception('Position argument not recognized')
    
    if height is None:
        height = width

    ax.set_clip_on(False)

    x = get_pos_in_log_space(c, width)
    lax = ax.inset_axes([x[0],ypos,x[1]-x[0],height],transform=ax.transData)
    lax = plot_landscape(pop,conc=c,ax=lax,
                        **kwargs)
    
    if vert_lines:
        if vert_lines_ydata is None:
            yl = ax.get_ylim()
            ydata = np.arange(yl[0],yl[1],0.1)
        else:
            ydata = vert_lines_ydata
        xdata = np.ones(len(ydata))*c
        ax.plot(xdata,ydata,'--',color='black',**vert_lines_kwargs)        
    
    return ax,lax

def plot_population_count(pop,
                        c,
                        ax=None,
                        thresh=None,
                        normalize=False,
                        max_cells=None,
                        logscale=True,
                        **kwargs):
    if ax is None:
        fig,ax = plt.subplots(figsize=(6,4))
    if thresh is None:
        thresh = pop.max_cells/10
    c1 = [245,100,100]
    c1 = [c/255 for c in c1]
    c2 = [100,100,245]
    c2 = [c/255 for c in c2]
    if c[-1] < thresh:
        if normalize:
            c = c/pop.max_cells
        ax.plot(c,color=c2,label='extinct',**kwargs)
    else:
        if normalize:
            c = c/pop.max_cells
        ax.plot(c,color=c1,label='resistant',**kwargs)
    
    xticks = ax.get_xticks()
    xlabels = xticks
    xlabels = xlabels*pop.timestep_scale
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    
    if logscale:
        ax.set_yscale('log')
    
    return ax

def plot_kaplan_meier(pop,
                    event_times,
                    label=None,
                    t_max=None,
                    t_vect = None,
                    n_sims=None,
                    ax=None,
                    mode='resistant',
                    errorband=True,
                    **kwargs):
    
    if t_max is None:
        t_max = int(max(event_times)) # hours
    if n_sims is None:
        n_sims = pop.n_sims
        
    survival_curve = np.ones(t_max)*100
    
    for t in range(len(survival_curve)-1):
        if t>0:
            survival_curve[t] = survival_curve[t-1]
        if any(event_times==t):
            num = np.argwhere(event_times==t)
            num = num.shape[0]
            perc = 100*num/n_sims
            survival_curve[t] = survival_curve[t]-perc

    survival_curve[-1] = survival_curve[-2]
    if ax is None:
        fig,ax = plt.subplot(figsize=(5,7))
    
    if mode == 'resistant':
        survival_curve = 100-survival_curve
        ylabel='% resistant'
    else:
        ylabel='% survival'
    
    if t_vect is None:
        t_vect = np.arange(t_max)
    ax.plot(t_vect,survival_curve,label=label,**kwargs)
    
    if errorband:
        # compute error bars
        # rule of succession explanation: https://en.wikipedia.org/wiki/Rule_of_succession
        err = np.zeros(t_max)
        for t in range(t_max):
            p = (np.array(survival_curve[t]/100) + 1)/(n_sims + 2) # uniform prior (rule of succession) 
            n = n_sims
            q = 1-p
            # standard deviation of the estimator of the parameter of a binomial distribution
            
            err[t] = 100*((p*q)/n)**0.5 #
    t = np.arange(t_max)
    
    ax.fill_between(t_vect,survival_curve-err,survival_curve+err,alpha=0.4)
    
    xticks = ax.get_xticks()
    xlabels = xticks
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    
    xl = [0,len(survival_curve)]
    ax.set_xlim(xl)
    ax.set_ylim([0,100])
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Days')
    return ax

def find_zero_crossing(v):

    indx = None

    for i in range(len(v)-1):
        if v[i]*v[i+1] < 0:
            indx = i + 1
            break

    return indx

def get_msw(ref_fitness_curve,cur_fitness_curve):
    """Computes selection windows for two dose-response curves

    Args:
        ref_fitness_curve (array-like): reference dose response curve
        cur_fitness_curve (array-like): mutant dose response curve

    Returns:
        dict: dict of selection windows
    """

    chunks = {}

    left_indx = 0

    right_indx = len(ref_fitness_curve)

    ref_fitness_curve = np.array(ref_fitness_curve)
    cur_fitness_curve = np.array(cur_fitness_curve)

    # find the death window

    if any(ref_fitness_curve<0):
        ref_death = np.argwhere(ref_fitness_curve<0)
        ref_death = ref_death[0][0]
        if any(cur_fitness_curve<0):
            cur_death = np.argwhere(cur_fitness_curve<0)
            cur_death = cur_death[0][0]

            left_death_window_indx = np.max([cur_death,ref_death])
            right_death_window_indx = len(ref_fitness_curve)
            chunks['net loss'] = [(left_death_window_indx,right_death_window_indx)]
        else:
            # left_death_window = max(conc)+1
            left_death_window_indx = len(ref_fitness_curve)
    else:
        # left_death_window = max(conc)+1
        left_death_window_indx = len(ref_fitness_curve)

    while left_indx < len(ref_fitness_curve):

        rfc_t = ref_fitness_curve[left_indx:]
        cfc_t = cur_fitness_curve[left_indx:]

        indx = find_zero_crossing(rfc_t-cfc_t)

        if indx is not None:
        # if any(rfc_t-cfc_t == 0): # this means there is a crossing
            right_indx = indx
        else:
            # this must be the last chunk before net loss
            right_indx = left_death_window_indx - left_indx -1

            if rfc_t[right_indx-1] > cfc_t[right_indx-1]:
                label = 'reference'
            else:
                label = 'mutant'
            
            if label not in chunks.keys():
                chunks[label] = []
            # right_conc = conc[right_indx+left_indx]
            # left_conc = conc[left_indx]
            chunks[label].append((left_indx,right_indx+left_indx))
            left_indx = right_indx + left_indx + 1

            break

        # check to make shure this crossing is before there is net death
        if right_indx+left_indx > left_death_window_indx:
            break
        else: # figure out if ref selection or mutant selection
            if rfc_t[right_indx-1] > cfc_t[right_indx-1]:
                label = 'reference'
            else:
                label = 'mutant'
            
            if label not in chunks.keys():
                chunks[label] = []
            # right_conc = conc[right_indx+left_indx]
            # left_conc = conc[left_indx]
            chunks[label].append((left_indx,right_indx+left_indx))
            left_indx = right_indx + left_indx + 1

    return chunks

def plot_msw_to_ax(pop,ax,conc,wt_fitness_curve,mut_fitness_curve,wtlabel,mutlabel):
    
    ax.plot(conc,wt_fitness_curve,color='black',
                        label=wtlabel,linewidth=3)
                
    ax.plot(conc,mut_fitness_curve,color='white',
            label=mutlabel,linewidth=3)

    chunks = get_msw(wt_fitness_curve,mut_fitness_curve)

    
    r_d = pop.death_rate + 0.1

    for key in chunks:
        # label = key
        if key == 'net loss':
            color = '#e41a1c'
            label = 'net loss'
        elif key == 'reference':
            color = '#ff7f00'
            label = 'reference selection'
        else:
            color = 'tab:blue'
            label = 'mutant selection'

        for c in chunks[key]:
            # ax.axvspan(conc[c[0]],conc[c[1]-1],facecolor=color,alpha=0.7)
            start_x = c[0]
            end_x = c[1]-1

            if label == 'mutant selection':
                s = np.array((mut_fitness_curve+r_d)/(wt_fitness_curve+r_d))
                s = s[start_x:end_x]
                s = s-np.min(s)
                s = s/np.max(s)
            elif label == 'reference selection':
                s = np.array((wt_fitness_curve+r_d)/(mut_fitness_curve+r_d))
                s = s[start_x:end_x]
                s = s-np.min(s)
                s = s/np.max(s)
            else:
                s = np.ones(end_x-start_x)

            indx = 0

            for x in range(start_x,end_x):
                # x_rect = [conc[x],conc[x],conc[x+1],conc[x+1]]
                # if s[indx] == 1 and g == 0 and n == 1:
                #     ax.fill(x_rect,y,color,alpha=s[indx],label=label)
                # else:

                ax.axvspan(conc[x],conc[x+1],facecolor=color,alpha=s[indx])
                indx += 1

            # draw a box

    # ax.axvspan(msw_left, msw_right, 
    #         facecolor='#2ca02c',alpha=0.7)
    #         # label='MSW')
    # ax.axvspan(min(conc),msw_left, 
    #         facecolor='#ff7f00',alpha=0.7)
    #         # label='MSW')
    # ax.axvspan(msw_right,max(conc), 
    #         facecolor='#e41a1c',alpha=0.7)
    #         # label='MSW')
    return ax

def msw_grid(pop,genotypes,
            ax=None,
            legend=True,
            labelsize=10,
            ticklabelsize=10,
            comp_annotate_pos=10**-4.4,
            legendloc='best'):

    min_conc = pop.drug_conc_range[0]
    max_conc = pop.drug_conc_range[1]
    conc = np.logspace(min_conc,max_conc,1000)
    fc = fitness.gen_fitness_curves(pop,conc=conc)

    # get the row height
    n_rows = len(genotypes)*pop.n_allele + len(genotypes)
    h = 1/n_rows
    ylevel = 1 # top-down illustration

    # initialize a figure with appropriate aspect ratio
    width = 6 # inches
    row_height = 0.2 # inches

    r_d = pop.death_rate + 0.1

    if ax is None:
        fig,ax = plt.subplots(figsize=(width,row_height*n_rows))

    for g in genotypes:
        # g is the reference genotype
        # get neighbors
        neighbors = pop.gen_neighbors(g)
        # add a line for the reference label
        label = 'Reference = ' + pop.int_to_binary(g)

        pos = (10**3,ylevel-(0.8*h))
        ax.annotate(label,pos,xycoords='data',annotation_clip=True,fontsize=labelsize)

        ylevel += -h

        for n in neighbors:

            chunks = get_msw(fc[g],fc[n])

            y = [ylevel-h,ylevel,ylevel,ylevel-h]


            for key in chunks:
                # label = key
                if key == 'net loss':
                    color = '#e41a1c'
                    label = 'net loss'
                elif key == 'reference':
                    color = '#ff7f00'
                    label = 'reference selection'
                else:
                    color = 'tab:blue'
                    label = 'mutant selection'

                for c in chunks[key]:
                    start_x = c[0]
                    end_x = c[1]-1

                    if label == 'mutant selection':
                        s = np.array((fc[n]+r_d)/(fc[g]+r_d))
                        s = s[start_x:end_x]
                        s = s-np.min(s)
                        s = s/np.max(s)
                    elif label == 'reference selection':
                        s = np.array((fc[g]+r_d)/(fc[n]+r_d))
                        s = s[start_x:end_x]
                        s = s-np.min(s)
                        s = s/np.max(s)
                    elif label == 'net loss':
                        s = np.ones(end_x-start_x)

                    indx = 0

                    for x in range(start_x,end_x):
                        x_rect = [conc[x],conc[x],conc[x+1],conc[x+1]]
                        # if s[indx] == 1 and g == 0 and n == 1:
                        #     ax.fill(x_rect,y,color,alpha=s[indx],label=label)
                        # else:
                        ax.fill(x_rect,y,color,alpha=s[indx])
                        indx += 1

                    # draw a box

                    poly_coords = [(conc[start_x],y[0]),(conc[end_x],y[0]),
                                (conc[end_x],y[1]),(conc[start_x],y[1])]

                    ax.add_patch( Polygon(poly_coords,edgecolor='black',
                                        facecolor=None, closed=True, 
                                        clip_on=False,fill=False) )


            label = pop.int_to_binary(n)
            pos = (comp_annotate_pos,ylevel-(0.8*h))
            ax.annotate(label,pos,xycoords='data',annotation_clip=False,fontsize=labelsize)

            ylevel += -h
    
    ax.set_frame_on(False)
    ax.set_xscale('log')
    ax.set_xlim([10**pop.drug_conc_range[0],10**pop.drug_conc_range[1]])
    ax.tick_params(axis='x', which='major', labelsize=ticklabelsize)
    ax.set_yticks([])
    if legend:
        p1 = ax.plot([1,1], color = '#ff7f00',label = 'reference selection',linewidth=3)
        p2 = ax.plot([1,1], color = 'tab:blue',label = 'mutant selection',linewidth=3)
        p0 = ax.plot([1,1], color = '#e41a1c',label = 'net loss',linewidth=3)
        ax.legend(frameon=False,ncol=3,loc=legendloc)

        lines = ax.get_lines()
        lines[-1].remove()
        lines[-2].remove()
        lines[-3].remove()
    return ax


# Helper methods

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

def x_ticks_to_days(pop,ax):
    
    xticks = ax.get_xticks()
    xlabels = xticks
    xlabels = xlabels*pop.timestep_scale
    xlabels = xlabels/24
    xlabels = np.array(xlabels).astype('int')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    
    x_max = pop.n_timestep
    ax.set_xlim(0,x_max)
    
    return ax

def shiftx(ax,pad):
    """
    shiftx: shifts ax in the x direction by pad units

    Parameters
    ----------
    ax : matplotlib axes
        axes whose position is shifted
    pad : float
        shift amount

    Returns
    -------
    ax : matplotlib axes
        shifted axes

    """
    pos = ax.get_position()
    pos.x0 = pos.x0+pad
    pos.x1 = pos.x1+pad
    ax.set_position(pos)
    
    return ax

def shifty(ax,pad):
    """
    shifty: shifts ax in the y direction by pad units

    Parameters
    ----------
    ax : matplotlib axes
        axes whose position is shifted
    pad : float
        shift amount

    Returns
    -------
    ax : matplotlib axes
        shifted axes

    """
    pos = ax.get_position()
    pos.y0 = pos.y0+pad
    pos.y1 = pos.y1+pad
    ax.set_position(pos)
    
    return ax

def shrinky(ax,pad):
    """
    shrinky: shrinks axes by pad by moving the top border down

    Parameters
    ----------
    ax : matplotlib axes
        axes whose position is shifted
    pad : float
        shrink amount

    Returns
    -------
    ax : matplotlib axes
        shrunk axes

    """
    
    pos = ax.get_position()
    pos.y1 = pos.y1-pad
    ax.set_position(pos)    
    
    return ax