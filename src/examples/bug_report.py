import matplotlib.pyplot as plt
import numpy as np


# def make_twin_axes(ax, *args, **kwargs):
#     """Make a twinx axes of ax. This is used for twinx and twiny."""
#     # Typically, SubplotBase._make_twin_axes is called instead of this.
#     if 'sharex' in kwargs and 'sharey' in kwargs:
#         raise ValueError("Twinned Axes may share only one axis")
#     ax2 = ax.figure.add_axes(ax.get_position(True))
#     ax.set_adjustable('datalim')
#     ax2.set_adjustable('datalim')
#     ax._twinned_axes.join(ax, ax2)
#     return ax2

fig,ax = plt.subplots(figsize=(5,5))
# fig.draw_no_output()

# create a horizontal bar chart
barx = [1,2,3]
bary = np.array([1,2,3])/3

rects = ax.barh(barx,bary)

# create some data to plot
x = np.arange(50)

insets = [] # list to hold inset axes

# create inset axes at each level of the bar chart
for num in range(len(barx)):
    ypos = rects.patches[num].get_y()
    xpos = 1.2
    
    ax1 = ax.inset_axes([xpos,ypos,0.4,0.8],transform=ax.transData)
    
    ax1.plot(x,x**num,linewidth=2)
    insets.append(ax1)

# twin_ax = insets[0].twin_ax.set_in_layout(False)

# create twin axes for each inset axes
twins = []
for a in insets:
    twin_ax = a.twinx()
    twin_ax.set_in_layout(False)
    twins.append(twin_ax)