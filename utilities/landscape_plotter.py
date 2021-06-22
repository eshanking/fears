# plots a n-bit (2^n-element) fitness landscape without edges
import numpy as np
import matplotlib.pyplot as plt
import math

###############################################################################
# User inputs

# Input the fitness vector here.
# Must be a list of numbers of length 2^n. Will be normalized to 1.
fitness = [0.71,0.75,0.83,0,
           0.76,0.81,0.96,0.84,
           0.67,0.74,0.89,0.67,
           0.75,0.78,1,0.86]

# You may have to adjust these if the fitness vector is not size 16
circle_size = 600 # size of the circles representing genotypes
font_size = 8 # size of the genotype labels
xy_offset = [-.018,-0.01] # genotype label position offset (to cener them on the circles)

savename = 'landscape' # name for saving the figure

colormap='cividis'

###############################################################################

"""
Method for computing the layer in pascal's triangle that corresponds to the 
number of bits. This gives us the number of layers of the landscape and the 
number of mutants per row.
"""
def pascal(n):
    line = [1]
    for k in range(n):
        line.append(line[k]* (n-k) / (k+1))
    return line

# Method for converting integers to bit strings
def int_to_binary(num,n_allele):
    pad = int(math.log(n_allele,2))
    return bin(num)[2:].zfill(pad)

# Compute the number of possible mutations (length of bit string)
n = int(math.log(len(fitness),2))

# Calculate total number of alleles
n_allele = 2**n

# Compute the row of pascal's triangle that corresponds to the number of mutations 
p = pascal(n)

x_pos = 0
# Empty x and y position vectors
y = []
x = []

# Generate the x,y corrdinates of the circles
for i in range(len(p)):
    y_pos = (np.arange(p[i]) + 1) / (p[i] + 1)
    y.extend(y_pos)
    x.extend(x_pos*np.ones(int(p[i])))
    x_pos += 0.1
    
# Generate a figure and axis
fig,ax = plt.subplots()

# Create a scatter plot of the fitness landscape where the color of the spheres corresponds to the fitness
sc = ax.scatter(x,y,c=fitness,s=circle_size,cmap=colormap)
ax.axis('off')

# Annotate the spheres with the corresponding bit string
# x- and y- offsets are hacked for now, need a more elegant solution in the future
for i in range(2**n):
    ax.annotate(int_to_binary(i,n_allele), (x[i]+xy_offset[0],y[i]+xy_offset[1]),fontsize=font_size)

# Show the colorbar
cb = plt.colorbar(sc,pad=0)
cb.outline.set_visible(False)
cb.set_ticks([0,0.25,0.50,0.75,1])

# Save
savename_t = savename + '.svg'
plt.savefig(savename_t,bbox_inches="tight")