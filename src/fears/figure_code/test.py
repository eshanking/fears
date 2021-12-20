import matplotlib.pyplot as plt

fig,ax = plt.subplots(nrows=2,ncols=2,constrained_layout=True)

for i in range(2):
    for j in range(2):
        ax[i,j].plot([1,2,3])

for i in range(2):
    for j in range(2):
        ax[i,j].set_box_aspect(2)
       