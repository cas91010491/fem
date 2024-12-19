import imageio
from os import listdir, walk
from os.path import isfile, join
from pdb import set_trace


mypath = "2d_bfgs_ToVisualizeImages/OUTPUT_202411092128pseudo2d_plastic_BFGS_15/plots/"

# set_trace()

images = []
formats = ['png', 'jpeg']

fig_numb = 61
filelist = range(1,819,4)

# for num in range(100):
for num in filelist:
    f= "fig"+str(fig_numb)+"-"+str(num)+".png"
    images.append(imageio.imread(mypath+f))
# images = [imageio.imread(mypath+f) for f in listdir(mypath) if f[-3:] in formats]

imageio.mimsave(mypath + f"a_minimz_inc{fig_numb}.gif", images, fps=48)  # Try fps=1000 for a faster video
