import sys
import os
import pickle
cwd = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cwd+"/..")      #this will actually refer to the good FEAssembly.py file

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from pdb import set_trace
from PyClasses import FEModel
from PyClasses.Utilities import plot_coords


def plot4D(model,TIMES,UU,SED, ax=None, undef=False):
    from matplotlib.colors import Normalize

    sed = np.array(SED)
    smin,smax = [sed.min(),sed.max()]
    sed = (sed-smin)/(smax/4-smin)            # normalizing everything between 0 and 1

    if ax is None:
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider, Button
        from matplotlib import cm


        def update_inset_axis(elev, azim):
            axin2.view_init(elev=elev, azim=azim)
            fig_inset.canvas.draw_idle()  # Manually redraw the inset axis

        def update_view(event):
            elev, azim = ax.elev, ax.azim
            update_inset_axis(elev, azim)

        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')



        limits = np.zeros((3,2))
        for body in model.bodies:
            x_all = body.X + UU[:,body.DoFs]
            for i in range(3):
                if x_all[:,:,i].min() < limits[i,0]:
                    limits[i,0] = x_all[:,:,i].min()
                if x_all[:,:,i].max() > limits[i,0]:
                    limits[i,1] = x_all[:,:,i].max()
            
        limits[2][0] = 1.0

        ax.set_xlim3d(limits[0,0], limits[0,1])
        ax.set_ylim3d(limits[1,0], limits[1,1])
        ax.set_zlim3d(limits[2,0], limits[2,1])
        ax.set_aspect("equal")
        ax.text2D(-0.05, 0.7, "time: "+str(round(0.0,8)), transform=ax.transAxes,backgroundcolor=(1.0,1.0,1.0,0.7),zorder=1000000,fontsize='x-large')
        ax.axis('off')

        axin2 = ax.inset_axes([-.08, -.022, 0.02, 0.02], projection='3d',transform = ax.transData,zorder=1000000)
        # axin2 = ax.inset_axes([-.08, 0.0, 0.02, 0.02], projection='3d',transform = ax.transData,zorder=1000000)
        axin2.axis('off')
        axin2.set_xlim3d(0, 1)
        axin2.set_ylim3d(0, 1)
        axin2.set_zlim3d(0, 1)


        if model.transform_2d is not None:
            ax.set_proj_type('ortho')
            ax.view_init(elev=0, azim=-90, roll=0)


        else:
            ax.view_init(elev=5, azim=-60)

        axin2.view_init(elev=ax.elev, azim=ax.azim)

        fig_inset = plt.gcf()  # Get the figure for the inset axis
        fig_inset.canvas.mpl_connect('motion_notify_event', update_view)


    fixed_obj = model.bodies[1].surf.plot(ax,u= None,sed = None,ref=10,onlyMaster=True)
    # set_trace()
    
    init_time = 0
    u = UU[0]
    sedi = sed[0]
    # model.bodies[0].surf.plot(ax,u= u,sed = sedi,ref=1)
    model.bodies[0].plot(ax,u= u,sed = sedi,ref=10)
    model.bodies[0].plot(ax,u= u,sed = None,ref=1)

    # coords_obj = plot_coords(ax,orig=(2.0,-3.0,-5.0))
    coords_obj = plot_coords(axin2)

    cmap = cm.get_cmap('jet')
    # cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=smin,vmax=smax/4)),
    #                      ax=ax, shrink=0.5, aspect=15,extend='max')
    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=smin,vmax=smax/4)),
                         ax=ax, shrink=0.2, aspect=30,extend='max',orientation='horizontal',
                          location='bottom',pad = 0.0, anchor = (0.5,3.4))
    # cbar.set_label('Strain Energy Density')
    cbar.set_label('Plastic Multiplier')

    axtime = fig.add_axes([0.25, 0.1, 0.65, 0.03])

    time_slider = Slider(
        ax=axtime,
        label='time',
        valmin=0.0,
        valmax=0.99,
        valstep=0.001,
        valinit=init_time,
    )


    # The function to be called anytime a slider's value changes
    def update(UU,sed):
        xl = ax.get_xlim()
        yl = ax.get_ylim()
        zl = ax.get_zlim()
        az, el = ax.azim, ax.elev
        
        remove_plots(ax.collections)
        remove_plots(ax.texts)

        idx = round(len(UU)*time_slider.val)
        # print("idx =",idx)
        time = TIMES[idx]
        ax.text2D(-0.05, 0.7, "time: "+str(round(time,8)), transform=ax.transAxes,backgroundcolor=(1.0,1.0,1.0,0.7),zorder=1000000,fontsize='x-large')
        u = UU[idx]
        sedi = sed[idx]
        model.bodies[0].plot(ax,u= u,sed = sedi,ref=10)
        model.bodies[0].plot(ax,u= u,sed = None,ref=1)

        
        ax.set_xlim(xl)
        ax.set_ylim(yl)
        ax.set_zlim(zl)
        ax.azim = az
        ax.elev = el


    def remove_plots(obj_list):
        for obj in obj_list:
            if obj not in fixed_obj and obj not in coords_obj:
                obj.remove()

    # register the update function with each slider
    # time_slider.on_changed(update)
    time_slider.on_changed(lambda val: update(UU, sed))


    # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
    resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')


    def reset(event):
        time_slider.reset()
    button.on_clicked(reset)

    plt.tight_layout()

    plt.show()


# Load the FEModel object
# model = pickle.load(open(cwd + "/Model.dat", "rb"))

# model = pickle.load(open(cwd + "/Model.dat", "rb"))


with open(cwd + "/Model.dat", "rb") as file:
    model = pickle.load(file)

ptt = model.bodies[1]
ndofs = 3*(len(model.bodies[0].X)+len(ptt.X))
ptt.surf.ComputeGrgPatches(np.zeros(ndofs),model.contacts[0].masterNodes,exactNodesGiven=True)

data = np.loadtxt(cwd+'/data_u.csv', delimiter=',')
TIMES = data[:,0]
UU = data[:,1:]

SED = []
EPCUM = []
# epcum_gauss = np.loadtxt(cwd+'/EPcum.csv', delimiter=',').reshape(len(UU),-1,8)
for incr,u in enumerate(UU):
    SED.append(model.bodies[0].get_nodal_SED(u))
    # EPCUM.append(model.bodies[0].get_nodal_EPCUM(epcum_gauss[incr]))

plot4D(model,TIMES,UU,SED, ax=None, undef=False)
# plot4D(model,TIMES,UU,EPCUM, ax=None, undef=False)

