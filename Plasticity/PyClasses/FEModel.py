from PyClasses.BoundaryConditions import *
from PyClasses.Utilities import float_to_fraction, printif, checkVarsSize, flatList,plot_coords,quadratic_fit_min_zeros
import PyClasses.FEAssembly
from scipy import sparse
from scipy.sparse.linalg import spsolve, cgs, bicg, gmres, lsqr
from scipy.linalg import solve as slsolve
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import time, pickle, sys, os, csv

# import petsc4py
# petsc4py.init(sys.argv)
# from petsc4py import PETSc

from pdb import set_trace

class FEModel:
    def __init__(self,bodies, contacts, BoundaryConditions, UpdateOnIteration = True,transform_2d=None):
        self.Mainname = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        self.bodies = bodies
        self.contacts = contacts
        self.UOI = UpdateOnIteration
        self.transform_2d = transform_2d
        
        self.SED = [[] for _ in range(len(bodies))]     # stores nodal Strain-energy density for every successful increment
        X , DoFs, n0 = [], [], 0
        for bid, body in enumerate(bodies):
            X.extend(body.X)
            if not hasattr(self, "dim"):
                self.dim = len(X[0])
            nf = n0 + len(body.X)
            body.DoFs = np.array( [ [ self.dim*n + i for i in range(self.dim) ] for n in range(n0,nf) ] )
            body.id = bid
            # body.surf.DoFs = body.DoFs[body.surf.nodes]
            for sub_surf in body.surf.sub_surfs:
                # sub_surf.DoFs =  body.DoFs[sub_surf.nodes]
                sub_surf.DoFs =  body.DoFs
            body.surf.id = bid
            n0 = nf

        X = np.array( X, dtype=float)
        self.PlotBoundaries = [[min(X[:,0]), max(X[:,0])],
                               [min(X[:,1]), max(X[:,1])],
                               [min(X[:,2]), max(X[:,2])]]
        self.X = X.ravel()
        self.ndof = len(self.X)

        self.u = np.zeros_like(self.X,dtype=float)
        self.u_temp = np.array(self.u)
        self.fint = np.zeros_like(self.u,dtype=float)
        self.fext = np.zeros_like(self.u,dtype=float)

        self.BCs = InitializeBCs(BoundaryConditions)    # Creates BC objects with their attributes (dofs, vals, times, etc)


        current_time=time.strftime("%Y%m%d%H%M", time.localtime())
        dir =  "OUTPUT_"+current_time+self.Mainname+"/"
        if not os.path.exists(dir):
            os.mkdir(dir)
        self.output_dir = dir
        pickle.dump(self,open(dir+"Model.dat","wb"))


    def printContactStates(self,only_actives=True,veredict=False):
        for i_c, contact in enumerate(self.contacts):
            print("contact",i_c,":")
            contact.printStates(only_actives=only_actives,veredict=veredict)

    def plot(self, ax, ref = 10, specialPatches = None, almostSpecialPatches = None,
              specialNodes = None, time = None, fintC = False, plotHooks=False,OnlyMasterSurf=False):
        for body in self.bodies:
            if body in [contact.slaveBody for contact in self.contacts]:
                body.surf.plot(ax, self.u, specialNodes = specialNodes,ref=ref)
            elif body in [contact.masterBody for contact in self.contacts]:
                body.surf.plot(ax, self.u, specialPatches = specialPatches, almostSpecialPatches = almostSpecialPatches,ref=ref,onlyMaster=OnlyMasterSurf)
            else:
                body.surf.plot(ax, self.u,ref=ref)

        if fintC:
            for contact in self.contacts: contact.plotForces(self.u, ax)
        if plotHooks:
            for contact in self.contacts: contact.plotHooks(ax)


        if time != None:
            ax.text2D(0.05, 0.95, "time: "+str(round(time,6)), transform=ax.transAxes)

    def plotNow(self, fintC = False, as2D= False,OnlyMasterSurf=False):
        #PLOT OPTIONS
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_proj_type('ortho')
        # plt.axis('off')

        if as2D: ax.view_init(elev=0, azim=-90, roll=0)

        (minx,maxx),(miny,maxy),(minz,maxz) = self.PlotBoundaries
        hg = max([dij[1]-dij[0] for dij in self.PlotBoundaries])/2    #half-gap to be used in every direction
        cx,cy,cz = (minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2
        
        plot_coords(ax,orig=(minx-0.5,miny-0.5,minz-0.5))

        ax.set_xlim3d((cx-hg, cx+hg))
        ax.set_ylim3d((cy-hg, cy+hg))
        ax.set_zlim3d((cz-0.8*hg, cz+0.8*hg))

        # candidatePatches = sorted(set(flatList(self.contacts[0].candids))-set(self.contacts[0].actives))
        # activeNodes = [self.contacts[0].slaveNodes[idx] if self.contacts[0].proj[idx,3]!=0 else None for idx in range(self.contacts[0].nsn)]
        activeNodes = []

        # self.plot(ax, ref = 4,specialPatches = [(0,1,0,0.7),self.contacts[0].actives],almostSpecialPatches = [(1,0.6,0.0,0.7),candidatePatches],specialNodes = ["blue", activeNodes], fintC= fintC)
        self.plot(ax, ref = 10,specialPatches = None,almostSpecialPatches = None,specialNodes = {"blue": activeNodes}, fintC= fintC,OnlyMasterSurf=OnlyMasterSurf)
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.show()

    def plotContactNow(self, ctct_id=0, labels =True, block = True,SlaveQuads=False, context = False):
        """Enable this to visualize contact state"""
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        self.contacts[ctct_id].plotContact(self.u_temp,ax,labels=labels,SlaveQuads=SlaveQuads, context=context)
        self.contacts[ctct_id].plotHooks(ax)
        ax.set_zlim3d(-2, 2)
        plt.show(block = block)

    def savefig(self, increment, iteration = None,  elevation=[ 30 , 30], azimut=[-135 , -135], distance = [10, 10], times = None, dpi = 200, fintC = False, Hooks= False):
        #PLOT OPTIONS
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')


        if not os.path.exists(self.output_dir+"plots"):
            os.mkdir(self.output_dir+"plots")


        #for good appeareance try to keep a ratio  6range(x)=6range(y)=5range(z)
        (minx,maxx),(miny,maxy),(minz,maxz) = self.PlotBoundaries
        hg = max([dij[1]-dij[0] for dij in self.PlotBoundaries])/2    #half-gap to be used in every direction
        cx,cy,cz = (minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2

        ax.set_xlim3d((cx-hg, cx+hg))
        ax.set_ylim3d((cy-hg, cy+hg))
        ax.set_zlim3d((cz-0.8*hg, cz+0.8*hg))

        AlmostImportant = sorted(set(flatList(self.contacts[0].candids))-set(self.contacts[0].actives)) if len(self.contacts)!= 0 else []

        if times!= None:
            t0,t,tf = times
        else:
            t = None

        ctct = self.contacts[0]
        # activeNodes = [ctct.slaveNodes[idx] for idx in ctct.proj[:,3].nonzero()[0]]
        idxNodes = ctct.proj[:,3].nonzero()[0]
        activeNodes = [ctct.slaveNodes[idx] for idx in idxNodes if ctct.actives[idx] is not None ]
        # set_trace()
        activePatches = [pr[0] for pr in ctct.proj if pr[3]!=0]
        self.plot(ax, ref = 4,specialPatches = [(0,1,0,0.75),activePatches],
                              almostSpecialPatches = [(1.0,0.64,0.0,0.75), AlmostImportant],
                              specialNodes = {"red": activeNodes,}, 
                              time = t,
                              fintC = fintC,
                              plotHooks = Hooks,
                              OnlyMasterSurf=True
                              )

        #Camera effects
        tt = (t-t0)/(tf-t0)
        ax.view_init(elev=elevation[0]*(1-tt)+elevation[1]*tt,azim=azimut[0]*(1-tt)+azimut[1]*tt)
        ax.dist = distance[0]*(1-tt)+distance[1]*tt     #default dist = 10
        ax.set_proj_type('ortho')
        plot_coords(ax,orig=(minx-0.5+1.0,miny-0.5,minz-0.5))

        if iteration==None:
            plt.savefig(self.output_dir+"plots/fig"+str(increment)+".png",dpi = dpi)
        else:
            plt.savefig(self.output_dir+"plots/fig"+str(increment)+"-"+str(iteration)+".png",dpi = dpi)
        plt.close()



    def plot4D(self,UU,X,C, ax=None, undef=False):

        sed = np.array(self.SED)
        smin,smax = [min(sed),max(sed)]
        sed = (sed-smin)/(smax-smin)            # normalizing everything between 0 and 1

        if ax is None:
            import matplotlib.pyplot as plt
            from matplotlib.widgets import Slider, Button

            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')

        init_time = 0




        plotSphere(Circ,ax=ax)
        u = UU[0]
        plotTruss(X+u.reshape(-1,3),C,ax=ax,bar_types=btypes,show_nodes=True,showtypes=0)
        
        axtime = fig.add_axes([0.25, 0.1, 0.65, 0.03])

        time_slider = Slider(
            ax=axtime,
            label='time',
            valmin=0.0,
            valmax=0.99,
            valstep=0.01,
            valinit=init_time,
        )


        # The function to be called anytime a slider's value changes
        def update(val):
            xl = ax.get_xlim()
            yl = ax.get_ylim()
            zl = ax.get_zlim()
            az, el = ax.azim, ax.elev
            ax.clear()
            idx = round(100*time_slider.val)
            plotSphere(Circ,ax=ax)
            acts = np.asarray(ACTIVES[idx]).nonzero()[0]
            actives = [slaves[i] for i in acts]
            plotTruss(X+UU[idx].reshape(-1,3),C,ax=ax,bar_types=btypes,show_nodes=True,showtypes=0,actives=actives)
            ax.set_xlim(xl)
            ax.set_ylim(yl)
            ax.set_zlim(zl)
            ax.azim = az
            ax.elev = el

        # register the update function with each slider
        time_slider.on_changed(update)

        # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
        resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', hovercolor='0.975')


        def reset(event):
            time_slider.reset()
        button.on_clicked(reset)

        plt.show()





    def savedata(self,time,*args, dofs = "all", name = None, Sum=False):
        for arg in args:
            # if type(arg)==tuple:
            #     set_trace()
            filename = self.output_dir+"data_"+(arg if name is None else name)+".csv"
    
            # writing to csv file
            with open(filename, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                # creating a csv writer object
                csvwriter = csv.writer(csvfile)
                
                # writing the fields
                if type(dofs)==list or type(dofs)==np.ndarray:
                    if Sum:
                        csvwriter.writerow([time,sum(getattr(self,arg)[dofs])])
                    else:
                        csvwriter.writerow([time]+(getattr(self,arg)[dofs]).tolist())
                elif type(dofs)==str:
                    attr = getattr(self,arg)
                    if dofs == "all": 
                        if hasattr(attr,'__iter__'):
                            csvwriter.writerow([time]+attr.tolist())
                        else:
                            csvwriter.writerow([time]+[attr])
                    else: print("Unrecognized datatype to be stored from string. Not saving this data"), set_trace()
                else:
                    print("Unrecognized datatype to be stored. Not saving this data"), set_trace()
    
    def saveNodalData(self,body_id,time,*args, nodes = "all", name = None, Sum=False):
    
        if type(nodes)==list or type(nodes)==np.ndarray:
            dofs = self.bodies[body_id].DoFs[nodes]
        elif nodes == "all":
            dofs = self.bodies[body_id].DoFs
        for arg in args:
            filename = self.output_dir+(("sum_" if Sum else "")+(arg if name is None else name))+".csv"
    
            # writing to csv file
            with open(filename, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                # creating a csv writer object
                csvwriter = csv.writer(csvfile)
                # writing the fields
                if Sum:
                    sumX = sum(getattr(self,arg)[dofs.T[0]])
                    sumY = sum(getattr(self,arg)[dofs.T[1]])
                    sumZ = sum(getattr(self,arg)[dofs.T[2]])
                    csvwriter.writerow([time,sumX,sumY,sumZ])
                else:
                    csvwriter.writerow([time]+(getattr(self,arg)[dofs]).tolist())

    def saveContactData(self,ctct_id,time, slaveNodes = None):

        ctct = self.contacts[ctct_id]
        stringthings = ["gn", "kn"]
        slaveNodes = range(ctct.nsn) if slaveNodes is None else slaveNodes

        def states(idx):
            return [ctct.proj[idx,3], ctct.alpha_p[idx]*ctct.kn, ctct.hook[idx,3]]
        

        for stidx, stringthing in enumerate(stringthings):
            filename = self.output_dir+"ctct"+str(ctct_id)+"_"+stringthing+".csv"
    
            # writing to csv file
            with open(filename, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                # creating a csv writer object
                csvwriter = csv.writer(csvfile)
                
                # writing the fields
                csvwriter.writerow([time]+[states(idx)[stidx] for idx in slaveNodes])

    def applyBCs(self,t,tpre):

        print("Applying BCs from: ",round(tpre,6)," to: ",round(t,6))
        diri = []
        for bc in self.BCs:
            if bc.t0 <= t <= bc.tf or bc.t0 <= tpre <= bc.tf:
                # load factor
                if tpre < bc.t0:          # just entered in time interval
                    LF = (t - bc.t0) / (bc.tf - bc.t0)
                elif t > bc.tf:            # just exited time interval
                    LF = (bc.tf - tpre) / (bc.tf - bc.t0)
                else:                   # going through time interval (normal case)
                    LF = (t - tpre ) / (bc.tf - bc.t0)

                if bc.BCtype == "dirichlet":          # Dirichlet
                    diri.extend(bc.DoFs)
                    if LF!=0:
                        for dof in bc.DoFs:
                            self.u[dof] += bc.vals[dof%self.dim]*LF
                elif LF!=0:                           # Neumann
                    for dof in bc.DoFs:
                        self.fext[dof] += bc.vals[dof%self.dim]*LF

        di = list(set(diri))      #Gets rid of eventually repeated DoFs
        self.diri = di
        self.free = np.delete(np.arange(len(self.X)), di)       # free DOFs 'b'

    def IsContactAssumptionCorrect(self):
        set_trace()
        for contact in self.contacts:

            if len(contact.activePairs)!=len(contact.GNs):
                return False

            for pair, gn in zip(contact.candidatePairs,contact.GNs):
                if ((gn<0) != (pair in contact.activePairs)): return False
        return True

    def setReferences(self, justOutput = False):
        u_ref = self.u.copy()
        u_temp_ref = self.u_temp.copy()
        fint_ref = self.fint.copy()
        alpha_ref = list([contact.alpha_p.copy() for contact in self.contacts])
        new_ref   = list([contact.new.copy() for contact in self.contacts])
        xs_ref = list([contact.xs.copy() for contact in self.contacts])
        act_ref = list([contact.actives.copy() for contact in self.contacts])
        cand_ref = list([contact.candids.copy() for contact in self.contacts])

        FPconv_ref =  list([body.FPconv.copy() if not body.isRigid else None for body in self.bodies])
        EPcum_ref =  list([(body.EPcum).copy() if not body.isRigid else None for body in self.bodies])

        REF = [u_ref,u_temp_ref, fint_ref, alpha_ref, new_ref, xs_ref,act_ref,cand_ref,FPconv_ref,EPcum_ref]

        if justOutput:
            return REF

        self.REF = REF

    def getReferences(self, alphas=True,actives=False):
        self.u = self.REF[0].copy()
        self.u_temp = self.REF[1].copy()
        self.fint = self.REF[2].copy()
        alpha_ref, new_ref, xs_ref,act_ref,cand_ref,FPconv_ref,EPcum_ref = self.REF[3:]
        # self.u_temp = np.array(self.u)
        for icont, contact in enumerate(self.contacts):
            if alphas: contact.alpha_p = alpha_ref[icont]
            contact.new  = new_ref[icont].copy()
            contact.xs   = xs_ref[icont].copy()
            if actives: 
                contact.actives = act_ref[icont].copy()
                contact.candids = cand_ref[icont].copy()

        for i_b,body in enumerate(self.bodies):
            if not body.isRigid:
                body.FPconv = FPconv_ref[i_b].copy()
                body.FPtemp = FPconv_ref[i_b].copy()
                body.EPcum = EPcum_ref[i_b].copy()


    def get_fint(self, DispTime = False, temp = False):
        t0_fint = time.time()
        self.fint = np.zeros_like(self.u,dtype=float)
        for body in self.bodies: body.get_fint(self, temp = temp)
        if DispTime: print("Getting fint : ",time.time()-t0_fint,"s")

        u = self.u if temp==False else self.u_temp
        for ctct in self.contacts:
            if temp: ctct.patch_changes = []                 # to keep track of the changes during iterations
            sDoFs  = ctct.slaveBody.DoFs[ctct.slaveNodes]
            ctct.xs = np.array(ctct.slaveBody.X )[ctct.slaveNodes ] + np.array(u[sDoFs ])      # u_temp
            ctct.getfintC(self, DispTime=DispTime,useANN=False)      # uses xs
            # ctct.getfintC_unilateral(self, DispTime=DispTime,useANN=False)      # uses xs


    def get_K(self, DispTime = False):
        t0_K = time.time()
        self.K = sparse.coo_matrix((self.ndof,self.ndof),dtype=float)
        for body in self.bodies:
            body.get_K(self)    # uses model.u_temp

        printif(DispTime,"Getting K : ",time.time()-t0_K,"s")

        for contact in self.contacts:
            contact.getKC(self, DispTime=DispTime)     #uses model.u_temp
            # contact.getKC_unilateral(self, DispTime=DispTime)     #uses model.u_temp


    def solve_NR(self,tol=1e-10,maxiter=10,plotIters=False):

        if self.transform_2d is not None:
            return self.solve_NR_2d(tol=tol,maxiter=maxiter,plotIters=plotIters)

        # NEWTON-RAPHSON
        di, fr = self.diri, self.free
        TimeDisp = False
        RES, niter = 1+tol, 0

        self.get_fint(DispTime=TimeDisp,temp=True)   # uses self.u_temp     (no dirichlet)

        dua = self.u[di] - self.u_temp[di]      # difference due to dirichlet BCs applied
        
        RES_prev = 1e100
        RES = 1e99
        minRES = float(RES)

        # while RES>tol and RES<5*RES_prev and niter<maxiter and not np.isnan(RES):
        while RES>tol and niter<maxiter and not np.isnan(RES):
            RES_prev = RES

            # getting K and KC
            self.get_K(DispTime=TimeDisp)  # <-- uses self.u_temp           (no dirichlet)

            # Linear System

            Kaa=self.K[np.ix_(di,di)]
            Kab=self.K[np.ix_(di,fr)]
            Kba=self.K[np.ix_(fr,di)]
            Kbb=self.K[np.ix_(fr,fr)]

            fina=self.fint[di]
            finb=self.fint[fr]
            fexb=self.fext[fr]      # always zero (?)

            dub =spsolve(Kbb,fexb-finb-Kba.dot(dua))     # Uses all available processors
            self.fext[di] = fina + Kaa.dot(dua) + Kab.dot(dub)

            dua = np.zeros_like(self.u[di])

            self.u[fr] += dub

            self.get_fint(DispTime=TimeDisp)   # uses self.u

            RES = self.residual(printRes=True)

            if RES<minRES:
                minRES = float(RES)
            
            for ic, ctct in enumerate(self.contacts):
                pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
                with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerow([self.t if niter==0 else None,niter,RES]+ctct.patch_changes)

            if plotIters: self.plotContactNow()

            niter += 1
            self.u_temp = np.array(self.u)

            # # if (niter == maxiter and RES>tol) or np.isnan(RES) or RES>1e+15:
            # if (niter == maxiter and RES>tol) or np.isnan(RES) or RES>1e+15:
            #     print('Increment failed to converge !!!!! Redoing half')
            #     return False

        return RES<tol, RES

    def solve_NR_2d(self,tol=1e-10,maxiter=10,plotIters=True):

        # N = self.transform_2d*np.sqrt(0.5)
        Ns,Nt = self.transform_2d

        # Give (S)ymmetry to u, u_temp
        self.u = Ns@self.u
        self.u_temp = Ns@self.u_temp


        # NEWTON-RAPHSON
        bool_di = np.zeros(self.ndof)

        bool_di[self.diri] = 1
        di = np.where(Nt.T@bool_di!=0)[0]
        fr = list(set(range(Nt.shape[1]))-set(di))

        TimeDisp = False
        RES, niter = 1+tol, 0


        # set_trace()

        self.get_fint(DispTime=TimeDisp,temp=True)   # uses self.u_temp     (no dirichlet)
        fint_r = Nt.T@Ns.T@self.fint
        fext_r = np.zeros_like(fint_r)
        u_r = Nt.T@self.u
        u_r_temp = Nt.T@self.u_temp
        
        RES_prev = 1e100
        RES = 1e99
        minRES = float(RES)

        Nta = Nt[:,di]
        Ntb = Nt[:,fr]
        # dua = (Nt.T@(self.u - self.u_temp))[di]      # difference due to dirichlet BCs applied
        dua = (u_r - u_r_temp)[di]      # difference due to dirichlet BCs applied

        # set_trace() 
        # while RES>tol and RES<5*RES_prev and niter<maxiter and not np.isnan(RES):
        while RES>tol and  niter<maxiter and not np.isnan(RES):
            RES_prev = RES

            # getting K and KC
            self.get_K(DispTime=TimeDisp)  # <-- uses self.u_temp           (no dirichlet)
            Kr = Nt.T@Ns.T@self.K@Ns@Nt

            # Linear System
            Kaa=Kr[np.ix_(di,di)]
            Kab=Kr[np.ix_(di,fr)]
            Kba=Kr[np.ix_(fr,di)]
            Kbb=Kr[np.ix_(fr,fr)]

            fina=fint_r[di]
            finb=fint_r[fr]
            fexb=fext_r[fr]     # always zero (?)


            dub =spsolve(Kbb,fexb-finb-Kba.dot(dua))     # Uses all available processors
            fext_r[di] = (fina + Kaa.dot(dua) + Kab.dot(dub))
            # self.fext = Ns@self.fext


            dua[:] = 0.0

            u_r[fr] += dub

            self.u[self.free] = (Ns@Nt@u_r)[self.free]
            self.u = Ns@self.u

            # self.plotNow(as2D=True,OnlyMasterSurf=True)

            self.get_fint(DispTime=TimeDisp)   # uses self.u
            fint_r = Nt.T@Ns.T@self.fint


            # RES = self.residual(printRes=True)
            RES = norm(fint_r-fext_r)
            print("resid :",RES)


            if RES<minRES:
                minRES = float(RES)
            
            for ic, ctct in enumerate(self.contacts):
                pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
                with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerow([self.t if niter==0 else None,niter,RES]+ctct.patch_changes)

            if plotIters: self.plotContactNow()

            niter += 1
            self.u_temp = np.array(self.u)

            # # if (niter == maxiter and RES>tol) or np.isnan(RES) or RES>1e+15:
            # if (niter == maxiter and RES>tol) or np.isnan(RES) or RES>1e+15:
            #     print('Increment failed to converge !!!!! Redoing half')
            #     return False

        return RES<tol, RES

    def BFGS(self,FUN, JAC, u0, tol = 1e-8, tol2 = 1e-10):
        alpha_init = 1
        c_par = 1e-4
        c_par2 = 0.9
        r_par = 0.5
        
        free_ind = self.free
        nfr = len(free_ind)
        f = JAC(u0)
        m_new = FUN(u0)
        self.write_m_and_f(m_new,norm(f),0)
        u = u0.copy()
        f_2 = f[free_ind]
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)
        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        while np.linalg.norm(f_2 - f_new) > tol and np.linalg.norm(f_2) > 0:
            printalphas = True
            iter += 1
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old
            
            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)
            
            alpha3 = alpha_init
            ux = u.copy()
            ux[free_ind] = u[free_ind] + alpha3 * h_new
            f = JAC(ux)
            m_3 = FUN(ux)
            f_3 = f[free_ind]
            self.write_m_and_f(m_3,norm(f_3),iter)
            
            signal1 = 0
            if m_3 <= m_new + c_par * alpha3 * np.dot(h_new, f_new) and abs(np.dot(h_new, f_3)) <= c_par2 * abs(np.dot(h_new, f_new)):
                signal1 = 1
            
            iter2a = 0
            while m_3 < m_new + c_par * alpha3 * np.dot(h_new, f_new) and np.dot(h_new, f_3) < -c_par2 * np.dot(h_new, f_new) and signal1 == 0:
                alpha3 = alpha3 / r_par
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha3 * h_new
                f = JAC(ux)
                m_3 = FUN(ux)
                self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a)
                f_3 = f[free_ind]
                iter2a += 1
                if m_3 <= m_new + c_par * alpha3 * np.dot(h_new, f_new) and abs(np.dot(h_new, f_3)) <= c_par2 * abs(np.dot(h_new, f_new)):
                    signal1 = 1
            printif(printalphas,"alpha3:",alpha3)
            if signal1 == 0:
                alpha1 = 0
                alpha2 = alpha3 / 2
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha2 * h_new
                f = JAC(ux)
                m_2 = FUN(ux)
                self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a)
                f_2 = f[free_ind]
                signal2 = 0
                iter2 = 0
                while signal2 == 0:
                    iter2 += 1
                    if alpha3 - alpha1 < tol2:
                        signal2 = 1
                        m_2 = m_new
                        f_2 = f_new
                    elif m_2 > m_new + c_par * alpha2 * np.dot(h_new, f_new):
                        alpha3 = alpha2
                        m_3 = m_2
                        f_3 = f_2
                        alpha2 = 0.5 * (alpha1 + alpha2)
                        ux = u.copy()
                        ux[free_ind] = u[free_ind] + alpha2 * h_new
                        f = JAC(ux)
                        m_2 = FUN(ux)
                        self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a,iter2=iter2)
                        f_2 = f[free_ind]
                    else:
                        if np.dot(h_new, f_2) < c_par2 * np.dot(h_new, f_new):
                            alpha1 = alpha2
                            alpha2 = 0.5 * (alpha2 + alpha3)
                            ux = u.copy()
                            ux[free_ind] = u[free_ind] + alpha2 * h_new
                            f = JAC(ux)
                            m_2 = FUN(ux)
                            self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a,iter2=iter2)
                            f_2 = f[free_ind]
                        elif np.dot(h_new, f_2) > -c_par2 * np.dot(h_new, f_new):
                            alpha3 = alpha2
                            m_3 = m_2
                            f_3 = f_2
                            alpha2 = 0.5 * (alpha1 + alpha2)
                            ux = u.copy()
                            ux[free_ind] = u[free_ind] + alpha2 * h_new
                            f = JAC(ux)
                            m_2 = FUN(ux)
                            self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a,iter2=iter2)
                            f_2 = f[free_ind]
                        else:
                            signal2 = 1
                printif(printalphas,"alpha3:",alpha3, "  (after signal2=0)")
                            
            delta_u = ux[free_ind] - u[free_ind]
            u = ux
            if signal1 == 1:
                m_new = m_3
                f_2 = f_3.copy()

        return u, m_new, iter


    def BFGS_plastic_backup(self,FUNJAC, u0, tol = 1e-10, tol2 = 1e-15):
        # alpha_init = 0.001
        alpha_init = 1
        c_par = 1e-4
        c_par2 = 0.9
        r_par = 0.5
        
        free_ind = self.free
        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)

        u = u0.copy()
        f_2 = f[free_ind]
        m_new = norm(f_2)
        self.write_m_and_f(m_new,norm(f),0)
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        alpha = alpha_init
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > 1e-10:

            printalphas = False
            iter += 1
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old

            if iter>1:
                if not np.isfinite(norm(h_new)):
                    set_trace()


            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)

           
            if not np.isfinite(norm(h_new)):
                set_trace()

            m_new = abs(h_new@f_new)
            # m_new = norm(f_new)




            #TODO: plot values of m. Plot Armijo rule. And identify where Armijo is fulfilled

            n_points = 100
            ALPHAS = np.linspace(0.0, 100, n_points)
            MMs = np.zeros(n_points)
            HF = np.zeros(n_points)
            FF = np.zeros(n_points)
            ARM = np.zeros(n_points)
            CURV = np.zeros(n_points)
            
            for ai,alphai in enumerate(ALPHAS):
                ui = u.copy()
                ui[free_ind] = u[free_ind] + alphai * h_new
                print("BEFORE:\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp),"\talpha_i:",alphai)
                f_full , mi = FUNJAC(ui)
                print("AFTER :\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp))

                fi = f_full[free_ind]

                MMs[ai] = mi
                HF[ai] = h_new@fi
                FF[ai] = norm(fi)
                ARM[ai] = m_new + c_par*alphai*np.dot(h_new, f_new)
                # CURV[ai] = np.dot(h_new, fi) >= c_par2 * np.dot(h_new, f_new)
                
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)

            ax2.plot(ALPHAS,MMs,color = "blue",label="m")
            ax.plot(ALPHAS,HF,color = "green",label="h'f")
            # ax.plot(ALPHAS,FF,color = "yellow",label="|f|")
            # ax.plot(ALPHAS,ARM,color = "black",label="Armijo")
            plt.legend()
            plt.show()





            #alpha3 = alpha_init
            alpha3 = alpha
            ux = u.copy()
            ux[free_ind] = u[free_ind] + alpha3 * h_new

            f , m_3 = FUNJAC(ux)
            
            f_3 = f[free_ind]
            m_3 = abs(h_new@f_3)
            # m_3 = norm(f_3)
            self.write_m_and_f(m_3,norm(f_3),iter)
            
            signal1 = 0

            iter2a = 0
            if m_3 <= m_new + c_par * alpha3 * np.dot(h_new, f_new) and np.dot(h_new, f_3) >= c_par2 * np.dot(h_new, f_new):
                signal1 = 1

            elif  not (m_3 <= m_new + c_par * alpha3 * np.dot(h_new, f_new)) or m_3 is np.nan:
                pass
            else:
                while m_3 < m_new + c_par * alpha3 * np.dot(h_new, f_new) and signal1 == 0 and  m_3 is not np.nan:
                    iter2a += 1
                    alpha3 = alpha3 / r_par
                    ux = u.copy()
                    ux[free_ind] = u[free_ind] + alpha3 * h_new

                    f , m_3 = FUNJAC(ux)
                    
                    f_3 = f[free_ind]
                    m_3 = abs(h_new@f_3)
                    # m_3 = norm(f_3)
                    self.write_m_and_f(m_3,norm(f_3),iter,iter2a=iter2a)

                    if m_3 <= m_new + c_par * alpha3 * np.dot(h_new, f_new) and np.dot(h_new, f_3) >= c_par2 * np.dot(h_new, f_new):
                        signal1 = 1

                
            
            printif(printalphas,"alpha3:",alpha3)
            if signal1 == 0:
                alpha1 = 0

                alpha2 = alpha3 / 2
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha2 * h_new

                f , m_2 = FUNJAC(ux)
                
                f_2 = f[free_ind]
                m_2 = abs(h_new@f_2)
                # m_2 = norm(f_2)
                self.write_m_and_f(m_2,norm(f_2),iter,iter2a=iter2a)
                signal2 = 0
                iter2 = 0
                while signal2 == 0:
                    iter2 += 1
                    if alpha3 - alpha1 < tol2:
                        signal2 = 1
                        m_2 = m_new
                        f_2 = f_new
                    elif m_2 > m_new + c_par * alpha2 * np.dot(h_new, f_new) or m_2 is np.nan:
                        alpha3 = alpha2
                        alpha2 = 0.5 * (alpha1 + alpha2)
                        m_3 = m_2
                        f_3 = f_2
                        ux = u.copy()
                        ux[free_ind] = u[free_ind] + alpha2 * h_new

                        f , m_2 = FUNJAC(ux)
                        
                        f_2 = f[free_ind]
                        m_2 = abs(h_new@f_2)
                        # m_2 = norm(f_2)
                    
                        self.write_m_and_f(m_2,norm(f_2),iter,iter2a=iter2a,iter2=iter2,case = "a")

                    elif np.dot(h_new, f_2) < c_par2 * np.dot(h_new, f_new):
                        alpha1 = alpha2
                        alpha2 = 0.5 * (alpha2 + alpha3)

                        ux = u.copy()
                        ux[free_ind] = u[free_ind] + alpha2 * h_new
        
                        f , m_2 = FUNJAC(ux)
        
                        f_2 = f[free_ind]
                        m_2 = abs(h_new@f_2)
                        # m_2 = norm(f_2)
                        self.write_m_and_f(m_2,norm(f_2),iter,iter2a=iter2a,iter2=iter2,case = "b")
                    else:
                        signal2 = 1

                printif(printalphas,"alpha3:",alpha3, "  (after signal2=0)")
                            
            delta_u = ux[free_ind] - u[free_ind]
            u = ux
            if signal1 == 1:
                m_new = m_3
                f_2 = f_3.copy()
                alpha = alpha3
            else:
                alpha = alpha2

            print("|f2|:",norm(f_2))

        return u, m_new, iter

    def BFGS_plastic(self,FUNJAC, u0, tol = 1e-10, tol2 = 1e-15):
        # alpha_init = 0.001
        alpha_init = 1
        c_par = 1e-4
        c_par2 = 0.9
        r_par = 0.5
        
        free_ind = self.free
        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)

        u = u0.copy()
        f_2 = f[free_ind]
        m_new = norm(f_2)
        self.write_m_and_f(m_new,norm(f),0)
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        alpha = alpha_init
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > 1e-10:

            printalphas = False
            iter += 1
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old

            if iter>1:
                if not np.isfinite(norm(h_new)):
                    set_trace()


            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)

           
            if not np.isfinite(norm(h_new)):
                set_trace()

            m_new = abs(h_new@f_new)
            # m_new = norm(f_new)


            # Plotting lineseach...
            if iter>71:

                n_points = 200
                ALPHAS = np.linspace(0.0, 0.5*alpha_init, n_points)
                MMs = np.zeros(n_points)
                HF = np.zeros(n_points)
                FF = np.zeros(n_points)
                ARM = np.zeros(n_points)
                CURV = np.zeros(n_points)
                
                for ai,alphai in enumerate(ALPHAS):
                    ui = u.copy()
                    ui[free_ind] = u[free_ind] + alphai * h_new


                    # print("BEFORE:\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp),"\talpha_i:",alphai)
                    f_full , mi = FUNJAC(ui,show=True)
                    # print("AFTER :\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp))

                    fi = f_full[free_ind]

                    MMs[ai] = mi


                    HF[ai] = h_new@fi
                    FF[ai] = norm(fi)
                    ARM[ai] = m_new + c_par*alphai*np.dot(h_new, f_new)
                    # CURV[ai] = np.dot(h_new, fi) >= c_par2 * np.dot(h_new, f_new)
                    
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(111)

                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)

                ax.plot(ALPHAS,HF,color = "green",label="h'f")
                ax.plot(ALPHAS,FF,color = "yellow",label="|f|")
                ax2.plot(ALPHAS,MMs,color = "blue",label="m")
                ax2.plot(ALPHAS,ARM,color = "black",label="Armijo")
                plt.legend()
                plt.show()




            alpha3 = 2*alpha_init
            m_3=np.nan
            while m_3 is np.nan:
                alpha3 /= 2
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha3 * h_new

                f , m_3 = FUNJAC(ux)
                
                f_3 = f[free_ind]
                m_3 = h_new@f_3

            a = h_new.T@(f_3-f_new)/alpha3
            alpha = -(h_new.T@f_new)/a

            


            if iter==71:
                # alpha = 715
                alpha = 701.6221041531118
                print("iter 71")


            # if a<0:
            #     dalphax = alpha3
            #     alphax = alpha3
            #     f_x = f_3
            #     k_x = a
            #     kx_prev = k_x*10

            #     # while (abs(k_x)<abs(kx_prev)) and k_x*kx_prev>0 and not (np.dot(h_new, f_x) >= c_par2*np.dot(h_new,f_new)):
            #     while not (np.dot(h_new, f_x) >= c_par2*np.dot(h_new,f_new)):
            #         alphax += dalphax
            #         ux = u.copy()
            #         fx_prev = f_x
            #         kx_prev = k_x

            #         ux[free_ind] = u[free_ind] + alphax * h_new
            #         f , mx = FUNJAC(ux)
            #         f_x = f[free_ind]
            #         k_x = h_new.T@(f_x-fx_prev)/dalphax

            #         print("alphax:",alphax, "\tmx:",mx,"\th'fx:",h_new.T@f_x,"\tkx:",k_x)

            #         if (h_new.T@f_x)*(h_new.T@f_3)<0:       # if crossed the line...
            #             dalphax = -abs(dalphax)/2           # go back a bit
            #             k_x = kx_prev
            #             f_x = fx_prev
            #         else:
            #             dalphax = abs(k_x*(dalphax)/(k_x-kx_prev))

            #     alpha = alphax

            #     set_trace()



            ux = u.copy()
            ux[free_ind] = u[free_ind] + alpha * h_new
            f , m_3 = FUNJAC(ux)
            f_3 = f[free_ind]
            m_3 = h_new@f_3

            if np.dot(h_new, f_3) >= c_par2 * np.dot(h_new, f_new):
                f_2 = f_3.copy()
                m_2 = m_3
            else:
                f_2 = f_new.copy()
                m_2 = m_new


            # m_3 = norm(f_3)
            self.write_m_and_f(m_2,norm(f_2),iter)
                       
            delta_u = ux[free_ind] - u[free_ind]
            u = ux
            

            print("alpha:",alpha,"\t|f2|:",norm(f_2))

        return u, m_new, iter

    def BFGS_plastic_diego(self,FUNJAC, u0, tol = 1e-10, tol2 = 1e-15, free_ind = None):
        # alpha_init = 0.01
        alpha_init = 1
        # c_par2 = 0.9
        c_par2 = 0.95

        if free_ind is None:
            free_ind = self.free

        
        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)

        u = u0.copy()
        f_2 = f[free_ind]
        m_new = norm(f_2)
        self.write_m_and_f(m_new,norm(f),0)
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        alpha = alpha_init
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > tol:

            iter += 1
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old

            if iter>1:
                if not np.isfinite(norm(h_new)):
                    set_trace()


            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)

           
            if not np.isfinite(norm(h_new)):
                set_trace()

            m_new = abs(h_new@f_new)
            # m_new = norm(f_new)


            # Plotting lineseach...
            if iter<-1:

                n_points = 50
                ALPHAS = np.linspace(0.0, 10*alpha_init, n_points)
                MMs = np.zeros(n_points)
                HF = np.zeros(n_points)
                FF = np.zeros(n_points)
                ARM = np.zeros(n_points)
                CURV = np.zeros(n_points)
                
                for ai,alphai in enumerate(ALPHAS):
                    ui = u.copy()
                    ui[free_ind] = u[free_ind] + alphai * h_new
                    # print("BEFORE:\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp),"\talpha_i:",alphai)
                    f_full , mi = FUNJAC(ui)
                    # print("AFTER :\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp))

                    fi = f_full[free_ind]

                    MMs[ai] = mi
                    HF[ai] = h_new@fi
                    FF[ai] = norm(fi)
                    # ARM[ai] = m_new + c_par*alphai*np.dot(h_new, f_new)
                    # CURV[ai] = np.dot(h_new, fi) >= c_par2 * np.dot(h_new, f_new)
                    
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(111)

                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)

                ax2.plot(ALPHAS,MMs,color = "blue",label="m")
                ax.plot(ALPHAS,HF,color = "green",label="h'f")
                # ax.plot(ALPHAS,FF,color = "yellow",label="|f|")
                # ax.plot(ALPHAS,ARM,color = "black",label="Armijo")
                plt.legend()
                plt.show()

            a1 = 0
            f1 = h_new@f_new

            # compute f3 ensuring that it is a finite number
            f3 = np.nan
            a3 = 2*alpha_init
            while not np.isfinite(norm(f3)):
                a3 /= 2
                ux = u.copy()
                ux[free_ind] = u[free_ind] + a3*h_new
                f , _ = FUNJAC(ux)
                f_3 = f[free_ind]
                f3 = h_new@f_3


            # TODO: currently 3rd point leads to parabilic interpolation. Try to see if f1<f3<0, and the line formed
            # cuts the OX axis not too far, use a linear interpolation to compute the 3rd point intead of the average alpha 



            a2 = 0.5*(a1+a3)
            ux = u.copy()
            ux[free_ind] = u[free_ind] + a2*h_new
            f , _ = FUNJAC(ux)
            f_2 = f[free_ind]
            f2 = h_new@f_2

            # In case of concavity, push from the left until convex
            # while not f1<f3:
            while not f1<f2<f3 and max([abs(f1),abs(f2),abs(f3)])>1e-25:
                parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])

                if parab["a"] <0:
                    # print("\t...but parab.a <0")
                    break
                

                if parab['zeros'] is None:      # this is the begining of 'f' is concave itself. We move to the right
                    delta = a2-a1

                    a1 = a2
                    f1 = f2

                    a2 = a3 
                    f2 = f3

                    a3 += 1.5*delta
                    ux = u.copy()
                    ux[free_ind] = u[free_ind] + a3*h_new
                    f , _ = FUNJAC(ux)
                    f_3 = f[free_ind]
                    f3 = h_new@f_3

                    # print("\tleft: parab.zeros is None, all moved to the right:",[a1,a2,a3],"\f3:",f3)

                    continue

                alpha_to_min  = parab['minimum'][0]
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha_to_min*h_new
                f , _ = FUNJAC(ux)
                f_m = f[free_ind]
                fm = h_new@f_m

                print("\tFrom Left: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3],"\tam:",alpha_to_min, "\tf0:",fm)


                if alpha_to_min>a3:

                    a1 = a2
                    f1 = f2

                    a2 = a3
                    f2 = f3

                    a3 = alpha_to_min
                    f3 = fm

                elif alpha_to_min>a2:

                    a1 = a2
                    f1 = f2

                    a2 = alpha_to_min
                    f2 = fm

                else:
                    a1 = alpha_to_min
                    f1 = fm




            f0 = 100
            f_0 = f_3

            # print("\tbefore (right): \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])
            # while f0>0 or not (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new)):
            while f0>1e-12 or not (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new)):


                # To avoid cases like: [0.0, 8.0e-4, 8.1e-4] [-600000, 1.0e-6,1.1e-6]. Solution: put a2 more in the middle 
                slope_change = (f2-f1)/(f3-f2)
                while slope_change>1000 and max([abs(f1),abs(f2),abs(f3)])>1e-15:
                    a2=(a1+a2)/2
                    ux = u.copy()
                    ux[free_ind] = u[free_ind] + a2*h_new
                    f , _ = FUNJAC(ux)
                    f_2 = f[free_ind]
                    f2 = h_new@f_2
                    print("Too much slope change!!\ta2:",a2,"\tf2:",f2)

                    slope_change = (f3-f2)/(f2-f1)


                try:
                    parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])

                    if parab['zeros'] is None:      # this is the begining of 'f' is concave itself. We move to the right
                        delta = a2-a1

                        a1 = a2
                        f1 = f2

                        a2 = a3
                        f2 = f3

                        a3 += delta
                        ux = u.copy()
                        ux[free_ind] = u[free_ind] + a3*h_new
                        f , _ = FUNJAC(ux)
                        f_3 = f[free_ind]
                        f3 = h_new@f_3

                        continue
                    
                    if parab["a"]>0:
                        alpha_to_zero = max(parab['zeros'])
                    else:
                        alpha_to_zero = min(parab['zeros'])

                except:
                    alpha_to_zero = (f1*a3-f3*a1)/(f1-f3)

                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha_to_zero*h_new
                f , m0 = FUNJAC(ux)
                f_0 = f[free_ind]
                f0 = h_new@f_0

                print("\tFrom Right: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3],"\ta0:",alpha_to_zero, "\tf0:",f0)


                if alpha_to_zero>a3:    # In this case f1,f2,f3 are all negative. since f0=0 and f is increasing in the interval

                    a1 = a2
                    f1 = f2

                    a2 = a3
                    f2 = f3

                    a3 = alpha_to_zero
                    f3 = f0

                elif alpha_to_zero>a2:      # In this case f1,f2<0,  f3>0

                    # if f2<0:
                    a1 = a2
                    f1 = f2

                    a2 = alpha_to_zero
                    f2 = f0

                elif alpha_to_zero>a1:      # Here f2,f3>0 but f0 could be positive and in that case it should NOT replace f1
                    # if f0<0 and abs(a2-alpha_to_zero)/abs(a1-alpha_to_zero)<100:
                    if f0<0:
                        a1 = alpha_to_zero
                        f1 = f0
                    else:
                        a3 = a2
                        f3 = f2
                        
                        a2 = alpha_to_zero
                        f2 = f0

                else:
                    # It shouldn't even reach here because f1<0 and f increases
                    a3 = a2
                    f3 = f2
                    
                    a2 = a1
                    f2 = f1

                    a1 = alpha_to_zero
                    f1 = f0

                
            alpha = alpha_to_zero
            f_2 = f_0
            m_2 = f0

            print("\tFinally: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])


            # m_3 = norm(f_3)
            self.write_m_and_f(m_2,norm(f_2),iter)
                       
            delta_u = ux[free_ind] - u[free_ind]
            u = ux

            if alpha>580:
                set_trace()                

            print("alpha:",alpha,"\t|f2|:",norm(f_2))

        return u, m0, iter , norm(f_2)

    def BFGS_plastic_diego2(self,FUNJAC, u0, tol = 1e-10, tol2 = 1e-15, free_ind = None):
        # alpha_init = 0.01
        alpha_init = 1
        # c_par2 = 0.9
        c_par2 = 0.95

        if free_ind is None:
            free_ind = self.free

        
        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)

        u = u0.copy()
        f_2 = f[free_ind]
        m_new = norm(f_2)
        self.write_m_and_f(m_new,norm(f),0)
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        alpha = alpha_init
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > tol:

            iter += 1
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old

            if iter>1:
                if not np.isfinite(norm(h_new)):
                    set_trace()


            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)

           
            if not np.isfinite(norm(h_new)):
                set_trace()

            m_new = abs(h_new@f_new)
            # m_new = norm(f_new)


            # Plotting lineseach...
            if iter<-1:

                n_points = 50
                ALPHAS = np.linspace(0.0, 10*alpha_init, n_points)
                MMs = np.zeros(n_points)
                HF = np.zeros(n_points)
                FF = np.zeros(n_points)
                ARM = np.zeros(n_points)
                CURV = np.zeros(n_points)
                
                for ai,alphai in enumerate(ALPHAS):
                    ui = u.copy()
                    ui[free_ind] = u[free_ind] + alphai * h_new
                    # print("BEFORE:\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp),"\talpha_i:",alphai)
                    f_full , mi = FUNJAC(ui)
                    # print("AFTER :\tDEP_cum:", norm(self.bodies[0].DELTA_EPcum),"\tFPtemp:",norm(self.bodies[0].FPtemp))

                    fi = f_full[free_ind]

                    MMs[ai] = mi
                    HF[ai] = h_new@fi
                    FF[ai] = norm(fi)
                    # ARM[ai] = m_new + c_par*alphai*np.dot(h_new, f_new)
                    # CURV[ai] = np.dot(h_new, fi) >= c_par2 * np.dot(h_new, f_new)
                    
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(111)

                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)

                ax2.plot(ALPHAS,MMs,color = "blue",label="m")
                ax.plot(ALPHAS,HF,color = "green",label="h'f")
                # ax.plot(ALPHAS,FF,color = "yellow",label="|f|")
                # ax.plot(ALPHAS,ARM,color = "black",label="Armijo")
                plt.legend()
                plt.show()

            a1 = 0
            f1 = h_new@f_new

            # compute f3 ensuring that it is a finite number
            f3 = np.nan
            a3 = 2*alpha_init
            while not np.isfinite(norm(f3)):
                a3 /= 2
                ux = u.copy()
                ux[free_ind] = u[free_ind] + a3*h_new
                f , _ = FUNJAC(ux)
                f_3 = f[free_ind]
                f3 = h_new@f_3


            # TODO: currently 3rd point leads to parabilic interpolation. Try to see if f1<f3<0, and the line formed
            # cuts the OX axis not too far, use a linear interpolation to compute the 3rd point intead of the average alpha 



            a2 = 0.5*(a1+a3)
            ux = u.copy()
            ux[free_ind] = u[free_ind] + a2*h_new
            f , _ = FUNJAC(ux)
            f_2 = f[free_ind]
            f2 = h_new@f_2

            # In case of concavity (in 'm'), push from the left until convex
            # while not f1<f3:
            while not f1<f2<f3 and max([abs(f1),abs(f2),abs(f3)])>1e-25:
                parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])

                if parab["a"] <0:
                    # print("\t...but parab.a <0")
                    break
                
                if parab['zeros'] is None:      # this is the begining of 'f' is concave itself. We move to the right
                    delta = a2-a1

                    a1 = a2
                    f1 = f2

                    a2 = a3 
                    f2 = f3

                    a3 += 1.5*delta
                    ux = u.copy()
                    ux[free_ind] = u[free_ind] + a3*h_new
                    f , _ = FUNJAC(ux)
                    f_3 = f[free_ind]
                    f3 = h_new@f_3

                    # print("\tleft: parab.zeros is None, all moved to the right:",[a1,a2,a3],"\f3:",f3)
                    continue

                alpha_to_min  = parab['minimum'][0]
                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha_to_min*h_new
                f , _ = FUNJAC(ux)
                f_m = f[free_ind]
                fm = h_new@f_m

                print("\tFrom Left: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3],"\tam:",alpha_to_min, "\tf0:",fm)

                if alpha_to_min>a3:

                    a1 = a2
                    f1 = f2

                    a2 = a3
                    f2 = f3

                    a3 = alpha_to_min
                    f3 = fm

                elif alpha_to_min>a2:

                    a1 = a2
                    f1 = f2

                    a2 = alpha_to_min
                    f2 = fm

                else:
                    a1 = alpha_to_min
                    f1 = fm


            f0 = 100
            f_0 = f_3

            # print("\tbefore (right): \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])
            # while f0>0 or not (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new)):
            while f0>1e-12 or not (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new)):


                # To avoid cases like: [0.0, 8.0e-4, 8.1e-4] [-600000, 1.0e-6,1.1e-6]. Solution: put a2 more in the middle 
                slope_change = (f2-f1)/(f3-f2)
                while slope_change>1000 and max([abs(f1),abs(f2),abs(f3)])>1e-15:
                    a2=(a1+a2)/2
                    ux = u.copy()
                    ux[free_ind] = u[free_ind] + a2*h_new
                    f , _ = FUNJAC(ux)
                    f_2 = f[free_ind]
                    f2 = h_new@f_2
                    print("Too much slope change!!\ta2:",a2,"\tf2:",f2)

                    slope_change = (f3-f2)/(f2-f1)


                try:
                    parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])

                    if parab['zeros'] is None:      # this is the begining of 'f' is concave itself. We move to the right
                        delta = a2-a1

                        a1 = a2
                        f1 = f2

                        a2 = a3
                        f2 = f3

                        a3 += delta
                        ux = u.copy()
                        ux[free_ind] = u[free_ind] + a3*h_new
                        f , _ = FUNJAC(ux)
                        f_3 = f[free_ind]
                        f3 = h_new@f_3

                        continue
                    
                    if parab["a"]>0:
                        alpha_to_zero = max(parab['zeros'])
                    else:
                        alpha_to_zero = min(parab['zeros'])

                except:
                    alpha_to_zero = (f1*a3-f3*a1)/(f1-f3)

                ux = u.copy()
                ux[free_ind] = u[free_ind] + alpha_to_zero*h_new
                f , m0 = FUNJAC(ux)
                f_0 = f[free_ind]
                f0 = h_new@f_0

                print("\tFrom Right: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3],"\ta0:",alpha_to_zero, "\tf0:",f0)


                if alpha_to_zero>a3:    # In this case f1,f2,f3 are all negative. since f0=0 and f is increasing in the interval

                    a1 = a2
                    f1 = f2

                    a2 = a3
                    f2 = f3

                    a3 = alpha_to_zero
                    f3 = f0

                elif alpha_to_zero>a2:      # In this case f1,f2<0,  f3>0

                    if (a3-alpha_to_zero)/(alpha_to_zero-a2)>20:
                        a3 = alpha_to_zero
                        f3 = f0
                        
                    else:
                        a1 = a2
                        f1 = f2

                        a2 = alpha_to_zero
                        f2 = f0




                elif alpha_to_zero>a1:      # Here f2,f3>0 but f0 could be positive and in that case it should NOT replace f1
                    # # if f0<0 and abs(a2-alpha_to_zero)/abs(a1-alpha_to_zero)<100:
                    # if f0<0:
                    #     a1 = alpha_to_zero
                    #     f1 = f0
                    # else:
                    #     a3 = a2
                    #     f3 = f2
                        
                    #     a2 = alpha_to_zero
                    #     f2 = f0
                    if (a2-alpha_to_zero)/(alpha_to_zero-a1)>20 or f0>0:
                        a3 = a2
                        f3 = f2
                        
                        a2 = alpha_to_zero
                        f2 = f0
                    else:
                        a1 = alpha_to_zero
                        f1 = f0

                else:
                    # It shouldn't even reach here because f1<0 and f increases
                    a3 = a2
                    f3 = f2
                    
                    a2 = a1
                    f2 = f1

                    a1 = alpha_to_zero
                    f1 = f0

                
            alpha = alpha_to_zero
            f_2 = f_0
            m_2 = f0

            print("\tFinally: \talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])


            # m_3 = norm(f_3)
            self.write_m_and_f(m_2,norm(f_2),iter)
                       
            delta_u = ux[free_ind] - u[free_ind]
            u = ux

            # if alpha<1:
            #     set_trace()    



            print("alpha:",alpha,"\t|f2|:",norm(f_2))

        return u, m0, iter , norm(f_2)

    def TR_plastic(self,FUNJAC,HESS, u0, tol = 1e-10):

        free_ind = self.free

        TR_rad = 0.001
        TR_rad_min = 0
        TR_rad_max = 1

        FPtemp_accepted =  list([body.FPtemp.copy() if not body.isRigid else None for body in self.bodies])
        delta_EPcum_accepted =  list([(body.DELTA_EPcum).copy() if not body.isRigid else None for body in self.bodies])

        ff , m_new = FUNJAC(u0)
        f_new = ff[free_ind]

        for i_b,body in enumerate(self.bodies):
            if not body.isRigid:
                body.FPtemp = FPtemp_accepted[i_b].copy()
                body.DELTA_EPcum = delta_EPcum_accepted[i_b].copy()

        

        rho = 1
        iter = 0
        update_signal = 1

        h_full = np.zeros_like(u0)

        u = u0
        # set_trace()
        # while rho>0 and norm(f_new)>tol:
        # while rho!=0 and norm(f_new)>tol:
        while norm(f_new)>tol:

            iter += 1
            if update_signal==1:
                K = HESS(u)
                K = K[np.ix_(free_ind,free_ind)]

            h = (-(f_new.T@f_new)/float(f_new.T@K@f_new)) * f_new      # minimizer in the steepest decent direction


            if norm(h)>TR_rad:
                h = h/norm(h)*TR_rad
            else:
                if update_signal==1:
                    newt =np.linalg.solve(K,-f_new).ravel()
                if norm(newt)<TR_rad:
                    h = newt
                else:
                    cauchy = h.copy()
                    a = cauchy.T@(cauchy - 2*newt)+ newt.T@newt
                    b = 2*cauchy.T@(newt-cauchy)
                    c = cauchy.T@cauchy - TR_rad**2

                    tau = 0
                    eq = c

                    while abs(eq)>1e-12:
                        stiff = 2*a*tau + b
                        tau -= eq/stiff

                        eq = a*tau**2 + b*tau + c

                    h = cauchy + tau*(newt-cauchy)

            
            h_full[free_ind] = h
            ff , m_h = FUNJAC(u + h_full)
            f_h = ff[free_ind]
        
            rho = -(m_new-m_h)/(0.5*h.T@K@h + f_new.T@h)
            # rho = abs(-(m_new-m_h)/(0.5*h.T@K@h + f_new.T@h))

            if rho<0.25:
                TR_rad = max(0.25*TR_rad,TR_rad_min)
            else:
                
                if (rho>0.75) and (norm(h)<TR_rad+1e-5) and (norm(h)>TR_rad-1e-5):
                    TR_rad = min(2*TR_rad,TR_rad_max)
            

            update_signal = 0

            # if rho>0.25 or  (TR_rad == TR_rad_min and m_h<m_new):
            if m_h<m_new:
                u += h_full
                m_new = m_h
                f_new = f_h
                update_signal = 1

                FPtemp_accepted =  list([body.FPtemp.copy() if not body.isRigid else None for body in self.bodies])
                delta_EPcum_accepted =  list([(body.DELTA_EPcum).copy() if not body.isRigid else None for body in self.bodies])

            else:
                # TODO: reset prpojections and plastics
                for i_b,body in enumerate(self.bodies):
                    if not body.isRigid:
                        body.FPtemp = FPtemp_accepted[i_b].copy()
                        body.DELTA_EPcum = delta_EPcum_accepted[i_b].copy()

            print("h:",norm(h), "\tf:",norm(f_new),"\trho:",rho,"\tTR:",TR_rad, "\tm_h:",m_h, "\tdm:",m_h-m_new)

        return u, m_new, iter



    def write_m_and_f(self,m,f,iter,iter2a=0,iter2=0,case=0):
        for ic, ctct in enumerate(self.contacts):
            pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
            with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow([self.t if iter==0 else None,[iter,iter2a,iter2,case],[m,f]]+ctct.patch_changes)



    def Energy(self,u):
        En,EnC = 0.0, 0.0
        for body in self.bodies: 
            En+=body.compute_m_plastic(u)
            # En+=body.compute_m(u)
        for ctct in self.contacts:
            EnC+=ctct.compute_m(u)
        # print('En',En,"\tEnC",EnC,'\t\tEn_tot',En+EnC)
        return En+EnC
    
    def Force(self,u):
        force = np.zeros_like(self.fint)
        for body in self.bodies: 
            force+=body.compute_f_plastic(u,self)  #called like that to differentiate from get_fint which takes u=Model.u and modifies Model.fint
            # force+=body.compute_f(u,self)  #called like that to differentiate from get_fint which takes u=Model.u and modifies Model.fint
        for ctct in self.contacts:
            force+=ctct.compute_f(u,self)
        # print('Force',np.linalg.norm(force))
        force[self.diri]=0.0
        return force

    def Energy_and_Force(self,u,show=False):
        if self.transform_2d is not None:
            Ns,Nt = self.transform_2d
            u = (Ns@Nt@u)

        En,EnC = 0.0, 0.0
        force = np.zeros_like(self.fint)

        for body in self.bodies: 

            force_bi, Ebi = body.compute_mf_plastic(u,self)
            En += Ebi
            force+=force_bi

            # En+=body.compute_m(u)
        for ctct in self.contacts:
            # mCi,fCi = ctct.compute_mf(u,self)     # Bilateral constraint
            mCi,fCi = ctct.compute_mf_unilateral(u,self)
            EnC+=mCi
            force+=fCi

        if show:
            print("Eb:",En,"\tEc:",EnC)

        if self.transform_2d is not None:
            return Nt.T@Ns.T@force, En+EnC


        return force, En+EnC

    def Hessian(self,u):
        # K = sparse.coo_matrix((self.ndof,self.ndof),dtype=float)
        K = np.zeros((self.ndof,self.ndof),dtype=float)

        for body in self.bodies: 
            K += body.compute_k_plastic(u,self)    # uses model.u_temp

        for ctct in self.contacts:
            K += ctct.compute_k(u)     #uses model.u_temp
        # print('En',En,"\tEnC",EnC,'\t\tEn_tot',En+EnC)
        return K

        # t0_K = time.time()
        # for body in self.bodies:

        # printif(DispTime,"Getting K : ",time.time()-t0_K,"s")

        # for contact in self.contacts:
        #     contact.getKC(self, DispTime=DispTime)     #uses model.u_temp



    def solve_BFGS(self,tol=1e-10,maxiter=10,plotIters=False):
        u0=np.array(self.u)
        # MIN = minimize(self.Energy,u0,method='BFGS',jac=self.Force,options={'disp':False})
        # self.u = MIN.x
        self.u, m_new, iter = self.BFGS(self.Energy,self.Force,u0)
        return True

    def solve_BFGS_plastic(self,tol=1e-10,maxiter=10,plotIters=False):

        if self.transform_2d is not None:
            Ns,Nt = self.transform_2d
            u0 = Nt.T@Ns@np.array(self.u)                                 # This reduces the vector u0
            bool_di = np.zeros(self.ndof)
            bool_di[self.diri] = 1
            di = np.where(Nt.T@bool_di!=0)[0]
            fr = list(set(range(Nt.shape[1]))-set(di))      # free dofs in the reduced model

        else:
            u0=np.array(self.u)
            fr = None
    
        # MIN = minimize(self.Energy,u0,method='BFGS',jac=self.Force,options={'disp':False})
        # self.u = MIN.x
        # self.u, m_new, iter = self.BFGS_plastic(self.Energy_and_Force,u0)
        self.u, m_new, iter,res = self.BFGS_plastic_diego2(self.Energy_and_Force,u0,free_ind = fr)

        if self.transform_2d is not None:
            self.u = (Ns@Nt@self.u)


        return True,res




    def solve_TR_plastic(self,tol=1e-10,maxiter=10,plotIters=False):
        u0=np.array(self.u)
        self.u, m_new, iter = self.TR_plastic(self.Energy_and_Force,self.Hessian,u0)
        return True




    def residual(self, printRes = False):
        RES = np.linalg.norm(self.fint - self.fext)
        if printRes: print("resid :",RES)
        return RES


    def Solve(self, t0 = 0, tf = 1, TimeSteps = 1, max_iter=10 , recover = False, ForcedShift = False):
        self.ntstps = TimeSteps
        
        dt_base = (tf-t0)/TimeSteps; tolerance = 1e-10; ndof = len(self.X)
        t = t0 ; ti = 0
        # MaxBisect = 2**20   # Here we don't want to use minimization. Let's just see what NR can do
        MaxBisect = 2**10
        dt = dt_base; u_ref = np.array(self.u)
        num = 0
        den = 1
        self.bisect = int(np.log2(den))
        self.Redos_due_to_kn_adjustment = 0


        pickle.dump({body.name:[body.surf.quads,body.surf.nodes] for body in self.bodies},open(self.output_dir+"RecoveryOuters.dat","wb"))
        
        self.setReferences()
        skipBCs = False
        tmp = 0

        tracing = False

        if recover:
            self.REF,t, dt, ti,[num,den] = pickle.load(open("OUTPUT_202410141316pseudo2d/"+"RecoveryData.dat","rb"))
            # self.REF,t, dt, ti,[num,den] = pickle.load(open("OUTPUT_202410081833pseudo2d/"+"RecoveryData.dat","rb"))
            self.bisect = int(np.log2(den))
            
            self.getReferences(actives=True)
            for contact in self.contacts:   contact.getCandidates(self.u)
            # for contact in self.contacts:   contact.getCandidatesANN(self.u)

            if ForcedShift:
                rem = t%dt_base     # remaining time from previous time increment
                t_pre = t-rem
                dt = t_pre+dt_base-t
                ForcedShift = False

        DoMinimization = False

        redo_count = 0
        residuals_incr = []
        objectives_found = []
        u_found_incr = []

        while t+dt < tf+1e-4:
            if ForcedShift:
                rem = t%dt_base     # remaining time from previous time increment
                t_pre = t-rem
                dt = t_pre+dt_base-t
                ForcedShift = False
                num = den


            print(" ")
            t += dt
            num = min(num+1,den)

            self.t = t
            print("time: ",round(t,6),"s  -----------------------------------")
            print("("+str(num-1)+"/"+str(den)+" --> "+str(num)+"/"+str(den)+" in the increment)")

            for body in self.bodies:
                body.DELTA_EPcum = np.zeros((len(body.hexas),8))


            self.u_temp = np.array(self.u)  # copy of 'u' wihout dirichlet (to be used in Newton ONLY)

            if not skipBCs:
                self.applyBCs(t,t-dt)
            else:
                skipBCs = False

            print("du_diri:",norm(self.u[self.diri]-self.u_temp[self.diri]))

            ########################
            ### Solver Algorithm ###
            ########################
            actives_before_solving = list(self.contacts[0].actives)

            if DoMinimization:
                converged, res = self.solve_BFGS_plastic(tol=tolerance)
                # converged = self.solve_TR_plastic(tol=tolerance)
                self.u_temp = np.array(self.u)  # copy of 'u' so that solution is directly used in NR
            else:
                converged, res = self.solve_NR(tol=tolerance,maxiter=max_iter)
            # converged, res = self.solve_NR(tol=tolerance,maxiter=max_iter)
            actives_after_solving = list(self.contacts[0].actives)

            print("ACTIVE NODES:")
            print("Before solving:",actives_before_solving)
            print("After solving :",actives_after_solving)


            print("delta_epcum:",norm(self.bodies[0].DELTA_EPcum))


            if not converged:
                print('Increment failed to converge !!!!! Redoing half')
                t -= dt
                num -= 1

                self.getReferences(actives=True)    # resets u=u_temp <- u_ref

                pchfile = self.output_dir+"ctct"+str(0)+"iters_details.csv"
                with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                    csvwriter = csv.writer(csvfile)

                    # if int(round(dt_base/dt)) >= MaxBisect:
                    if den >= MaxBisect:
                        print("MAXIMUM BISECTION LEVEL (%s) REACHED:"%MaxBisect)
                        # set_trace()
                        DoMinimization=True
                        csvwriter.writerow(['','','MINIMIZATION'])
                    else:
                        # dt = dt/2

                        num = 2*num
                        den = 2*den
                        self.bisect = int(np.log2(den))

                        dt = dt_base/den

                        print("BISECTION LEVEL: %s:"%int(round(dt_base/dt)))
                        # if int(round(dt_base/dt))==64:set_trace()
                        print("dt/dt_base=",dt/dt_base)
                        csvwriter.writerow(['','','Cutback'])
                tmp += 1


            else:

                #############################
                #### CONVERGED INCREMENT ####
                #############################
                # self.savefig(redo_count,azimut=[-90, -90],elevation=[0,0],distance=[0,0],times=[t0,t,tf],fintC=False,Hooks=True)

                Redo = False
                for ic, ctct in enumerate(self.contacts):
                    ctct.actives_prev.append(list(ctct.actives))
                    ctct.getCandidates(self.u, CheckActive = True, TimeDisp=False,tracing=tracing)    # updates Patches->BSs (always) -> candidatePairs (on choice)
                    # ctct.getCandidatesANN(self.u, CheckActive = True, TimeDisp=False,tracing=tracing)    # updates Patches->BSs (always) -> candidatePairs (on choice)

                    acts_bool = [True if el is not None else False for el in ctct.actives]


                    if ctct.actives_prev[-1]!=ctct.actives: 
                        print("Actives have changed! Will Redo increment...")
                        Redo = True

                        pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
                        with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                            csvwriter = csv.writer(csvfile)
                            entered = list(set(ctct.actives).difference(set(ctct.actives_prev[-1])))
                            exited  = list(set(ctct.actives_prev[-1]).difference(set(ctct.actives)))

                            csvwriter.writerow(['','','RedoAct','IN']+entered+['OUT']+exited)

                        if ctct.actives in ctct.actives_prev:
                            print("CYCLE DETECTED!!! --> Do Minimization")
                            DoMinimization = True
                    # else:
                    #     ctct.actives_prev.append(list(ctct.actives))

                # # CYCLE DETECTION. TODO: Detect based on active set, not in residual. So that it doesn't computed again something repeated
                # if res in residuals_incr:
                #     ctct.SolveCycle()
                # else:
                #     residuals_incr.append(res)
                #     u_found_incr.append(self.u.copy())



                self.printContactStates(veredict=True)
                
                # if not Redo:
                #     for ctct in self.contacts:
                #         if ctct.checkGNs(): 
                #             Redo = True
                #             self.Redos_due_to_kn_adjustment += 1

                #             pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
                #             with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                #                 csvwriter = csv.writer(csvfile)
                #                 csvwriter.writerow(['','','RedoPen'])



                if Redo:        # don't write "else". It could have changed inside previous block
                    redo_count += 1
                    t-=dt
                    num -= 1
                    # skipBCs=True
                    tmp += 1
                    # u_ref_redo = self.u.copy()
                    self.getReferences()
                    # self.u = u_ref_redo.copy()

                    continue


                ti += 1
                for ctct in self.contacts:
                    ctct.actives_prev = []

                for body in self.bodies:
                    if not body.isRigid:
                        body.FPconv = body.FPtemp.copy()
                        body.EPcum += body.DELTA_EPcum

                DoMinimization = False

                redo_count = 0
                residuals_incr = []
                objectives_found = []
                u_found_incr = []

                for ib, body in enumerate(self.bodies):
                    sed_incr = body.get_nodal_SED(self.u)
                    self.SED[ib].append(sed_incr)


                print("##################")
                self.savefig(ti,azimut=[-90, -90],elevation=[0,0],distance=[10,10],times=[t0,t,tf],fintC=False,Hooks=True)
                # self.savefig(ti,azimut=[-90, -90],elevation=[0,0],distance=[10,10],times=[t0,t,tf],fintC=False,Hooks=True)

                self.setReferences()      

                sBody = self.contacts[0].slaveBody
                sNodes_diri = []
                for node in range(len(sBody.X)):
                    if sBody.DoFs[node][0] in self.diri:
                        sNodes_diri.append(node)
                self.saveNodalData(sBody.id, t, "fint", nodes = sNodes_diri, Sum=True)
                # self.savedata(t,'u','bisect','Redos_due_to_kn_adjustment')
                self.savedata(t,'u','bisect')




                self.saveContactData(0,t)

                self.Redos_due_to_kn_adjustment = 0

                if dt != dt_base:
                    # ratio = t/dt_base - int((t-t0)/dt_base)
                    # dt = dt_base/(float_to_fraction(ratio)[1])
                    fin = False
                    while fin==False:
                        if num%2==0:
                            num /= 2
                            den /= 2
                        else:
                            fin = True
                    num = round(num)    # could have changed to float during simplification but it's still an integer
                    den = round(den)    # could have changed to float during simplification but it's still an integer
                    self.bisect = int(np.log2(den))
                    dt = dt_base/den

                # saving:
                pickle.dump([self.REF,t, dt, ti,[num,den]],open(self.output_dir+"RecoveryData.dat","wb"))

        else:
            print(self.u)
            print("Increments Redone : ",tmp)

