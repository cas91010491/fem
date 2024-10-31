from PyClasses.BoundaryConditions import *
from PyClasses.Utilities import float_to_fraction, printif,brents_method, checkVarsSize, flatList,plot_coords,quadratic_fit_min_zeros
import PyClasses.FEAssembly
from scipy import sparse
from scipy.sparse.linalg import spsolve, cgs, bicg, gmres, lsqr
from scipy.linalg import solve as slsolve
from scipy.optimize import minimize, brentq
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import time, pickle, sys, os, csv

# import petsc4py
# petsc4py.init(sys.argv)
# from petsc4py import PETSc

from pdb import set_trace

class FEModel:
    def __init__(self,bodies, contacts, BoundaryConditions, UpdateOnIteration = True,transform_2d=None,subname=""):
        self.Mainname = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        self.bodies = bodies
        self.contacts = contacts
        self.UOI = UpdateOnIteration
        self.transform_2d = transform_2d

        # main counter during the simulation
        self.COUNTS_NAMES = ["incr_accptd", "incr_rjctd", "incr_mnmzd", "NR_iters", "mnmzn_iters", "mnmzn_fn_evals"]
        self.COUNTS = np.zeros((6,),dtype = int)

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

        dir =  "OUTPUT_"+current_time+self.Mainname+subname+"/"
        if not os.path.exists(dir):
            os.mkdir(dir)
        self.output_dir = dir
        pickle.dump(self,open(dir+"Model.dat","wb"))


    def printContactStates(self,only_actives=True,veredict=False):
        for i_c, contact in enumerate(self.contacts):
            print("contact",i_c,":")
            contact.printStates(only_actives=only_actives,veredict=veredict)

    def plot(self, ax, ref = 10, specialPatches = None, almostSpecialPatches = None,
              specialNodes = None, text = None, fintC = False, plotHooks=False,OnlyMasterSurf=False, u=None):
        
        if u is None:
            u = self.u

        for body in self.bodies:
            if body in [contact.slaveBody for contact in self.contacts]:
                # body.surf.plot(ax, u, specialNodes = specialNodes,ref=ref)
                body.surf.plot(ax, u, specialNodes = specialNodes,ref=1)
            elif body in [contact.masterBody for contact in self.contacts]:
                body.surf.plot(ax, u, specialPatches = specialPatches, almostSpecialPatches = almostSpecialPatches,ref=ref,onlyMaster=OnlyMasterSurf)
            else:
                body.surf.plot(ax, u,ref=ref)

        if fintC:
            for contact in self.contacts: contact.plotForces(u, ax)
        if plotHooks:
            for contact in self.contacts: contact.plotHooks(ax)

        if text is not None:
            ax.text2D(0.17, 0.7, text, transform=ax.transAxes)

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

    def savefig(self, increment, iteration = None,  elevation=[ 30 , 30], azimut=[-45 , -45], distance = [10, 10], times = None, dpi = 200, fintC = False, Hooks= False, u = None,simm_time=None):
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
        # ax.set_zlim3d((cz-0.8*hg, cz+0.8*hg))
        ax.set_zlim3d((0, 5))

        # # 3d Example:
        # ax.set_xlim3d((-2, 14))
        # ax.set_ylim3d((-4, 4))
        # ax.set_zlim3d((-2, 6))

        ax.set_aspect("equal")
        ax.xaxis._axinfo["grid"].update({"linewidth":0.1, "color" : "gray"})
        ax.yaxis._axinfo["grid"].update({"linewidth":0.1, "color" : "gray"})
        ax.zaxis._axinfo["grid"].update({"linewidth":0.1, "color" : "gray"})
        ax.tick_params(axis='z', which='both', direction='out', length=0)
        ax.set_yticks([])

        AlmostImportant = sorted(set(flatList(self.contacts[0].candids))-set(self.contacts[0].actives)) if len(self.contacts)!= 0 else []

        if times!= None:
            t0,t,tf = times
            tt = (t-t0)/(tf-t0)
        else:
            t = simm_time if simm_time is not None else 0
            tt = 0

        ctct = self.contacts[0]
        # activeNodes = [ctct.slaveNodes[idx] for idx in ctct.proj[:,3].nonzero()[0]]
        idxNodes = ctct.proj[:,3].nonzero()[0]
        activeNodes = [ctct.slaveNodes[idx] for idx in idxNodes if ctct.actives[idx] is not None ]
        # set_trace()
        activePatches = [pr[0] for pr in ctct.proj if pr[3]!=0]
        self.plot(ax, ref = 10,specialPatches = [(0,1,0,0.75),activePatches],
                              almostSpecialPatches = [(1.0,0.64,0.0,0.75), AlmostImportant],
                              specialNodes = {"red": activeNodes,}, 
                              text = "time: "+str(round(t,12)) + ((", iter:"+str(iteration)) if iteration is not None else ""),
                              fintC = fintC,
                              plotHooks = Hooks,
                              OnlyMasterSurf=True,
                              u = u
                              )

        #Camera effects
        ax.view_init(elev=elevation[0]*(1-tt)+elevation[1]*tt,azim=azimut[0]*(1-tt)+azimut[1]*tt)
        ax.dist = distance[0]*(1-tt)+distance[1]*tt     #default dist = 10
        ax.set_proj_type('ortho')
        # plot_coords(ax,orig=(minx-0.5+1.0,miny-0.5,minz-0.5))
        plot_coords(ax,orig=(-5.3,0.0,0.15))

        plt.tight_layout()

        if iteration==None:
            plt.savefig(self.output_dir+"plots/fig"+str(increment)+".png",dpi = dpi)
        else:
            plt.savefig(self.output_dir+"plots/fig"+str(increment)+"-"+str(iteration)+".png",dpi = dpi)
        plt.close()

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
            if self.IterUpdate:
                ctct.getfintC_unilateral(self, DispTime=DispTime,useANN=False)      # uses xs
            else:
                ctct.getfintC(self, DispTime=DispTime,useANN=False)      # uses xs


    def get_K(self, DispTime = False):
        t0_K = time.time()
        self.K = sparse.coo_matrix((self.ndof,self.ndof),dtype=float)
        for body in self.bodies:
            body.get_K(self)    # uses model.u_temp

        printif(DispTime,"Getting K : ",time.time()-t0_K,"s")

        for contact in self.contacts:
            if self.IterUpdate:
                contact.getKC_unilateral(self, DispTime=DispTime)     #uses model.u_temp
            else:
                contact.getKC(self, DispTime=DispTime)     #uses model.u_temp


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
            self.COUNTS[3] += 1

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
            self.COUNTS[3] += 1


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

    def  BFGS(self,FUNJAC, u0, tol = 1e-10, free_ind = None,ti = None,simm_time = None, plot = False):
        
        if free_ind is None:
            free_ind = self.free
        
        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)
        u = u0.copy()
        f_2 = f[free_ind]
        K_new_inv = np.eye(nfr)
        f_new = np.zeros(nfr)
        m0 = 0

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > tol:
            self.COUNTS[4] += 1

            iter += 1
            print("ITER:",iter)
            f_old = f_new.copy()
            f_new = f_2.copy()
            K_old_inv = K_new_inv.copy()
            delta_f = f_new - f_old

            if iter == 1:
                h_new = -np.dot(K_old_inv, f_new)
            else:
                K_new_inv = K_old_inv + ((np.inner(delta_u, delta_f) + np.inner(delta_f, np.dot(K_old_inv, delta_f)))*(np.outer(delta_u,delta_u)))/ (np.dot(delta_u, delta_f) ** 2)- (np.outer(np.dot(K_old_inv, delta_f),delta_u) + np.inner(np.outer(delta_u, delta_f),K_old_inv)) / np.dot(delta_u, delta_f)
                h_new = -np.dot(K_new_inv, f_new)

           
            if not np.isfinite(norm(h_new)):
                set_trace()

            a2,f2,f_2,ux = self.linesearch(FUNJAC, u, h_new, f_new, free_ind, alpha_init=1, c_par2=0.9)

            self.write_m_and_f(0.0,norm(f_2),iter)
                       
            delta_u = ux[free_ind] - u[free_ind]
            u = ux

            if plot:
                if self.transform_2d is None:
                    # 3D case
                    self.savefig(ti,iter,distance=[10,10],u = Ns@Nt@ux,simm_time=simm_time)
                else:
                    # for 2D case
                    Ns,Nt = self.transform_2d
                    self.savefig(ti,iter,azimut=[-90, -90],elevation=[0,0],distance=[10,10],u = Ns@Nt@ux,simm_time=simm_time)
    
            print("\talpha:",a2,"\tf2:",f2,"\t\t|f_2|:",norm(f_2))

        return u, m0, iter , norm(f_2)

    def LBFGS(self,FUNJAC, u0, tol = 1e-10,nli=10, free_ind = None,ti = None,simm_time = None, plot=False):

        if free_ind is None:
            free_ind = self.free

        nfr = len(free_ind)

        f , m_new = FUNJAC(u0)

        u = u0.copy()
        f_2 = f[free_ind]
        m_new = norm(f_2)
        self.write_m_and_f(m_new,norm(f),0)
        f_new = np.zeros(nfr)
        m0 = 0
        DF = np.zeros((nfr,nli))
        DU = np.zeros((nfr,nli))
        rho = np.zeros((nli,1))

        
        for ctct in self.contacts:
            ctct.patch_changes = []

        iter = 0
        while np.linalg.norm(f_2 - f_new) > 0 and np.linalg.norm(f_2) > tol:
            self.COUNTS[4] += 1

            iter += 1
            print("ITER:",iter)
            f_old = f_new.copy()
            f_new = f_2.copy()
            delta_f = f_new - f_old

            if iter == 1:
                h_new = -f_new
            else:
                gamma = np.zeros((nli,1))
                h_new = f_new.copy()
                for j in reversed(range(max(0,nli-iter+1),nli)):
                    rho_j = rho[j]
                    delta_u = DU[:,j]
                    delta_f = DF[:,j]
                    gamma_j = rho_j*delta_u@h_new
                    h_new -= gamma_j*delta_f
                    gamma[j] = gamma_j


                for j in range(max(0,nli-iter+1),nli):
                    rho_j = rho[j]
                    delta_u = DU[:,j]
                    delta_f = DF[:,j]
                    gamma_j = gamma[j]
                    eta = rho_j*delta_f@h_new
                    h_new += (gamma_j - eta)*delta_u

                h_new = -h_new

            if not np.isfinite(norm(h_new)):
                set_trace()

            a2,f2,f_2,ux = self.linesearch(FUNJAC, u, h_new, f_new, free_ind)

            self.write_m_and_f(0.0,norm(f_2),iter)
                       
            delta_u = ux[free_ind] - u[free_ind]
            u = ux.copy()

            DU[:,:nli-1] = DU[:,1:nli].copy()
            DU[:,-1] = delta_u.copy()
            DF[:,:nli-1] = DF[:,1:nli].copy()
            DF[:,-1] = f_2-f_new
            rho[:nli-1] = rho[1:nli].copy()
            rho[-1] = 1/((f_2-f_new)@delta_u)

            if plot:
                if self.transform_2d is None:
                    # 3D case
                    self.savefig(ti,iter,distance=[10,10],u = Ns@Nt@ux,simm_time=simm_time)
                else:
                    # for 2D case
                    Ns,Nt = self.transform_2d
                    self.savefig(ti,iter,azimut=[-90, -90],elevation=[0,0],distance=[10,10],u = Ns@Nt@ux,simm_time=simm_time)


            print("\talpha:",a2,"\tf2:",f2,"\t\t|f_2|:",norm(f_2))

        return u, m0, iter , norm(f_2)

    def linesearch(self, FUNJAC, u, h_new, f_new, free_ind, alpha_init=1, c_par2=0.9):
        ###################
        ### LINE SEARCH ###
        ###################
        # from joblib import Memory
        # cache_dir = './cache_directory'
        # memory = Memory(cache_dir, verbose=0)
        
        # @memory.cache
        def func(alpha,FUNJAC,u,free_ind,h_new,checkSecant=False):
            ux = u.copy()
            ux[free_ind] = u[free_ind] + alpha*h_new
            if checkSecant:
                fb,fc,ft,_,_ = FUNJAC(ux,split=True)
                fb = fb[free_ind]
                fc = fc[free_ind]
                ft = ft[free_ind]
                fb = h_new@fb
                fc = h_new@fc
                ft = h_new@ft

                return fb,fc,ft

            f , _ = FUNJAC(ux)
            f_3 = f[free_ind]
            f3 = h_new@f_3

            # self.ALs.append(alpha)
            # self.FFs.append(f3)
            return f3, f_3, ux

        # self.ALs = []   # Stores the values found during linesearch
        # self.FFs = []   # Stores the values found during linesearch

        a1 = 0
        f1 = h_new@f_new
        f_1 = f_new

        tol2 = min(abs(f1)/100 , 1e-12)

        a2 = 2*alpha_init
        f2 = np.nan
        while np.isnan(f2):
            a2 /= 2
            f2,f_2,ux = func(a2,FUNJAC,u,free_ind,h_new)

        ready_to_bisect = f2>0
        min_found = abs(f2)<tol2 and (np.dot(h_new, f_2) >= c_par2 * np.dot(h_new, f_new))

        print("\talphas:",[a1,a2],"\tf",[f1,f2])

        if not ready_to_bisect and not min_found:

            # Compute (a valid) a3
            delta = 2*(a2-a1)
            f3 = np.nan
            while np.isnan(f3):                
                delta /= 2
                a3 = a2 + delta                
                f3,f_3,ux = func(a3,FUNJAC,u,free_ind,h_new)

            print("\talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])

            ready_to_bisect = f3>0
            if ready_to_bisect:
                a1,f2,f_1 = a2,f2,f_2   
                a2,f2,f_2 = a3,f3,f_3       # from here it will go directly to bisection
            min_found = abs(f3)<tol2 and (np.dot(h_new, f_3) >= c_par2 * np.dot(h_new, f_new))

            while not ready_to_bisect and not min_found:

                try:
                    parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])  # if singular or no-zeros, goes to 'except'
                    if parab["a"]>0:
                        a0 = max(parab["zeros"])
                    else:
                        a0 = min(parab["zeros"])    # >0 already guaranteed by previous while loop
                    is_alpha_pos = a0>0
                except:
                    parab = {'zeros':None}
                    is_alpha_pos = False # This will force to search to the right (in while loop below)

                # Make sure the parabola crosses X in ascending manner and that cross is at alpha>0.
                while parab["zeros"] is None or not is_alpha_pos:
                    print("\tparabola not ready. Moving to right...")
                    # Let's move to the right
                    delta = 1.5*2*(a2-a1) # step 50% bigger each time and initially double because it will do at least one bisection
                    a1 = a2
                    f1 = f2

                    a2 = a3
                    f2 = f3

                    # Compute (a valid) a3
                    f3 = np.nan
                    while np.isnan(f3):                
                        delta /= 2
                        a3 = a2 + delta                
                        f3,f_3,ux = func(a3,FUNJAC,u,free_ind,h_new)

                    ready_to_bisect = f3>0
                    min_found = abs(f3)<tol2 and (np.dot(h_new, f_3) >= c_par2 * np.dot(h_new, f_new))
                    if ready_to_bisect or min_found:
                        a1,f2,f_1 = a2,f2,f_2   
                        a2,f2,f_2 = a3,f3,f_3
                        break

                    print("\talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])
                    try:
                        parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])  # if singular or no-zeros, goes to 'except'
                        if parab["a"]>0:
                            a0 = max(parab["zeros"])
                        else:
                            a0 = min(parab["zeros"])    # >0 already guaranteed by previous while loop
                        is_alpha_pos = a0>0
                    except:
                        parab = {'zeros':None}
                        is_alpha_pos = False # This will force to search to the right (in while loop below)

                if ready_to_bisect or min_found:
                    # this would mean they were found in the previous while loop
                    break

                delta = 2*(a0-a3)
                f0 = np.nan
                while np.isnan(f0):
                    delta /= 2
                    a0 = a3 + delta
                    f0,f_0,ux = func(a0,FUNJAC,u,free_ind,h_new)
                
                ready_to_bisect = f0>0
                min_found = abs(f0)<tol2 and (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new))
                print("\talphas:",[a1,a2,a3,a0],"\tf",[f1,f2,f3,f0])

                if min_found or ready_to_bisect:
                    a1,f2,f_1 = a3,f3,f_3
                    a2,f2,f_2 = a0,f0,f_0
                    break
                
                a1,f1,f_1 = a2,f2,f_2
                a2,f2,f_2 = a3,f3,f_3
                a3,f3,f_3 = a0,f0,f_0


        ####################
        # Quadratic Search #
        ####################
        secant = False
        bisection = False
        if not min_found:

            print("\t## Performing quadratic secant search ##")
            a3,f3,f_3 = a2,f2,f_2
            
            a2 = (a1+a3)/2
            f2,f_2,ux = func(a2,FUNJAC,u,free_ind,h_new)

            f0,f_0 = f1,f_1 # to enter the loop

            qbs_iter = 0
            while not (abs(f0)<tol2 and (np.dot(h_new, f_0) >= c_par2 * np.dot(h_new, f_new))):
                qbs_iter += 1

                if qbs_iter>20 and f2>0:
                    secant = True
                    break

                try:
                    parab = quadratic_fit_min_zeros([[a1,f1],[a2,f2],[a3,f3]])  # if singular or no-zeros, goes to 'except'
                    if parab["a"]>0:
                        a0 = max(parab["zeros"])
                    else:
                        a0 = min(parab["zeros"])    # >0 already guaranteed by previous while loop

                except:
                    secant = True
                    break
                
                f0,f_0,ux = func(a0,FUNJAC,u,free_ind,h_new)
                print("\talphas:",[a1,a2,a3],[a0],"\tf",[f1,f2,f3],[f0])

                if a0>a3:    # In this case f1,f2,f3 are all negative. since f0~0 and f is increasing in the interval
                    a1,f1,f_1 = a2,f2,f_2
                    a2,f2,f_2 = a3,f3,f_3
                    a3,f3,f_3 = a0,f0,f_0

                elif a0>a2:      # In this case f1,f2<0,  f3>0

                    d01 = np.sqrt((a0-a1)**2+(f0-f1)**2)
                    d03 = np.sqrt((a0-a3)**2+(f0-f3)**2)

                    if (f0>0 or abs(f3/f0)>4) and d03/d01>1:
                        # in this case f0 is better to be kept than f3
                        a3,f3,f_3 = a0,f0,f_0
                        
                    else:
                        a1,f1,f_1 = a2,f2,f_2
                        a2,f2,f_2 = a0,f0,f_0

                elif a0>a1:      # Here f2,f3>0 but f0 could be positive and in that case it should NOT replace f1

                    d01 = np.sqrt((a0-a1)**2+(f0-f1)**2)
                    d03 = np.sqrt((a0-a3)**2+(f0-f3)**2)

                    if (f0<0 or abs(f1/f0)>4) and d01/d03>1:
                        # in this case f0 is better to be kept than f1
                        a1,f1,f_1 = a0,f0,f_0
                        
                    else:
                        a3,f3,f_3 = a2,f2,f_2
                        a2,f2,f_2 = a0,f0,f_0

                else:
                    # Less likely to reach this part since f1<0 and f increases
                    a3,f3,f_3 = a2,f2,f_2
                    a2,f2,f_2 = a1,f1,f_1
                    a1,f1,f_1 = a0,f0,f_0


            if secant:

                ###################
                # Secant (linear) #
                ###################

                print("\t## parabola failed. Finishing off with Secant method ##")
                print("\tInitially, we have:")
                print("\t\talphas:",[a1,a2,a3],"\tf",[f1,f2,f3])

                sec_iter = 0
                while not (abs(f2)<tol2 and (np.dot(h_new, f_2) >= c_par2 * np.dot(h_new, f_new))) and f2 not in [f1,f3]:
                    if f2>0:
                        a3,f3,f_3 = a2,f2,f_2.copy()
                    else:
                        a1,f1,f_1 = a2,f2,f_2.copy()

                    a2 = (a1*f3 - a3*f1)/(f3-f1)
                    f2,f_2,ux = func(a2,FUNJAC,u,free_ind,h_new)
                    print("\talphas:",[a1,a3],[a2],"\tf",[f1,f3],[f2])

                    sec_iter +=1 
                    if sec_iter>10:
                        bisection = True
                        break
            else:
                a2,f2,f_2 = a0,f0,f_0

            if bisection:
                #############
                # Bisection #
                #############

                print("\t## secant failed. Finishing off with bisection method ##")
                print("\tInitially, we have:")
                print("\t\talphas:",[a1,a3],"\tf",[f1,f3])

                a2,f2,f_2 = a1,f1,f_1.copy()    # just to enter the loop
                bisect_iter = 0
                while not (abs(f2)<tol2 and (np.dot(h_new, f_2) >= c_par2 * np.dot(h_new, f_new))):
                    if f2>0:
                        a3,f3,f_3 = a2,f2,f_2.copy()
                    else:
                        a1,f1,f_1 = a2,f2,f_2.copy()

                    a2 = (a1+a3)/2
                    f2,f_2,ux = func(a2,FUNJAC,u,free_ind,h_new)
                    print("\talphas:",[a1,a3],[a2],"\tf",[f1,f3],[f2])

                    bisect_iter +=1 
                    if bisect_iter>40:
                        print("\tmaxiter reached! returning current value")
                        break

        # Here, one could plot all the values stored in self.ALs and self.FFs

        return a2,f2,f_2,ux

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

    def Energy_and_Force(self,u,split=False,show=False):
        self.COUNTS[5] += 1

        if self.transform_2d is not None:
            Ns,Nt = self.transform_2d
            u = (Ns@Nt@u)

        En,EnC = 0.0, 0.0
        force_body = np.zeros_like(self.fint)
        force_contact = np.zeros_like(self.fint)

        for body in self.bodies: 

            force_bi, Ebi = body.compute_mf_plastic(u,self)
            En += Ebi
            force_body+=force_bi

            # En+=body.compute_m(u)
        for ctct in self.contacts:
            # mCi,fCi = ctct.compute_mf(u,self)     # Bilateral constraint
            mCi,fCi = ctct.compute_mf_unilateral(u,self)
            EnC+=mCi
            force_contact+=fCi


        force = force_body + force_contact

        if show:
            print("Eb:",En,"\tEc:",EnC)

        if self.transform_2d is not None:
            if split:
                return Nt.T@Ns.T@force_body,Nt.T@Ns.T@force_contact,Nt.T@Ns.T@force, En,EnC

            return Nt.T@Ns.T@force, En+EnC

        if split:
            return force_body,force_contact,force, En,EnC
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


    def minimize(self,tol=1e-10,maxiter=10,plotIters=False,ti=None,simm_time=None,method = "BFGS", plot = False):

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

        if method == "BFGS":
            self.u, m_new, iter,res = self.BFGS(self.Energy_and_Force,u0,free_ind = fr,ti=ti,simm_time=simm_time,plot=plot)
        elif "LBFGS" in method:
            nli = int(method.replace("LBFGS",""))
            self.u, m_new, iter,res = self.LBFGS(self.Energy_and_Force,u0,nli=nli,free_ind = fr,ti=ti,simm_time=simm_time,plot=plot)


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


    def Solve(self, t0 = 0, tf = 1, TimeSteps = 1, max_iter=10 , recover = False, ForcedShift = False,
               IterUpdate = False,minimethod = "BFGS",plot=1):
        self.ntstps = TimeSteps
        self.IterUpdate = IterUpdate
        
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

        tracing = False

        if recover:
            self.REF,t, dt, ti,[num,den],self.COUNTS = pickle.load(open(recover))
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
                converged, res = self.minimize(tol=tolerance,ti=ti,simm_time=t,method=minimethod,plot=plot-1)
                self.COUNTS[2] += 1
                self.u_temp = np.array(self.u)  # copy of 'u' so that solution is directly used in NR
            else:
                converged, res = self.solve_NR(tol=tolerance,maxiter=max_iter)
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

                    if den >= MaxBisect:
                        print("MAXIMUM BISECTION LEVEL (%s) REACHED:"%MaxBisect)
                        # break
                        DoMinimization=True
                        csvwriter.writerow(['','','MINIMIZATION'])
                    else:

                        num = 2*num
                        den = 2*den
                        self.bisect = int(np.log2(den))

                        dt = dt_base/den

                        print("BISECTION LEVEL: %s:"%int(round(dt_base/dt)))
                        print("dt/dt_base=",dt/dt_base)
                        csvwriter.writerow(['','','Cutback'])

                for ctct in self.contacts:
                    ctct.actives_prev = []

                self.COUNTS[1] += 1

            else:

                #############################
                #### CONVERGED INCREMENT ####
                #############################
                Redo = False
                cycle_found = False
                for ic, ctct in enumerate(self.contacts):
                    # ctct.actives_prev.append(list(ctct.actives))    # At this point, 'actives' is the observed state
                    actives_prev = list(ctct.actives)    # At this point, 'actives' is the observed state
                    ctct.getCandidates(self.u, CheckActive = True, TimeDisp=False,tracing=tracing)    # updates Patches->BSs (always) -> candidatePairs (on choice)
                    # ctct.getCandidatesANN(self.u, CheckActive = True, TimeDisp=False,tracing=tracing)    # updates Patches->BSs (always) -> candidatePairs (on choice)

                    # acts_bool = [True if el is not None else False for el in ctct.actives]


                    # if ctct.actives_prev[-1]!=ctct.actives: 
                    if actives_prev!=ctct.actives: 
                        print("Actives have changed! Will Redo increment...")
                        Redo = True

                        pchfile = self.output_dir+"ctct"+str(ic)+"iters_details.csv"
                        with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                            csvwriter = csv.writer(csvfile)
                            entered = list(set(ctct.actives).difference(set(actives_prev)))
                            exited  = list(set(actives_prev).difference(set(ctct.actives)))

                            csvwriter.writerow(['','','RedoAct','IN']+entered+['OUT']+exited)


                        # CYCLE CHECK
                        if ctct.actives in ctct.actives_prev:
                            print("CYCLE DETECTED!!!")
                            t -= dt
                            num -= 1

                            self.getReferences(actives=True)    # resets u=u_temp <- u_ref

                            pchfile = self.output_dir+"ctct"+str(0)+"iters_details.csv"
                            with open(pchfile, 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                                csvwriter = csv.writer(csvfile)

                                if den >= MaxBisect:
                                    print("MAXIMUM BISECTION LEVEL (%s) REACHED:"%MaxBisect)
                                    # break
                                    DoMinimization=True
                                    csvwriter.writerow(['','','MINIMIZATION'])
                                else:

                                    num = 2*num
                                    den = 2*den
                                    self.bisect = int(np.log2(den))

                                    dt = dt_base/den

                                    print("BISECTION LEVEL: %s:"%int(round(dt_base/dt)))
                                    print("dt/dt_base=",dt/dt_base)
                                    csvwriter.writerow(['','','Cutback'])

                            self.COUNTS[1] += 1
                            cycle_found = True
                            for ctct in self.contacts:
                                ctct.actives_prev = []

                            break
                    

                if cycle_found:
                    continue

                ctct.actives_prev.append(list(ctct.actives))    # At this point, 'actives' is the observed state


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
                    # u_ref_redo = self.u.copy()
                    self.getReferences()
                    # self.u = u_ref_redo.copy()
                    self.COUNTS[1] += 1

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


                print("##################")

                if plot:
                    if self.transform_2d is None:
                        # 3D case
                        self.savefig(ti,distance=[10,10],times=[t0,t,tf],fintC=False,Hooks=True)
                    else:
                        # for 2D case
                        self.savefig(ti,azimut=[-90, -90],elevation=[0,0],distance=[10,10],times=[t0,t,tf],fintC=False,Hooks=True)


                self.setReferences()   

                self.COUNTS[0] += 1   

                sBody = self.contacts[0].slaveBody
                sNodes_diri = []
                for node in range(len(sBody.X)):
                    if sBody.DoFs[node][0] in self.diri:
                        sNodes_diri.append(node)
                self.saveNodalData(sBody.id, t, "fint", nodes = sNodes_diri, Sum=True)
                self.savedata(t,'u','bisect','Redos_due_to_kn_adjustment')




                with open(self.output_dir+"SED.csv", 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerow(self.bodies[0].get_nodal_SED(self.u).ravel().tolist())

                
                if self.bodies[0].plastic_param[0]<1e20:        # If it is actually plastic
                    with open(self.output_dir+"EPcum.csv", 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
                        csvwriter = csv.writer(csvfile)
                        csvwriter.writerow(self.bodies[0].EPcum.ravel().tolist())



                self.saveContactData(0,t)

                self.Redos_due_to_kn_adjustment = 0

                if dt != dt_base:
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
                pickle.dump([self.REF,t, dt, ti,[num,den],self.COUNTS],open(self.output_dir+"RecoveryData.dat","wb"))

        else:
            print("\n\n############################")
            print("### SIMULATION COMPLETED ###")
            print("############################\n")

        for val,count in zip(self.COUNTS,self.COUNTS_NAMES):
            print(count+":\t",val)


        with open(self.output_dir+"COUNTERS.csv", 'a') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(self.COUNTS_NAMES)
            csvwriter.writerow(self.COUNTS.tolist())




