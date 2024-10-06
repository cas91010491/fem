from unittest.mock import patch
from PyClasses.Utilities import *
from PyClasses.FEAssembly import *
from PyClasses.BoundingSpheres import *
from scipy import sparse
import concurrent.futures


from pdb import set_trace


class Contact:
    def __init__(self,slave, master, kn=1.0, kt=None,
                    cubicT=None, C1Edges=True, OutPatchAllowance=0.0, maxGN=None, mu=0, f0=0):          # initialized either with bodies or surfaces. Finally internally handled slave/master pair are surfaces.

        if type(slave) == list:        # pair [ body/surf , subsurf_nodes ]
            self.slaveBody = slave[0] if type(slave[0])==FEAssembly else slave[0].body
            self.slaveSurf = self.slaveBody.surf
            self.slaveNodes= slave[1]
        else:
            self.slaveSurf  = slave  if type(slave )==Surface else slave .surf
            self.slaveBody  = self.slaveSurf.body
            self.slaveNodes = self.slaveSurf.nodes

        if type(master) == list:        # pair [ body/surf , subsurf_nodes ]
            self.masterBody = master[0] if type(master[0])==FEAssembly else master[0].body
            self.masterSurf = self.masterBody.surf
            self.masterNodes= master[1]
            self.ipatches = self.quadsInNodeSelection(master[1])
        else:
            self.masterSurf = master if type(master)==Surface else master.surf
            self.masterBody = self.masterSurf.body
            self.masterNodes= self.masterSurf.nodes
            self.ipatches = range(len(self.masterSurf.quads))

        nsn = len(self.slaveNodes)
        self.nsn = nsn
        self.xs = np.array(self.slaveBody.X[self.slaveNodes])
        self.kn = kn
        self.kt = kn if kt is None else kt   # change to self.kn for updated stiffness based on gn_max
        self.f0 = f0*np.ones(nsn)
        self.mu = mu
        self.cubicT = cubicT
        self.C1Edges = C1Edges
        self.OPA = OutPatchAllowance            # extended patch surface to avoid chattering effect at edges. (In surface coordinates). currently not in use
        self.maxGN = maxGN
        self.minGN = None if maxGN is None else maxGN/2.0

        # projeccion search 
        self.candids = [[] for _ in range(nsn)]     # each node can have several candidates, hence we need DIFFERENT/UNIQUE empty lists
        self.actives = [None]*nsn

        # Normal state
        self.proj = np.zeros((nsn,4))    # keeps track of (patch, tx, ty, gn)  for each slave node  (t = (tx,ty) ) for current settings at each iteration
        self.alpha_p = np.ones(nsn)

        # Tangent state
        self.hook = np.zeros((nsn,4))    # keeps track of (patch,t0x,t0y,gn0)  for each slave node  (t0=(t0x,t0y)) for incremental initial settings
        self.slaveGTs   = np.zeros(nsn)
        self.Stick      = np.full (nsn,True)
        self.new      = np.full (nsn,True)  # This remains True until the iterarion where the spring starts acting (inclusive)

    def quadsInNodeSelection(self,nodes):
        ipatches = []
        for iquad, quad in enumerate(self.masterSurf.quads):
            quadInSelection = False
            for node in quad:
                if node not in self.masterNodes:
                    break
            else:
                quadInSelection = True
            if quadInSelection:
                ipatches.append(iquad)
        return ipatches


    def getCandidates(self, u, method = 'avgXmaxR', JustCP = False, TimeDisp=False,tracing=False):
        t0 = time.time()
        all_patches = self.masterSurf.patches
        self.candids = [[] for _ in range(self.nsn)]     # each node can have several candidates, hence we need DIFFERENT/UNIQUE empty lists

        # if tracing: set_trace()

        # updated positions of all slave nodes
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])
        # self.xs_pre = np.array(self.xs)
        self.xs = xs

        # Initiallizing GrgPatches if needed
        if all_patches[self.ipatches[0]] == None:
            MasterDoFs = self.masterBody.DoFs[self.masterNodes]
            xm = np.array(self.masterBody.X)[self.masterNodes] + np.array(u[MasterDoFs])
            BS_slave  = BS(xs, method = method)
            BS_master = BS(xm, method = method)
            if BS_slave.CollidesWithBS(BS_master):
                for ip in self.ipatches:
                    if not self.C1Edges:
                        self.masterSurf.patches[ip] = GrgPatch(self.masterSurf,ip, self.masterNodes)    # Takes SOME neighbours to avoid unwanted C1 at edges
                    else:
                        self.masterSurf.patches[ip] = GrgPatch(self.masterSurf,ip)  # Performs C1 smoothing at all edges

        else:
            for ipatch in self.ipatches:
                patch_obj = all_patches[ipatch]
                patch_obj.getCtrlPts(u)
                for ii in range(self.nsn):
                    if patch_obj.BS.ContainsNode(xs[ii]):
                        self.candids[ii].append(ipatch)

        printif(TimeDisp,"collisions checked in "+str(time.time()-t0)+ " s")

    def IsActive(self, idx, patch_id, OPA=None, tracing = False):
        """Checks whether there is penetration in a node-patch pair based on 't' and 'gn'. Returns bool. """
        surf = self.masterSurf
        patch = surf.patches[patch_id]
        xs = self.xs[idx]

        if tracing: set_trace()

        if patch == None:
            return 1

        t = patch.findProjection(xs, seeding=10)
        eps = self.OPA if OPA is None else OPA
        if not (0-eps<=t[0]<=1+eps and 0-eps<=t[1]<=1+eps):
            return False
        normal = patch.D3Grg(t)
        xc = patch.Grg(t)
        gn = (xs-xc) @ normal             # np.dot(a,b) <=> a @ b


        # if self.hook[idx,3] ==0:
        # else:
        #     if int(self.hook[idx,0])!=patch_id: return False
        #     t = self.hook[idx,1:3]
        #     xc0 = patch.Grg(t)
        #     n0  = patch.D3Grg(t)
        #     gn = (xs-xc0)@n0



        if gn >= 0:
            return False

        self.proj[idx] = np.array([patch_id,t[0],t[1],gn])
        return True


    def getfintC(self, Model, seeding = 10, DispTime = False, tracing = False, niter= 3):

        if DispTime: ti = time.time()
        self.actives = [None]*self.nsn
        self.proj = np.zeros((self.nsn,4))
        surf = self.masterSurf
        sBody = self.slaveBody

        for idx in range(self.nsn):
            for patch_id in self.candids[idx]:

                if self.IsActive(idx,patch_id, tracing=tracing):
                    self.actives[idx] = patch_id
                    node_id = self.slaveNodes[idx]
                    dofSlave = sBody.DoFs[node_id]
                    patch = surf.patches[patch_id]
                    dofPatch = surf.body.DoFs[patch.squad]
                    dofC = np.append(dofSlave,dofPatch)
                    xs = self.xs[idx]
                    kn  = self.alpha_p[idx]*self.kn

                    new = self.new[idx]

                    if self.hook[idx,3]==0 or self.mu==0:         # case where the node has just entered into contact (frictionless state is assumed)
                        # fintCN,gn,t_proj = patch.fintC_tot(xs,kn,self.kt,self.mu, cubicT=self.cubicT,disp=True)
                        fintCN = patch.fintC_tot2(xs,kn,self.kt,self.mu, cubicT=self.cubicT,disp=False)
                    else:
                        # fintCN, fintCT, gn, gt, t_proj  = patch.fintC_tot(xs,kn,self.kt,self.mu,hook=self.hook[idx],Stick=self.Stick[idx], cubicT=self.cubicT,disp=True, niter=niter, new = new)
                        fintCN, fintCT, self.slaveGTs[idx], self.Stick[idx] = patch.fintC_tot2(xs,kn,self.kt,self.mu,hook=self.hook[idx],Stick=self.Stick[idx], cubicT=self.cubicT,disp=False, niter=niter, new = new, tracing=False)
                        # FT = patch.torchC(xs,Model.u,self.kn,self.kt,self.mu,hook=self.hook[idx],Stick=self.Stick[idx], cubicT=self.cubicT, case=1,thing='f', tracing = tracing)
                        # if tracing: set_trace()
                        hook_patch = surf.patches[int(self.hook[idx,0])]
                        dofHPatch = surf.body.DoFs[hook_patch.squad]
                        dofCT = np.append(dofSlave,dofHPatch)
                        Model.fint[np.ix_(dofCT)] += fintCT

                    Model.fint[np.ix_(dofC)] += fintCN

        if DispTime: print("Getting fintC: ",time.time()-ti," s")


    def getKC(self, Model, seeding = 10, DispTime = False,tracing = False, niter=3):
        if DispTime: ti = time.time()
        surf = self.masterSurf
        sBody = self.slaveBody

        if tracing: set_trace()

        for idx in range(self.nsn):
            xs = self.xs[idx]
            if self.actives[idx] is None: continue
            patch_id = self.actives[idx]
            patch = surf.patches[patch_id]

            if self.proj[idx,3] != 0:
                node_id = self.slaveNodes[idx]
                dofSlave = sBody.DoFs[node_id]
                dofPatch = surf.body.DoFs[patch.squad]
                dofC = np.append(dofSlave,dofPatch)
                rn = np.repeat(dofC,len(dofC))
                cn = np.  tile(dofC,len(dofC))
                kn = self.alpha_p[idx]*self.kn

                new = self.new[idx]

                if self.hook[idx,3]==0 or self.mu==0:         # case where the node has just entered into contact (frictionless state is assumed)
                    KCN = patch.KC_tot2(xs,kn,self.kt,self.mu, cubicT=self.cubicT, niter=niter,new=new)
                else:
                    KCN,KCT,KCNT = patch.KC_tot2(xs,kn,self.kt,self.mu,hook=self.hook[idx],Stick=self.Stick[idx], cubicT=self.cubicT, niter=niter, new = new)

                    if KCT is not None:
                        hook_patch = surf.patches[int(self.hook[idx,0])]
                        dofHPatch = surf.body.DoFs[hook_patch.squad]
                        dofCT = np.append(dofSlave,dofHPatch)
                        rt = np.repeat(dofCT,len(dofCT))
                        ct = np.  tile(dofCT,len(dofCT))
                        sKCT  = sparse.coo_matrix((KCT.ravel()  ,(rt,ct)),shape=Model.K.shape)
                        Model.K += sKCT

                        if KCNT is not None:
                            rtn = np.repeat(dofCT,len(dofC ))
                            ctn = np.  tile(dofC ,len(dofCT))

                            sKCNT = sparse.coo_matrix((KCNT.ravel() ,(rtn,ctn)),shape=Model.K.shape)
                            Model.K += sKCNT


                sKC = sparse.coo_matrix((KCN.ravel(),(rn,cn)),shape=Model.K.shape)

                Model.K += sKC

        if DispTime: print("Getting KC: ",time.time()-ti," s")

    def checkStickSlip(self,impose=True,trace=False):
        """Updates 'Stick' atttribute for each slave node."""
        Redo = False
        for idx in range(self.nsn):
            kn,mu,gn = self.kn*self.alpha_p[idx], self.mu, self.proj[idx,3]
            N = abs(kn*gn)
            self.f0[idx] = mu*N
            if self.proj[idx,3]!=0:     # if there is penetration...
                node = self.slaveNodes[idx]
                if   (self.kt*self.slaveGTs[idx]<self.f0[idx]) and not self.Stick[idx]:
                    self.Stick[idx] = True
                    print("node ",node, " changing from Slip to Stick" if impose else "")
                    Redo = True
                elif (self.kt*self.slaveGTs[idx]>self.f0[idx]) and self.Stick[idx]:
                    self.Stick[idx] = 0
                    print("node ",node, " changing from Stick to Slip" if impose else "")
                    Redo = True
                print("T = ",self.kt*self.slaveGTs[idx],"\tmu*N = ",self.f0[idx])        
        if impose: return Redo

    def checkGNs(self):
        Redo = False
        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                node = self.slaveNodes[idx]
                gn = self.proj[idx,3]
                if abs(gn)>self.maxGN:
                    Redo = True
                    self.alpha_p[idx] *= (2.0*abs(gn)/self.maxGN)
                    print("too much penetration in node ",node)
                elif abs(gn)<self.minGN and self.alpha_p[idx]>1:
                    Redo = True
                    self.alpha_p[idx] = max(((2/3)*abs(gn)/self.minGN)*self.alpha_p[idx],1)       # alpha_p
                    print("too little penetration in node ",node)
            else: self.alpha_p[idx] = 1.0
        return Redo

    def UpdateHooks(self, fixedhook=False):
        for idx in range(self.nsn):
            if self.actives[idx] is not None:       # "actives" is the penetrated patch, updated in each iteration.
                if self.hook[idx,3]==0:    # if node JUST entered into contact...
                    self.hook[idx,0] = self.actives[idx]    #... create hook with current master patch
                    self.hook[idx,1:3] = self.proj[idx,1:3]
                    self.hook[idx,3] = -self.proj[idx,3]
                    self.new[idx]=True
                    # self.hook[idx,4] = not self.Stick[idx]

                elif not self.Stick[idx] and not fixedhook:  # if plastic limit is exceeded...
                    # compute new critical point within the elongated spring
                    p0_id = int(self.hook[idx,0])
                    p_id  = int(self.proj[idx,0])
                    t0 = self.hook[idx,1:3]
                    gn0 = self.hook[idx,3]
                    patch_hook = self.masterSurf.patches[p0_id]
                    patch_proj = self.masterSurf.patches[ p_id]
                    xc0 = patch_hook.Grg(t0)
                    n0 = patch_hook.D3Grg(t0)
                    xs = self.xs[idx]
                    # xs0 = xc0 - gn0*n0          #old method
                    xs0 = xc0 + np.dot(xs-xc0,n0)*n0
                    L = norm(xs0-xs)
                    gt_crit = self.f0[idx]/self.kt
                    # compute new patch,t0, gn0
                    xsm = xs - gt_crit*(xs-xs0)/L
                    midpatch = patch_proj
                    tm  = midpatch.findProjection(xsm)
                    if (not 0<=tm[0]<=1) or (not 0<=tm[1]<=1):
                        set_trace()
                        midpatch, p_id = patch_hook, p0_id
                        tm  = midpatch.findProjection(xsm)
                    xcm = midpatch.Grg(tm)
                    gnm = norm(xsm-xcm)

                    self.hook[idx] = [p_id,tm[0],tm[1],gnm]
                    # set (/verify?) Stick[idx] = True
                    self.Stick[idx] = True
                    self.new[idx] = False 
                else:
                    self.new[idx] = False 
            else:
                self.hook[idx] = 0.0   # if the node is (now) not active. There should be no data


    def plotForces(self,u,ax, factor = 40):
        if len(self.activePairs)==0:
            return None

        for (node,patch_id),fintC in zip(self.activePairs,self.fC):
            patch = self.masterSurf.patches[patch_id]
            dofs_node  = self.slaveBody.DoFs[node]
            dofs_patch = self.masterBody.DoFs[patch.squad]
            x_node = self.slaveBody.X[node] + u[dofs_node]
            x_patch = self.masterBody.X[patch.squad] + u[dofs_patch]
            for i, (fx,fy,fz) in enumerate(factor*fintC):
                x,y,z = x_node if i==0 else x_patch[i-1]
                ax.quiver(x, y, z, fx, fy, fz)

    def plotHooks(self,ax):
        for idx in range(self.nsn):
            if self.hook[idx,3]!=0:
                patch_obj = self.masterSurf.patches[int(self.hook[idx,0])]
                xc0 = patch_obj.Grg(self.hook[idx,1:3])
                xs = self.xs[idx]
                ax.scatter(xc0[0],xc0[1],xc0[2],color="black",marker='x',s=0.4)
                line = np.array([xc0,xs]).T
                ax.plot(line[0],line[1],line[2], color="black" if self.Stick[idx] else "red", lw = 0.2)

    def plotContact(self,u,ax,labels = True, SlaveQuads = False, context=False):
        "now I wanna plot all slave nodes + CANDIDATE patches"

        candidatePatches = flatList(self.candids)

        patch_transparency = 0.2
        for patch in self.masterSurf.patches:
            if (patch is None) or (not hasattr(patch,"CtrlPts")): continue
            # if not hasattr(patch,"CtrlPts"):continue
            pid = patch.iquad
            if pid in candidatePatches:
                color = (0,1,0,patch_transparency) if pid in self.actives else (1,0.5,0,patch_transparency)
                patch.plot(ax,color=color, label=labels)
            elif context:
                patch.plot(ax,color=(0.5,0.5,0.5,patch_transparency/2), label=labels)


        # for node in self.slaveNodes:
        for idx in range(self.nsn):
            # xs = self.slaveBody.X[node] + u[self.slaveBody.DoFs[node]]
            xs = self.xs[idx]
            color = "red" if self.proj[idx,3]!=0 else "blue"
            ax.scatter(xs[0],xs[1],xs[2],color=color)
            if labels:
                ax.text(xs[0],xs[1],xs[2],str(self.slaveNodes[idx]), color = color)

        if SlaveQuads:
            slaveQuads = self.slaveBody.SelectQuadsByNodes(self.slaveNodes)
            for quad in slaveQuads:
                self.slaveBody.surf.plotQuad(ax,u,quad,color=(0.5,0.5,0.5,0.5))


    def printStates(self):
        for idx in range(self.nsn):
            if len(self.candids[idx])>0:
                slavenode = self.slaveNodes[idx]
                candids = self.candids[idx]
                active = self.actives[idx]
                print("node ",slavenode,": \tcandidates :",candids,("\tactive: "+str(active)+
                ([" (slip)"," (stick)"][self.Stick[idx]] if self.hook[idx,3]!=0 else "")) if active is not None else "",
                "\tnew: "+str(self.new[idx]))
