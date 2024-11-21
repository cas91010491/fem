from unittest.mock import patch
from PyClasses.Utilities import *
from PyClasses.FEAssembly import *
from PyClasses.BoundingSpheres import *
from scipy import sparse
import concurrent.futures


from pdb import set_trace

# from thesis_sourcecode.src.model_training.patch_model_settings import PatchClassificationModel
# from thesis_sourcecode.src.model_training.surface_points_model_settings import SurfacePointsModelForOnePatch
# from thesis_sourcecode.src.model_training.signed_distance_model_settings import SignedDistanceModel


class Contact:
    def __init__(self,slave, master, kn=1.0, kt=None,
                    cubicT=None, C1Edges=True, OutPatchAllowance=0.0, maxGN=None, mu=0, f0=0, ANNmodel = None):          # initialized either with bodies or surfaces. Finally internally handled slave/master pair are surfaces.

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
        self.actives_prev = []
        self.candids_sparse = sparse.csr_matrix((self.nsn,len(self.masterBody.surf.patches)))
        self.actives_sparse = sparse.csr_matrix((self.nsn,len(self.masterBody.surf.patches)))
        # self.candids_sparse = np.zeros((self.nsn,len(self.masterBody.surf.patches)),dtype=bool)
        # self.actives_sparse = np.zeros((self.nsn,len(self.masterBody.surf.patches)),dtype=bool)
        self.candidpairs = np.array((10000,2),dtype=int)
        self.paired_T1T2 = np.array((10000,2),dtype=np.float64)

        # Normal state
        self.proj = np.zeros((nsn,4))    # keeps track of (patch, tx, ty, gn)  for each slave node  (t = (tx,ty) ) for current settings at each iteration
        self.alpha_p = np.ones(nsn)

        # Tangent state
        self.hook = np.zeros((nsn,4))    # keeps track of (patch,t0x,t0y,gn0)  for each slave node  (t0=(t0x,t0y)) for incremental initial settings
        self.slaveGTs   = np.zeros(nsn)
        self.Stick      = np.full (nsn,True)
        self.new      = np.full (nsn,True)  # This remains True until the iterarion where the spring starts acting (inclusive)

        self.patch_changes = []

        self.ANNmodel = ANNmodel

        patch_classifier_name = "final_patch_model_edges-shape-512-512-bs-64"
        # self.patch_classifier = PatchClassificationModel(name=patch_classifier_name)



    def quadsInNodeSelection(self,nodes):
        ipatches = []
        inclusive_selection = True     # if one node is selected, the whole quad is selected
        for iquad, quad in enumerate(self.masterSurf.quads):

            if inclusive_selection:
                quadInSelection = False
                for node in quad:
                    if node in self.masterNodes:
                        quadInSelection = True
                        break

            else:
                quadInSelection = False
                for node in quad:
                    if node not in self.masterNodes:
                        break
                else:
                    quadInSelection = True
            if quadInSelection:
                ipatches.append(iquad)
        return ipatches

    def getCandidatesANN(self, u,CheckActive = False,TimeDisp = False,tracing = False):
       
        if tracing:
            set_trace()

        t0 = time.time()

        self.candids = [[] for _ in range(self.nsn)]
        if CheckActive:
            self.actives = [None]*self.nsn


        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])
        self.xs = xs
        xs_forANN = xs + np.array([-6,0,0])

        # 1. get posibble active nodes using signed_distance_model
        signed_distance_regressor = SignedDistanceModel(name="testing_second_sd_model-shape-512-512-bs-64")
        predicted_signed_distance = signed_distance_regressor.predict(points=xs_forANN)
        proximal_candidates = np.where(predicted_signed_distance<0.05)[0]   # Slave nodes close to the master surface

        all_patches_obj = self.masterBody.surf.patches
        # self.candids_sparse = sparse.csr_matrix((self.nsn,len(all_patches_obj)))
        # self.actives_sparse = sparse.csr_matrix((self.nsn,len(all_patches_obj)))
        # self.candids_sparse.data.fill(0)
        # self.candids_sparse.eliminate_zeros()
        # self.actives_sparse.data.fill(0)
        # self.actives_sparse.eliminate_zeros()
        self.candids_sparse[:,:]=False
        self.actives_sparse[:,:]=False

        # 2. from possible active nodes get their candidates using model_predictions
        min_num_cands = 2       # 3 caused problems in a specific increment. So better 4 to cover corner cases (structured mesh)
        max_num_cands = 9
        min_categorization_output = 1e-6
        if len(proximal_candidates) > 0:
            model_predictions = self.patch_classifier.Predict(points=xs_forANN[proximal_candidates], n=max_num_cands)
            
            # Iterate through each node's predictions
            count = 0
            for ii, idx in enumerate(proximal_candidates):

                # Extract predictions for the current node
                node_predictions = model_predictions[ii]
                selected_patches = node_predictions[node_predictions[:, 1]>min_categorization_output]
                n_patches = max(min(len(selected_patches),max_num_cands),min_num_cands)
                if len(selected_patches) < min_num_cands:
                    selected_patches = node_predictions[:min_num_cands]       # this guarantees a minimum ammount of candidates

                # # Extend the candidates list with the selected patches
                # candidates[idx].extend(selected_patches[:, 0].astype(int).tolist())
                CandPatches_idx = selected_patches[:,0].astype(int).tolist()
                self.candids_sparse[idx,CandPatches_idx] = True
                # self.candidpairs[count:count+n_patches] = [[idx]*n_patches,CandPatches_idx]
                self.candids[idx] = CandPatches_idx


        # # T1 = sparse.csr_matrix((self.nsn,len(self.masterBody.surf.patches)))
        # # T2 = sparse.csr_matrix((self.nsn,len(self.masterBody.surf.patches)))
        T1 = -1*np.ones((self.nsn,len(self.masterBody.surf.patches)),dtype=np.float64)
        T2 = -1*np.ones((self.nsn,len(self.masterBody.surf.patches)),dtype=np.float64)

        for i_patch, patch in enumerate(all_patches_obj):
            nodes_in_patch = self.candids_sparse[:,i_patch].nonzero()[0]
            if len(nodes_in_patch)>0:
                predictions_for_t1t2 = patch.MinDistANN( xs_forANN[nodes_in_patch] , verbose=0 )
                T1[nodes_in_patch,i_patch] = predictions_for_t1t2[:,0]
                T2[nodes_in_patch,i_patch] = predictions_for_t1t2[:,1]



            if CheckActive:
                for ii in nodes_in_patch:
                        # self.candids[ii].append(i_patch)      # Unnecessary. Already done some lines above
                        t1 = T1[ii,i_patch]
                        t2 = T2[ii,i_patch]
                        active = self.IsActive(ii,i_patch,useANN=True,t0=None)
                        # active = self.IsActive(ii,i_patch,useANN=True,t0=None)
                        if active and self.actives[ii] is not None:
                            set_trace()
                            # If this ever happens, it means this node 'ii' seems to be
                            # active from the point of view of more than one patch.
                            # TODO:Make sure it chooses the one with significantly 
                            # higher prob. according to model_predictions[ii,1]
                            pass
                        elif active:
                            self.actives[ii]=i_patch
                            self.actives_sparse[ii,i_patch] = True


        printif(TimeDisp,"collisions checked in "+str(time.time()-t0)+ " s")

    def getCandidates(self, u, method = 'avgXmaxR', CheckActive = False, TimeDisp=False,tracing=False):
        # if tracing:
        #     set_trace()
        t0 = time.time()
        all_patches = self.masterSurf.patches
        self.candids = [[] for _ in range(self.nsn)]     # each node can have several candidates, hence we need DIFFERENT/UNIQUE empty lists
        if CheckActive:
            self.actives = [None]*self.nsn

        # updated positions of all slave nodes
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])
        self.xs = xs

        # Initiallizing GrgPatches if needed
        if all(element is None for element in all_patches):
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
                if all_patches[ipatch] is None:
                    continue
                patch_obj = all_patches[ipatch]
                if not self.masterSurf.body.isRigid:
                    patch_obj.getCtrlPts(u)

                for ii in range(self.nsn):
                    if patch_obj.BS.ContainsNode(xs[ii]):
                        self.candids[ii].append(ipatch)
                        if CheckActive:
                            if self.IsActive(ii,ipatch):
                                self.actives[ii]=ipatch


        for ii in range(self.nsn):
            if self.actives[ii] is None:
                self.proj[ii] = np.zeros((1,4))
        
        printif(TimeDisp,"collisions checked in "+str(time.time()-t0)+ " s")

    def updateActive(self):
        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                new_active=None
                for cand in self.candids[idx]:
                    if self.IsActive(idx,cand):
                        new_active=cand
                if new_active is None: 
                    print("No master patch anymore")
                    set_trace()
                    # if updatePatchChange:
                    #     if self.IsActive(ii,ipatch):
                            
                    # elif patch_obj.BS.ContainsNode(xs[ii]):
                    #     self.candids[ii].append(ipatch)
                    #     if CheckActive:
                    #         if self.IsActive(ii,ipatch):
                    #             self.actives[ii]=ipatch
    



    def IsActive(self, idx, patch_id, OPA=None,useANN=False,t0=None):
        """Checks whether there is penetration in a node-patch pair based on 't' and 'gn'. Returns bool. """
        surf = self.masterSurf
        patch = surf.patches[patch_id]
        xs = self.xs[idx]
        
        if patch == None:
            return 1

        if useANN and not patch.BS.ContainsNode(xs):   # This discards far nodes
            return False

        t = patch.findProjection(xs, seeding=10, ANNapprox = useANN, t0=t0)
        eps = self.OPA if OPA is None else OPA
        if not (0-eps<t[0]<1+eps and 0-eps<t[1]<1+eps):
            return False
        normal = patch.D3Grg(t)
        xc = patch.Grg(t)
        gn = (xs-xc) @ normal             # np.dot(a,b) <=> a @ b

        if gn >= 0:
            return False

        self.proj[idx] = np.array([patch_id,t[0],t[1],gn])
        return True



    def compute_mf(self, u, Model):
        surf = self.masterSurf

        m=0
        force=np.zeros(Model.fint.shape)
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_all = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])

        eventList_iter = []
        opa = self.OPA

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = xs_all[idx]
                kn  = self.alpha_p[idx]*self.kn
                
                patch_id = self.actives[idx]
                tried_updating_candidates = 0
                recursive_seeding = 1

                looper = 0
                is_patch_correct = False
                changed = False

                ANNapprox = False
                while not is_patch_correct:
                    patch = surf.patches[patch_id]

                    # if useANN:
                    #     t0 = [T1[idx,patch_id],T2[idx,patch_id]]
                    #     # If it's a decent candidate, evaluate
                    #     if (0-opaANN<t0[0]<1+opaANN and 0-opaANN<t0[1]<1+opaANN):
                    #         fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=useANN,t0=t0)
                    #         is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    # else:
                    mC,fintC,gn,t = patch.mf_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=ANNapprox,t0=None,recursive_seeding=recursive_seeding)
                    is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    # If not correct, try next candidate
                    if not is_patch_correct:
                        if looper==len(self.candids[idx]):  # No candidate is projecting well...

                            if tried_updating_candidates==1:
                                looper = 0
                                recursive_seeding = 100

                            elif tried_updating_candidates>1:
                                looper = 0
                                ANNapprox = True
                                
                                # TODO: add more elif<2,3,.. cases in which other methods are adopted to make sure we find projections for the slave
                                # for example include parameter "search_seeding" which in findProjections is 'recursive' to increase the search seeding 


                            elif tried_updating_candidates>2:
                                set_trace()
                                fintC = np.nan                      # <- this will force RedoHalf
                                break                               # ... and should be the last resource...
                            
                            self.getCandidates(u)
                            # set_trace()
                            # self.candids[idx] = self.patch_classifier.Predict(points=[xs+ np.array([-6,0,0])], n=9)[0,:,0].astype(int).tolist()
                            tried_updating_candidates += 1
                            looper = 0
                            # continue

                        changed = True
                        if len(self.candids[idx])>0:
                            patch_id = self.candids[idx][looper]    # this calls the next candidate patch
                            looper += 1
                        else:
                            # At this point the distortion brings the node too far from any patch
                            # patch_id = None
                            is_patch_correct = True
                            break

                if is_patch_correct:        # if patch changed
                    # node_id = self.slaveNodes[idx]
                    if changed:
                        eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->"+str(patch_id))

                    # if abs(gn)>10:
                    #     set_trace()

                    if gn>0:
                        if changed:
                            eventList_iter[-1]+=("out")     # ... if also changed patch ...
                        else:
                            eventList_iter.append(str(idx)+": out")     # ... or if only went out

                    self.actives[idx] = patch_id
                    # self.t1t2cache[idx] = t

                m += mC
                force[sDoFs[idx]] += fintC[:3]      # only for the slave node DoFs

        self.patch_changes=eventList_iter

        return m,force
    
    def compute_mf_unilateral(self, u, Model):
        surf = self.masterSurf

        m=0
        eventList_iter = []

        force=np.zeros(Model.fint.shape)
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        # xs_all = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])

        self.getCandidates(u)   # this updates self.xs


        opa = self.OPA
        for idx in range(self.nsn):
            # set_trace()
            changed = False
            xs = self.xs[idx]
            kn  = self.alpha_p[idx]*self.kn
            is_node_active = False

            # for patch_id, patch in enumerate(surf.patches):
            #     if patch is None: continue

            for patch_id in self.candids[idx]:
                patch = surf.patches[patch_id]

                recursive_seeding = 1
                mC,fintC,gn,t = patch.mf_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=False,t0=None,recursive_seeding=recursive_seeding)
                is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                if is_patch_correct:
                    if gn<0:
                        m += mC
                        force[sDoFs[idx]] += fintC[:3]      # only for the slave node DoFs

                        if patch_id != self.actives[idx]:
                            eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->"+str(patch_id))
                            changed = True

                        self.actives[idx] = patch_id
                        is_node_active = True
                        break
                    

            if not is_node_active:
                changed = self.actives[idx]!=None
                if changed :
                    eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->None")
                self.actives[idx] = None

        # print("actives:",self.actives)
        self.patch_changes=eventList_iter

        return m,force
    

    def compute_f(self, u, Model):
        surf = self.masterSurf

        force=np.zeros(Model.fint.shape)
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_all = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])

        eventList_iter = []
        opa = self.OPA

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = xs_all[idx]
                kn  = self.alpha_p[idx]*self.kn
                
                patch_id = self.actives[idx]
                tried_updating_candidates = 0

                looper = 0
                is_patch_correct = False
                changed = False

                ANNapprox = False
                while not is_patch_correct:
                    patch = surf.patches[patch_id]

                    # if useANN:
                    #     t0 = [T1[idx,patch_id],T2[idx,patch_id]]
                    #     # If it's a decent candidate, evaluate
                    #     if (0-opaANN<t0[0]<1+opaANN and 0-opaANN<t0[1]<1+opaANN):
                    #         fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=useANN,t0=t0)
                    #         is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    # else:
                    fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=ANNapprox,t0=None)
                    is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    # If not correct, try next candidate
                    if not is_patch_correct:
                        if looper==len(self.candids[idx]):  # No candidate is projecting well...
                            if tried_updating_candidates==1:
                                looper = 0
                                ANNapprox = True

                                # TODO: add more elif<2,3,.. cases in which other methods are adopted to make sure we find projections for the slave
                                # for example include parameter "search_seeding" which in findProjections is 'recursive' to increase the search seeding 


                            elif tried_updating_candidates>1:
                                set_trace()
                                fintC = np.nan                      # <- this will force RedoHalf
                                break                               # ... and should be the last resource...
                            self.candids[idx] = self.patch_classifier.Predict(points=[xs+ np.array([-6,0,0])], n=9)[0,:,0].astype(int).tolist()
                            tried_updating_candidates += 1
                            looper = 0
                            # continue

                        patch_id = self.candids[idx][looper]    # this calls the next candidate patch
                        looper += 1
                        changed = True

                if is_patch_correct:        # if patch changed
                    if changed:
                        eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->"+str(patch_id))

                    if abs(gn)>10:
                        set_trace()

                    if gn>0:
                        if changed:
                            eventList_iter[-1]+=("out")     # ... if also changed patch ...
                        else:
                            eventList_iter.append(str(idx)+": out")     # ... or if only went out

                    self.actives[idx] = patch_id
                    # self.t1t2cache[idx] = t

                
                force[sDoFs[idx]] += fintC[:3]      # only for the slave node DoFs

        self.patch_changes=eventList_iter

        return force
    

    def compute_f_backup(self, u, Model):
        surf = self.masterSurf

        force=np.zeros(Model.fint.shape)
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_all = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])

        eventList_iter = []
        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                changed = False

                xs = xs_all[idx]
                kn  = self.alpha_p[idx]*self.kn
                
                patch_id = self.actives[idx]
                is_patch_correct = False
                looper = 0

                while not is_patch_correct:
                    patch = surf.patches[patch_id]
                    fintC,gn,t = patch.fintC(xs,kn,cubicT=self.cubicT)
                    opa = self.OPA
                    is_patch_correct = 0-opa<=t[0]<=1+opa and 0-opa<=t[1]<=1+opa    #boolean
                    if not is_patch_correct:
                        if abs(looper)>len(self.candids[idx]):  # No candidate is projecting well...
                                                                # this means structure is distorted
                            fintC = np.zeros_like(force)
                            break
                        patch_id = self.candids[idx][looper]
                        looper -= 1
                        changed = True

                if is_patch_correct:        # if patch changed
                    if changed:
                        eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->"+str(patch_id))
                    if gn>0:
                        if not changed:
                            eventList_iter.append(str(idx)+": out")
                        else:
                            eventList_iter[-1]+=("out")

                    self.actives[idx] = patch_id
                    # self.t1t2cache[idx] = t

                
                force[sDoFs[idx]] += fintC[:3]      # only for the slave node DoFs

        self.patch_changes=eventList_iter
        
        return force

    def compute_m(self, u,t=None):
        surf = self.masterSurf
        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_temp = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(u[sDoFs ])

        mC = 0
        for idx in range(self.nsn):
            if self.actives[idx] is not None:

                xs = xs_temp[idx]
                patch_id = self.actives[idx]
                patch = surf.patches[patch_id]
                kn = self.alpha_p[idx]*self.kn

                # mC += patch.mC(xs,kn,cubicT=self.cubicT,t=self.t1t2cache[idx])
                mC += patch.mC(xs,kn,cubicT=self.cubicT,t=None)

        return mC

    def compute_k(self, u,t=None):
        ndofs = len(u)
        K = sparse.coo_matrix((ndofs,ndofs),dtype=float)

        surf = self.masterSurf
        sBody = self.slaveBody

        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_temp = np.array(self.slaveBody.X)[self.slaveNodes] + np.array(u[sDoFs ])

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = xs_temp[idx]
                patch_id = self.actives[idx]
                patch = surf.patches[patch_id]
                node_id = self.slaveNodes[idx]
                dofSlave = sBody.DoFs[node_id]
                dofPatch = surf.body.DoFs[patch.squad]
                dofC = np.append(dofSlave,dofPatch)
                rn = np.repeat(dofC,len(dofC))
                cn = np.  tile(dofC,len(dofC))
                kn = self.alpha_p[idx]*self.kn

                # KC = patch.KC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, t=self.t1t2cache[idx])
                KC = patch.KC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, t=None)
                sKC = sparse.coo_matrix((KC.ravel(),(rn,cn)),shape=K.shape)

                K += sKC

        return K



    def getfintC(self, Model, DispTime = False, useANN = False,tracing = False):

        # rct1t2 = np.zeros((10000,4))        # Here, I will store node, patch, t1,t2. To create the sparse matrices fast


        if tracing:
            set_trace()

        if DispTime: ti = time.time()
        surf = self.masterSurf
        sBody = self.slaveBody
        self.t1t2cache = -1*np.ones((self.nsn,2))

        if useANN:
            xs_forANN = self.xs + np.array([-6,0,0])
            all_patches_obj = surf.patches

            #  Get (fast) the initial guesses for t1t2 according to each patch
            # T1 = sparse.csr_matrix((self.nsn,len(all_patches_obj)))
            # T2 = sparse.csr_matrix((self.nsn,len(all_patches_obj)))
            # T1 = -1*np.ones((self.nsn,len(self.masterBody.surf.patches)),dtype=np.float64)
            # T2 = -1*np.ones((self.nsn,len(self.masterBody.surf.patches)),dtype=np.float64)

            rct_count = 0
            for i_patch, patch in enumerate(all_patches_obj):
                nodes_in_patch = self.candids_sparse[:,i_patch].nonzero()[0]
                num_nodes_in_patch = len(nodes_in_patch)
                if num_nodes_in_patch>0:
                    predictions_for_t1t2 = patch.MinDistANN( xs_forANN[nodes_in_patch] , verbose=0 )
                    # set_trace()
                    # T1[nodes_in_patch,i_patch] = predictions_for_t1t2[:,0]
                    # T2[nodes_in_patch,i_patch] = predictions_for_t1t2[:,1]
                    # rct1t2[rct_count:rct_count+num_nodes_in_patch]=np.array([nodes_in_patch,[i_patch]*num_nodes_in_patch])


        eventList_iter = []
        opa = self.OPA
        opaANN = 1e-2
        
        # print("Active Pairs:",end="")

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = self.xs[idx]
                kn  = self.alpha_p[idx]*self.kn
                
                patch_id = self.actives[idx]    # Current active patch

                tried_updating_candidates = False

                looper = 0      # To be used if the first candidate (self.candids[0] == self.active) is not correct
                is_patch_correct = False        # measures ONLY tangential correspondance
                changed = False
                while not is_patch_correct:
                    patch = surf.patches[patch_id]

                    if useANN:
                        t0 = [T1[idx,patch_id],T2[idx,patch_id]]
                        # If it's a decent candidate, evaluate
                        if (0-opaANN<t0[0]<1+opaANN and 0-opaANN<t0[1]<1+opaANN):
                            fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=useANN,t0=t0)
                            is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    else:
                        fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=False,t0=None)
                        is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                    # If not correct, try next candidate
                    if not is_patch_correct:
                        if looper==len(self.candids[idx]):  # No candidate is projecting well...
                            if tried_updating_candidates:
                                fintC = np.nan                      # <- this will force RedoHalf
                                break
                            #     self.actives[idx] = None
                            #     is_patch_correct = True # Just to stop search and accept node out
                            #     changed = False
                            #     gn = 1
                            #     patch_id = self.candids[0]
                            #                                     # this means structure is distorted

                            # for the 2D-case only!!
                            fintC = np.nan                      # <- this will force RedoHalf
                            break

                            self.candids[idx] = self.patch_classifier.Predict(points=[xs+ np.array([-6,0,0])], n=9)[0,:,0].astype(int).tolist()
                            tried_updating_candidates = True
                            looper = 0
                            # continue

                        patch_id = self.candids[idx][looper]    # this calls the next candidate patch
                        looper += 1
                        changed = True

                if is_patch_correct:        # if patch changed
                    node_id = self.slaveNodes[idx]
                    if changed:
                        eventList_iter.append(str(idx)+": "+str(self.actives[idx])+"-->"+str(patch_id))
                    if gn>0:
                        if changed:
                            eventList_iter[-1]+=("out")     # ... if also changed patch ...
                        else:
                            eventList_iter.append(str(idx)+": out")     # ... or if only went out

                    self.actives[idx] = patch_id
                    self.t1t2cache[idx] = t
                    # print("["+str(idx)+","+str(patch_id)+"],",end="")


                else:
                    print("Not correct patch found!!!")
                    return

                dofPatch = surf.body.DoFs[patch.squad]
                dofSlave = sBody.DoFs[node_id]
                dofC = np.append(dofSlave,dofPatch)

               
                Model.fint[np.ix_(dofC)] += fintC
        self.patch_changes=eventList_iter

        # print("")
        if DispTime: print("Getting fintC: ",time.time()-ti," s")

    def getfintC_unilateral(self, Model, DispTime = False, useANN = False,tracing = False):
        surf = self.masterSurf

        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_all = np.array(self.slaveBody.X )[self.slaveNodes ] + np.array(Model.u[sDoFs ])
        self.t1t2cache = -1*np.ones((self.nsn,2))
        eventList_iter = []


        opa = self.OPA
        for idx in range(self.nsn):
            # set_trace()
            xs = xs_all[idx]
            kn  = self.alpha_p[idx]*self.kn
            is_node_active = False
            changed = False
            for patch_id in self.candids[idx]:
                recursive_seeding = 1
                patch = surf.patches[patch_id]

                fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, ANNapprox=useANN,recursive_seeding=recursive_seeding)
                is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean

                if is_patch_correct:
                    changed = (patch_id!=self.actives[idx])
                    node_id = self.slaveNodes[idx]
                    if changed:
                        eventList_iter.append(str(node_id)+": "+str(self.actives[idx])+"-->"+str(patch_id))

                    if gn<0:
                        Model.fint[sDoFs[idx]] += fintC[:3]      # only for the slave node DoFs
                        self.actives[idx] = patch_id
                        is_node_active = True
                        self.t1t2cache[idx] = t

                        break

                    if changed:
                        eventList_iter[-1]+=("out")     # ... if also changed patch ...
                    else:
                        eventList_iter.append(str(node_id)+": out")     # ... or if only went out


            if not is_node_active:
                self.actives[idx] = None

        self.patch_changes=eventList_iter

        print("actives:",self.actives)


    def getfintC_backup(self, Model, DispTime = False):
        if DispTime: ti = time.time()
        surf = self.masterSurf
        sBody = self.slaveBody
        self.t1t2cache = -1*np.ones((self.nsn,2))


        eventList_iter = []
        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                node_id = self.slaveNodes[idx]
                changed = False
                dofSlave = sBody.DoFs[node_id]
                xs = self.xs[idx]
                kn  = self.alpha_p[idx]*self.kn
                
                patch_id = self.actives[idx]    # Current active patch
                is_patch_correct = False
                looper = 0

                while not is_patch_correct:
                    patch = surf.patches[patch_id]
                    fintC,gn,t = patch.fintC_fless_rigidMaster(xs,kn,cubicT=self.cubicT)
                    opa = self.OPA
                    is_patch_correct = 0-opa<=t[0]<=1+opa and 0-opa<=t[1]<=1+opa    #boolean
                    if not is_patch_correct:
                        if abs(looper)>len(self.candids[idx]):  # No candidate is projecting well...
                                                                # this means structure is distorted
                            fintC = np.nan                      # <- this will force RedoHalf
                            break
                        patch_id = self.candids[idx][looper]
                        looper -= 1
                        changed = True

                if is_patch_correct:        # if patch changed
                    if changed:
                        eventList_iter.append(str(node_id)+": "+str(self.actives[idx])+"-->"+str(patch_id))
                    if gn>0:
                        if changed:
                            eventList_iter[-1]+=("out")     # ... if also changed patch ...
                        else:
                            eventList_iter.append(str(node_id)+": out")     # ... or if only went out

                    self.actives[idx] = patch_id
                    self.t1t2cache[idx] = t

                else:
                    print("Not correct patch found!!!")
                    set_trace()

                dofPatch = surf.body.DoFs[patch.squad]
                dofC = np.append(dofSlave,dofPatch)

               
                Model.fint[np.ix_(dofC)] += fintC
        self.patch_changes=eventList_iter

        if DispTime: print("Getting fintC: ",time.time()-ti," s")

    def getKC(self, Model, DispTime = False):
        if DispTime: ti = time.time()
        surf = self.masterSurf
        sBody = self.slaveBody

        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_temp = np.array(self.slaveBody.X)[self.slaveNodes] + np.array(Model.u_temp[sDoFs ])

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = xs_temp[idx]
                patch_id = self.actives[idx]
                patch = surf.patches[patch_id]
                node_id = self.slaveNodes[idx]
                dofSlave = sBody.DoFs[node_id]
                dofPatch = surf.body.DoFs[patch.squad]
                dofC = np.append(dofSlave,dofPatch)
                rn = np.repeat(dofC,len(dofC))
                cn = np.  tile(dofC,len(dofC))
                kn = self.alpha_p[idx]*self.kn

                KC = patch.KC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, t=self.t1t2cache[idx])
                sKC = sparse.coo_matrix((KC.ravel(),(rn,cn)),shape=Model.K.shape)

                Model.K += sKC

        if DispTime: print("Getting KC: ",time.time()-ti," s")

    def getKC_unilateral(self, Model, DispTime = False):
        if DispTime: ti = time.time()
        surf = self.masterSurf
        sBody = self.slaveBody

        sDoFs  = self.slaveBody.DoFs[self.slaveNodes]
        xs_temp = np.array(self.slaveBody.X)[self.slaveNodes] + np.array(Model.u_temp[sDoFs ])

        for idx in range(self.nsn):
            if self.actives[idx] is not None:
                xs = xs_temp[idx]
                patch_id = self.actives[idx]
                patch = surf.patches[patch_id]
                node_id = self.slaveNodes[idx]
                dofSlave = sBody.DoFs[node_id]
                dofPatch = surf.body.DoFs[patch.squad]
                dofC = np.append(dofSlave,dofPatch)
                rn = np.repeat(dofC,len(dofC))
                cn = np.  tile(dofC,len(dofC))
                kn = self.alpha_p[idx]*self.kn

                KC = patch.KC_fless_rigidMaster(xs,kn,cubicT=self.cubicT, t=self.t1t2cache[idx])
                sKC = sparse.coo_matrix((KC.ravel(),(rn,cn)),shape=Model.K.shape)

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
                    # self.alpha_p[idx] *= (2.0*abs(gn)/self.maxGN)
                    self.alpha_p[idx] *= (1.1*abs(gn)/self.maxGN)
                    print("too much penetration in node ",node)
                elif abs(gn)<self.minGN and self.alpha_p[idx]>1:
                    Redo = True
                    # self.alpha_p[idx] = max(((2/3)*abs(gn)/self.minGN)*self.alpha_p[idx],1)       # alpha_p
                    self.alpha_p[idx] = max((0.9*abs(gn)/self.minGN)*self.alpha_p[idx],1)       # alpha_p
                    print("too little penetration in node ",node)
            else: 
                self.alpha_p[idx] = 1.0
                self.proj[idx,3] = 0.0
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

    def SolveCycle(self):

        print("Trying to exit the cycle...")

        def result(list_of_actives):
            import random
            from collections import Counter

            count = Counter(list_of_actives)
            max_count = max(count.values())
            most_frequent = [key for key, val in count.items() if val == max_count]
            result = random.choice(most_frequent)

            return result
        
        def asbool(list_of_actives):
            boollist = []
            for actives_i in self.actives_prev:
                boollist.append([act is not None for act in actives_i])

            return boollist


        nsn = len(self.slaveNodes)
        actives_prev_bool = asbool(self.actives_prev)
        final_state_bool = actives_prev_bool[-1]

        while final_state_bool in actives_prev_bool:
            final_state_bool = np.random.choice(a=[False, True],size=nsn).tolist()  # It must be in list to compare in 'while' condition

        final_state = [None]*nsn
        
        for i in range(nsn):
            if final_state_bool[i]:
                actives = []
                for actives_i in self.actives_prev:
                    if actives_i[i] is not None:
                        actives.append(actives_i[i])
                if len(actives)>0:
                    final_state[i] = result(actives)

        print("New active set proposed:", final_state)

        self.actives = final_state



    def printStates(self,only_actives=True,veredict=False):
        for idx in range(self.nsn):
            if len(self.candids[idx])>0:
                slavenode = self.slaveNodes[idx]
                active = self.actives[idx]
                if veredict:
                    changed = active != self.actives_prev[-1][idx]
                    if changed: 
                        print("node ",slavenode,":\t",self.actives_prev[-1][idx],"-->",active)
                else:
                    candids = self.candids[idx]
                    if not only_actives or active is not None:
                        print("node ",slavenode,": \tcandidates :",candids,
                            "" if active is None else ("\tactive: "+str(active)) )

