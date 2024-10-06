from PyClasses.Utilities import *
from PyClasses.GregoryPatches import *
import PyClasses

import time
import numpy as np
from numpy.linalg import norm
import pickle

from pdb import set_trace


class Surface:
    def __init__(self, *constructor, outers = None):   # gets quads and nodes on the surface and for each node their quads and neighbor nodes
        
        # Body-surface (automatically made by FEAssembly initiallizator)
        if len(constructor)==1 and type(constructor[0]) == PyClasses.FEAssembly.FEAssembly:
            self.body = constructor[0]
            self.X = constructor[0].X
            self.body.surf = self
            if  outers is not None:
                self.quads, self.nodes = outers
            else:
                self.quads, self.nodes = self.GetOuters()
                pickle.dump([self.quads,self.nodes],open("RecoveryOuters.dat","wb"))
            self.sub_surfs = []
            self.isSub = False
            self.patches = [None]*len(self.quads)

        # Surface out of pair (X, quads)
        elif type(constructor[0])in [list,np.ndarray] and len(constructor[1][0])==4:
            self.X = constructor[0]
            self.quads = constructor[1]
            self.nodes = list(set(flatList(self.quads)))
        else:
            set_trace()
        self.DoFs = []              # assigned by FEModel constructor
        self.NodeInQuads = self.GetQuadsPerNode()
        self.NeighborNodes = self.GetNeighborNodes()
        self.patches = [None]*len(self.quads)

    # def GetOuters(self):                        # Called by __init__
    #     quadss = []
    #     t0 = time.time()
    #     for nodes in self.body.hexas:
    #         bottom = list(nodes[0:4][::-1])
    #         top   = list(nodes[4:8])
    #         side1 = [nodes[0],nodes[1],nodes[5],nodes[4]]
    #         side2 = [nodes[1],nodes[2],nodes[6],nodes[5]]
    #         side3 = [nodes[2],nodes[3],nodes[7],nodes[6]]
    #         side4 = [nodes[3],nodes[0],nodes[4],nodes[7]]
    #         quadss.extend([bottom,top,side1,side2,side3,side4]) # one Quad per side (6 in total)
    #     iq1=0
    #     while iq1<len(quadss):
    #         iq2=iq1+1
    #         while iq2<len(quadss):
    #             if set(quadss[iq1]) == set(quadss[iq2]):
    #                 del quadss[iq1]
    #                 del quadss[iq2-1]
    #                 iq1 -= 1
    #                 break
    #             iq2 += 1
    #         iq1 += 1
    #     tf = time.time()
    #     print("Computing outer Quads : Done [ ",time.time()-t0," s ]")
    #     nodes = list(set(flatList(quadss)))

    #     return quadss, nodes

    def GetOuters(self):                        # Called by __init__
        t0 = time.time()
        from collections import defaultdict 

        track = dict()
        lookup = defaultdict(bool)
        for node in self.body.hexas:#<-remove # and it should work as a method 
            bottom = node[:4]; bottom.sort()
            top   = node[4:]; top.sort()
            side1 = [*node[:2],*node[4:6]]; side1.sort()
            side2 = [*node[1:3],*node[5:7]]; side2.sort()
            side3 = [*node[2:4],*node[6:8]]; side3.sort()
            side4 = [*node[3:5],node[0],node[7]]; side4.sort()

            lookup[bottom.__str__()] ^= True
            lookup[top.__str__()] ^= True
            lookup[side1.__str__()] ^= True
            lookup[side2.__str__()] ^= True
            lookup[side3.__str__()] ^= True
            lookup[side4.__str__()] ^= True

            track[bottom.__str__()] = [node[3], node[2], node[1], node[0]]
            track[top.__str__()]    = node[4:]
            track[side1.__str__()]  = [node[0], node[1], node[5], node[4]]
            track[side2.__str__()]  = [node[1], node[2], node[6], node[5]]
            track[side3.__str__()]  = [node[2], node[3], node[7], node[6]]
            track[side4.__str__()]  = [node[3], node[0], node[4], node[7]]

        quads = []
        print("Computing outer Quads : Done [ ",time.time()-t0," s ]")

        for (key,val) in lookup.items():
            val and quads.append(track[key])

        nodes = {*quads[0]}
        for quad in quads:
            nodes.update(quad)
        return quads, list(nodes)





    def GetQuadsPerNode(self):
        # Find quads related to each node
        NodeInQuads = [[] for _ in range(len(self.nodes))]   # list of quads in which each node is present
        t0 = time.time()
        for quadid, quad in enumerate(self.quads):      #quad is a list of 4 numbers (nodesIds)
            for nodeid in quad:
                loc_nid = self.nodes.index(nodeid)
                NodeInQuads[loc_nid].append(quadid)
        tf = time.time()
        print("relating Quads to each Node: Done [ "+str(tf-t0)+" s ]")
        return NodeInQuads

    def GetNeighborNodes(self):
        # Get NeighborNodes
        NeighborNodes = [[] for _ in range(len(self.nodes))]      # list of neighbouring nodes along the edges
        t0 = time.time()

        for nodeid_local, nodeid in enumerate(self.nodes):
            for quadid in self.NodeInQuads[nodeid_local]:
                OutQuad = self.quads[quadid] 
                idx = OutQuad.index(nodeid)
                NeighborNodes[nodeid_local].append([OutQuad[idx-3],OutQuad[idx-1]])     #TODO: Check the normals using this. Is this the right order?
        tf = time.time()
        print("Finding Pairs for each Node: Done [ "+str(tf-t0)+" s ]")
        return NeighborNodes

    def ComputeGrgPatches(self, u, nodes,compute_deriv=True):             # Computes nodal normals and all 20 Control points for each patch

        # Finding all quads involved according to inputted nodes (quadsExt) and the repective extension of nodes involved (nodesExt)
        quads = []              # quads where GP have to be (re)computed
        for node in nodes:
            quads.extend(self.NodeInQuads[node])
        quadsExt = list(set(quads))

        self.GrgPatch_idxs = quadsExt

        for iq, quad in enumerate(self.quads):

            # TODO: Change this in order to update instead of creating new patches all the time

            if iq in quadsExt:
                if self.patches[iq] == None:
                    self.patches[iq] = GrgPatch(self, iq)
                self.patches[iq].getCtrlPts(u, compute_deriv=compute_deriv)

            else:
                self.patches[iq] = None

    def QuadInterpolation(self,u,iquad,xv,yv):
        quad = np.array(self.quads[iquad] ) 
        X = np.array(self.X)[quad]

        if len(u) == len(self.body.X):
            u = np.reshape(u,(-1,3))[quad]  #CHANGE BACK!
        else:
            u = np.reshape(u[self.body.DoFs],(-1,3))[quad]  #CHANGE BACK!
        xa, xb, xc, xd = np.array(X+u)
        xab = xa*(1-xv) + xb*xv
        xdc = xd*(1-xv) + xc*xv

        return xab*(1-yv) + xdc*yv

    def plot(self, axis, u, specialPatches = None, almostSpecialPatches = None, specialNodes = None,ref = 1):
        t0 = time.time()

        for ipatch, patch in enumerate(self.patches):
            if patch!=None:
                if specialPatches!= None and ipatch in specialPatches[1]:
                # if specialPatches=="all":
                    patch.plot(axis, color = specialPatches[0], ref=ref)
                elif almostSpecialPatches!= None and ipatch in almostSpecialPatches[1]:
                # if specialPatches=="all":
                    patch.plot(axis, color = almostSpecialPatches[0], ref=ref)
                else:
                    patch.plot(axis, color = (.5,.5,.5,0.4), ref=ref)           # Inactive part of master
                    # self.plotQuad(axis, u, ipatch, color=(.5,.5,.5,0.7), ref=ref)
            else:
                self.plotQuad(axis, u, ipatch, color=(.5,.5,.5,0.6) ,ref=ref)   #Slave Body

            #     patch.plot(axis, color = "grey",ref=ref)      # <<<=== FIX this!
        self.plotNodes(axis, u, special = specialNodes)
            


        tf = time.time()
        print("Plotting " +str(self.body.name) +": Done [ "+str(tf-t0)+" s ]")
        
    def plotQuad(self, axis, u, iquad, color = "gray", surf=True,wire=False,ref=1):     #full of redundancy. FIX!!

        x = np.zeros((ref+1,ref+1),dtype = float)
        y = np.zeros((ref+1,ref+1),dtype = float)
        z = np.zeros((ref+1,ref+1),dtype = float)

        for i in range(ref+1):
            ti = i/(ref)
            for j in range(ref+1):
                tj = j/(ref)
                x[i,j], y[i,j], z[i,j] = self.QuadInterpolation(u,iquad,ti,tj)

        # set_trace()

        if surf:
            axis.plot_surface(x, y, z, color=color)
        if wire:
            axis.plot_wireframe(x, y, z, color="black")

    def plotNodes(self, axis, u,   special = None):
        if special == None: return None
        # Extremely slow. Slower than plotting surfaces
        for node in self.nodes:
            # if special != None and node in special[1]:
            #     x = np.array(u[self.body.DoFs[node]])+ np.array(self.X[node])
            #     axis.scatter(x[0],x[1],x[2], color = special[0], s = 0.25)    #redundancy between color and special (nothing too serious)
            for color in list(special.keys()):
                if special!=None and node in special[color]:
                    x = np.array(u[self.body.DoFs[node]])+ np.array(self.X[node])
                    axis.scatter(x[0],x[1],x[2], color=color, s = 0.25)    #redundancy between color and special (nothing too serious)


    def quadSize(self,iquad,u):
        # quad = np.array(self.quads[iquad] ) 
        # xa,xb,xc,xd = 
        pass

def HexaToFace(hexa):
    faces = np.array([list(hexa[0:4][::-1]),
            list(hexa[4:8]),
            [hexa[0],hexa[1],hexa[5],hexa[4]],
            [hexa[1],hexa[2],hexa[6],hexa[5]],
            [hexa[2],hexa[3],hexa[7],hexa[6]],
            [hexa[3],hexa[0],hexa[4],hexa[7]]])
    return faces
        
