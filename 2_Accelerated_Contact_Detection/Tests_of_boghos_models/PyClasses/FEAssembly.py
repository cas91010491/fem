from PyClasses.Utilities import flatList
from PyClasses.Surfaces import *

import numpy as np
import pickle, time
import concurrent.futures
import multiprocessing
from scipy import sparse
from numpy.linalg import norm as norm
from pdb import set_trace
from math import sqrt, pi, cos, sin
import random

from PyClasses import eigen_backend


class FEAssembly:
    def __init__(self,X,hexahedrons , name = "unnamed", recOuters = False):
        self.X = X
        self.hexas = hexahedrons.tolist() if type(hexahedrons)!= list else hexahedrons
        self.deleteReduntants()
        self.name = name
        self.DoFs = []              # assigned by FEModel constructor
        outers = pickle.load(open("RecoveryOuters.dat","rb"))[name] if recOuters else None
        self.surf = Surface(self, outers = outers)
        self.Youngsmodulus = 1
        self.Poissonsratio = 0.3
        self.isRigid = False

    def deleteReduntants(self):
        reds = self.redundantNodes()
        self.X = np.delete(self.X,reds,0)
        self.ReduceHexaOrders(reds)

    def redundantNodes(self):
        reds = []
        nodesinhexas = list(set(flatList(self.hexas)))
        for i in range(len(self.X)):
            if i not in nodesinhexas:
                reds.append(i)
        return reds

    def ReduceHexaOrders(self,idxs):
        sortedIdxs = sorted(idxs, reverse = True)
        new_hexas = []
        for hexa in self.hexas:
            new_hexa = []
            for node in hexa:
                for idx in sortedIdxs:
                    if node>idx:
                        node -= 1
                new_hexa.append(node)
            new_hexas.append(new_hexa)
        
        self.hexas = new_hexas

    # TRANSFORMATIONS
    def Translate(self, disp):
        newX = np.zeros_like(self.X)
        disp = np.array(disp)
        for i in range(len(self.X)):
            newX[i] = self.X[i] + disp

        self.X = newX
        self.surf.X = newX

    def Rotate(self, theta,dir, unit="rad"):
        import math as m

        if unit == "deg": theta=(pi/180)*theta       
        if dir in "xX":
            M = np.matrix([[ 1, 0           , 0           ],
                           [ 0, m.cos(theta),-m.sin(theta)],
                           [ 0, m.sin(theta), m.cos(theta)]])
        
        elif dir in "yY":
            M = np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                           [ 0           , 1, 0           ],
                           [-m.sin(theta), 0, m.cos(theta)]])

        else:
            M = np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                           [ m.sin(theta), m.cos(theta) , 0 ],
                           [ 0           , 0            , 1 ]])
        newX = np.array(self.X@M)
        self.X = newX
        self.surf.X = newX

    def RandDistort(self,d , nodes = None):
        """
            Distorts the position of the nodes in the provided list (local id for the current FEAssembly) an amount r_i such that 0<r_i<d
        """
        newX = np.array(self.X)
        if nodes == None:
            nodes = range(len(self.X))
        for node in nodes:
            H_ang = random.uniform(0,2*pi)
            V_ang = random.uniform(0,pi)
            r = random.uniform(0,d)
            newX[node] += np.array([ r*sin(V_ang)*cos(H_ang) , r*sin(V_ang)*sin(H_ang) , r*cos(V_ang) ])

            self.X = newX
            self.surf.X = newX

    def Resize(self, factor, center = None, dir = None):
        xc = np.array([0,0,0]) if center is None else np.array(center)
        newX = np.zeros_like(self.X)

        if dir is None:
            for idx, x in enumerate(self.X):
                newX[idx] = xc + factor*(x-xc)
        else:
            if dir in "Xx":
                dir = 0
            elif dir in "Yy":
                dir = 1
            elif dir in "Zz":
                dir = 2
            newX = self.X
            for idx, x in enumerate(self.X):
                newX[idx][dir] = xc[dir] + factor*(x[dir]-xc[dir])


        self.X = newX
        self.surf.X = newX


    def Bulge(self, factor, center=None, dir=None):
        xc = np.array([0, 0, 0]) if center is None else np.array(center)
        newX = np.zeros_like(self.X)

        if dir is None:
            for idx, x in enumerate(self.X):
                dist = np.linalg.norm(x - xc)
                newX[idx] = x + factor * (x - xc) / (dist + 1e-8)
        else:
            if dir in "Xx":
                dir = 0
            elif dir in "Yy":
                dir = 1
            elif dir in "Zz":
                dir = 2
            newX = self.X.copy()  # Create a copy of the original coordinates
            for idx, x in enumerate(self.X):
                dist = np.abs(x[dir] - xc[dir])
                if dist > 1e-8:
                    newX[idx][dir] = x[dir] + factor * (x[dir] - xc[dir]) / dist

        self.X = newX
        self.surf.X = newX



    # SELECTION OF NODES
    def SelectNodesByBox(self,xyz_a,xyz_b, OnSurface = False):
        nodes = []
        xa, ya, za = xyz_a
        xb, yb, zb = xyz_b
        for node_id, [x,y,z] in enumerate(self.X):
            if OnSurface and node_id not in self.surf.nodes:
                continue
            if xa<=x<=xb:
                if ya<=y<=yb:
                    if za<=z<=zb:
                        nodes.append(node_id)
        return nodes

    def SelectFlatSide(self, side, tol = 1e-6, OnSurface=False):
        minX = min(self.X[:,0])
        maxX = max(self.X[:,0])
        minY = min(self.X[:,1])
        maxY = max(self.X[:,1])
        minZ = min(self.X[:,2])
        maxZ = max(self.X[:,2])
        eps = tol
        if "x" in side:
            if "-" in side:
                return self.SelectNodesByBox([minX-eps,minY-eps,minZ-eps],[minX+eps,maxY+eps,maxZ+eps],OnSurface = OnSurface)
            else:
                return self.SelectNodesByBox([maxX-eps,minY-eps,minZ-eps],[maxX+eps,maxY+eps,maxZ+eps],OnSurface = OnSurface)
        elif "y" in side:
            if "-" in side:
                return self.SelectNodesByBox([minX-eps,minY-eps,minZ-eps],[maxX+eps,minY+eps,maxZ+eps],OnSurface = OnSurface)
            else:
                return self.SelectNodesByBox([minX-eps,maxY-eps,minZ-eps],[maxX+eps,maxY+eps,maxZ+eps],OnSurface = OnSurface)
        else:
            if "-" in side:
                return self.SelectNodesByBox([minX-eps,minY-eps,minZ-eps],[maxX+eps,maxY+eps,minZ+eps],OnSurface = OnSurface)
            else:
                return self.SelectNodesByBox([minX-eps,minY-eps,maxZ-eps],[maxX+eps,maxY+eps,maxZ+eps],OnSurface = OnSurface)

    def SelectHigherThan(self, dir, val = 0, Strict = True,OnSurface = False):
        nodes = []
        if type(dir) == int:
            None
        elif dir in "Xx":
            dir = 0
        elif dir in "Yy":
            dir = 1
        elif dir in "Zz":
            dir = 2
        for node_id, xyz in enumerate(self.X):
            if OnSurface and node_id not in self.surf.nodes:
                continue
            if Strict:
                if xyz[dir] > val:
                    nodes.append(node_id)
            elif xyz[dir] >= val:
                    nodes.append(node_id)


        return nodes

    def SelectLowerThan(self, dir, val = 0, Strict = True, OnSurface = False):
        nodes = []
        if type(dir) == int:
            None
        elif dir in "Xx":
            dir = 0
        elif dir in "Yy":
            dir = 1
        elif dir in "Zz":
            dir = 2
        for node_id, [x,y,z] in enumerate(self.X):
            if OnSurface and node_id not in self.surf.nodes:
                continue
            if Strict:
                if [x,y,z][dir] < val:
                    nodes.append(node_id)
            elif [x,y,z][dir] <= val:
                    nodes.append(node_id)
        return nodes

    def SelectNodesBySphere(self,xyz_center,radius, OnSurface = False):
        nodes = []
        xc = np.array(xyz_center,dtype=float)
        r = radius
        for node_id, x_node in enumerate(self.X):
            if OnSurface and node_id not in self.surf.nodes:
                continue
            x = np.array(x_node)
            if sqrt((x-xc)@(x-xc))<=r:
                nodes.append(node_id)
        return nodes

    def SelectAll(self):
        return range(len(self.X))

    # SELECTION OF QUADS:
    def SelectQuadsByNodes(self,nodes):
        quads = []
        nodeset = set(nodes)  # no surface needed. Hexas will be taken anyway
        for qid, quad in enumerate(self.surf.quads):
            common = set(quad).intersection(nodeset)
            if len(common)>0:
                quads.append(qid)
        return quads



    # SELECTION OF HEXAS:
    def SelectHexasByBox(self,xyz_a,xyz_b):
        nodes = self.SelectNodesByBox(xyz_a,xyz_b)  # no surface needed. Hexas will be taken anyway
        return self.SelectHexasByNodes(nodes)

    def SelectHexasByNodes(self,nodes):
        hexas = []
        nodeset = set(nodes)  # no surface needed. Hexas will be taken anyway
        for hid, hexa in enumerate(self.hexas):
            common = set(hexa).intersection(nodeset)
            if len(common)>0:
                hexas.append(hid)
        return hexas



    # MECHANICAL FUNCTIONS
    # Residuals at Finite Element level
    # 1. using cpp
    def fint_el(self,hexa,u):
        X_el = self.X[hexa]
        u_el = u[self.DoFs[hexa]]
        return eigen_backend.f_element_py(X_el, u_el, self.Youngsmodulus, self.Poissonsratio).ravel()
  
    def K_el(self,hexa,u):
        X_el = self.X[hexa]
        u_el = u[self.DoFs[hexa]]
        return eigen_backend.K_element_py(X_el, u_el, self.Youngsmodulus, self.Poissonsratio)

    # 2. using python (slower)
    def m_el(self,hexa,u):
        m=0
        X = self.X[hexa]
        u = u[self.DoFs[hexa]]

        d1 = self.Youngsmodulus*self.Poissonsratio/(2*(1+self.Poissonsratio)*(1-2*self.Poissonsratio))
        c10 = self.Youngsmodulus/(4*(1+self.Poissonsratio))

        gauss_points = 1/sqrt(3)*np.array([ [-1, -1, -1],
                                            [ 1, -1, -1],
                                            [ 1,  1, -1],
                                            [-1,  1, -1],
                                            [-1, -1,  1],
                                            [ 1, -1,  1],
                                            [ 1,  1,  1],
                                            [-1,  1,  1]])

        for (g1,g2,g3) in gauss_points:
            dNd_xi = 1/8*np.array([ [-(1-g2)*(1-g3) , -(1-g1)*(1-g3) , -(1-g1)*(1-g2) ],
                                    [ (1-g2)*(1-g3) , -(1+g1)*(1-g3) , -(1+g1)*(1-g2) ],
                                    [ (1+g2)*(1-g3) ,  (1+g1)*(1-g3) , -(1+g1)*(1+g2) ],
                                    [-(1+g2)*(1-g3) ,  (1-g1)*(1-g3) , -(1-g1)*(1+g2) ],
                                    [-(1-g2)*(1+g3) , -(1-g1)*(1+g3) ,  (1-g1)*(1-g2) ],
                                    [ (1-g2)*(1+g3) , -(1+g1)*(1+g3) ,  (1+g1)*(1-g2) ],
                                    [ (1+g2)*(1+g3) ,  (1+g1)*(1+g3) ,  (1+g1)*(1+g2) ],
                                    [-(1+g2)*(1+g3) ,  (1-g1)*(1+g3) ,  (1-g1)*(1+g2) ]])

            J = np.dot(dNd_xi.T,X)
            invJ = np.linalg.inv(J)
            detJ = np.linalg.det(J)

            dNdx = np.dot( dNd_xi , invJ.T )
            F = np.eye(len(dNdx.T)) + np.dot(dNdx.T,u).T
            detF = np.linalg.det(F)


            SED=c10*(np.trace(F.T@F)-3-2*np.log(detF))+d1*(np.log(detF))**2

            m += SED*detJ

        return m
    


    def m_el_extra(self, hexa, u):
        m = np.zeros(len(self.DoFs[hexa]))  # Initialize array to store nodal SED values
        X = self.X[hexa]
        u = u[self.DoFs[hexa]]

        d1 = self.Youngsmodulus * self.Poissonsratio / (2 * (1 + self.Poissonsratio) * (1 - 2 * self.Poissonsratio))
        c10 = self.Youngsmodulus / (4 * (1 + self.Poissonsratio))

        gauss_points = 1 / np.sqrt(3) * np.array([[-1, -1, -1],
                                                [1, -1, -1],
                                                [1, 1, -1],
                                                [-1, 1, -1],
                                                [-1, -1, 1],
                                                [1, -1, 1],
                                                [1, 1, 1],
                                                [-1, 1, 1]])

        for (g1, g2, g3) in gauss_points:
            dNd_xi = 1 / 8 * np.array([[-(1 - g2) * (1 - g3), -(1 - g1) * (1 - g3), -(1 - g1) * (1 - g2)],
                                    [(1 - g2) * (1 - g3), -(1 + g1) * (1 - g3), -(1 + g1) * (1 - g2)],
                                    [(1 + g2) * (1 - g3), (1 + g1) * (1 - g3), -(1 + g1) * (1 + g2)],
                                    [-(1 + g2) * (1 - g3), (1 - g1) * (1 - g3), -(1 - g1) * (1 + g2)],
                                    [-(1 - g2) * (1 + g3), -(1 - g1) * (1 + g3), (1 - g1) * (1 - g2)],
                                    [(1 - g2) * (1 + g3), -(1 + g1) * (1 + g3), (1 + g1) * (1 - g2)],
                                    [(1 + g2) * (1 + g3), (1 + g1) * (1 + g3), (1 + g1) * (1 + g2)],
                                    [-(1 + g2) * (1 + g3), (1 - g1) * (1 + g3), (1 - g1) * (1 + g2)]])

            J = np.dot(dNd_xi.T, X)
            invJ = np.linalg.inv(J)
            detJ = np.linalg.det(J)

            dNdx = np.dot(dNd_xi, invJ.T)
            F = np.eye(len(dNdx.T)) + np.dot(dNdx.T, u).T
            detF = np.linalg.det(F)

            SED = c10 * (np.trace(F.T @ F) - 3 - 2 * np.log(detF)) + d1 * (np.log(detF))**2

            # Extrapolate the SED values to the nodes
            N = 1 / 8 * np.array([(1 - g1) * (1 - g2) * (1 - g3),
                                (1 + g1) * (1 - g2) * (1 - g3),
                                (1 + g1) * (1 + g2) * (1 - g3),
                                (1 - g1) * (1 + g2) * (1 - g3),
                                (1 - g1) * (1 - g2) * (1 + g3),
                                (1 + g1) * (1 - g2) * (1 + g3),
                                (1 + g1) * (1 + g2) * (1 + g3),
                                (1 - g1) * (1 + g2) * (1 + g3)])

            m += SED * detJ * N

        return m


    def get_nodal_SED(self,u):
        if self.isRigid: return None
        sed_tot = np.zeros((len(self.X),2))

        for hexa in self.hexas:
            sed_tot[hexa,0] *= sed_tot[hexa,1]
            sed_tot[hexa,1] += 1
            sed_tot[hexa,0] += self.m_el_extra(hexa, u)
            sed_tot[hexa,0] /= sed_tot[hexa,1]

        return sed_tot[:,0]

    def fint_el_classic(self,hexa,u):
        fint_el_tens = np.zeros((8,3),dtype=float)
        X = self.X[hexa]
        u = u[self.DoFs[hexa]]

        d1 = self.Youngsmodulus*self.Poissonsratio/(2*(1+self.Poissonsratio)*(1-2*self.Poissonsratio))
        c10 = self.Youngsmodulus/(4*(1+self.Poissonsratio))

        gauss_points = 1/sqrt(3)*np.array([ [-1, -1, -1],
                                            [ 1, -1, -1],
                                            [ 1,  1, -1],
                                            [-1,  1, -1],
                                            [-1, -1,  1],
                                            [ 1, -1,  1],
                                            [ 1,  1,  1],
                                            [-1,  1,  1]])

        for (g1,g2,g3) in gauss_points:
            dNd_xi = 1/8*np.array([ [-(1-g2)*(1-g3) , -(1-g1)*(1-g3) , -(1-g1)*(1-g2) ],
                                    [ (1-g2)*(1-g3) , -(1+g1)*(1-g3) , -(1+g1)*(1-g2) ],
                                    [ (1+g2)*(1-g3) ,  (1+g1)*(1-g3) , -(1+g1)*(1+g2) ],
                                    [-(1+g2)*(1-g3) ,  (1-g1)*(1-g3) , -(1-g1)*(1+g2) ],
                                    [-(1-g2)*(1+g3) , -(1-g1)*(1+g3) ,  (1-g1)*(1-g2) ],
                                    [ (1-g2)*(1+g3) , -(1+g1)*(1+g3) ,  (1+g1)*(1-g2) ],
                                    [ (1+g2)*(1+g3) ,  (1+g1)*(1+g3) ,  (1+g1)*(1+g2) ],
                                    [-(1+g2)*(1+g3) ,  (1-g1)*(1+g3) ,  (1-g1)*(1+g2) ]])

            J = np.dot(dNd_xi.T,X)
            invJ = np.linalg.inv(J)
            detJ = np.linalg.det(J)

            dNdx = np.dot( dNd_xi , invJ.T )
            F = np.eye(len(dNdx.T)) + np.dot(dNdx.T,u).T
            detF = np.linalg.det(F)
            Finvc = np.linalg.inv(F).T

            P = 2*c10*(F-Finvc) + 2*d1*Finvc*np.log(detF)

            fint_el_tens += np.dot(dNdx,P.T)*detJ

        return np.reshape(fint_el_tens,24)
        
    def K_el_classic(self,hexa,u):
        K_el = np.zeros((8,3,8,3))
        X = self.X[hexa]
        u = u[self.DoFs[hexa]]
        Cb = self.Youngsmodulus*self.Poissonsratio/((1+self.Poissonsratio)*(1-2*self.Poissonsratio))
        Ca = self.Youngsmodulus/(2*(1+self.Poissonsratio))

        gauss_points = 1/sqrt(3)*np.array([ [-1, -1, -1],
                                            [ 1, -1, -1],
                                            [ 1,  1, -1],
                                            [-1,  1, -1],
                                            [-1, -1,  1],
                                            [ 1, -1,  1],
                                            [ 1,  1,  1],
                                            [-1,  1,  1]])

        for (g1,g2,g3) in gauss_points:
            dNd_xi = 1/8*np.array( [[-(1-g2)*(1-g3) , -(1-g1)*(1-g3) , -(1-g1)*(1-g2) ],
                                    [ (1-g2)*(1-g3) , -(1+g1)*(1-g3) , -(1+g1)*(1-g2) ],
                                    [ (1+g2)*(1-g3) ,  (1+g1)*(1-g3) , -(1+g1)*(1+g2) ],
                                    [-(1+g2)*(1-g3) ,  (1-g1)*(1-g3) , -(1-g1)*(1+g2) ],
                                    [-(1-g2)*(1+g3) , -(1-g1)*(1+g3) ,  (1-g1)*(1-g2) ],
                                    [ (1-g2)*(1+g3) , -(1+g1)*(1+g3) ,  (1+g1)*(1-g2) ],
                                    [ (1+g2)*(1+g3) ,  (1+g1)*(1+g3) ,  (1+g1)*(1+g2) ],
                                    [-(1+g2)*(1+g3) ,  (1-g1)*(1+g3) ,  (1-g1)*(1+g2) ]])

            J = np.dot(dNd_xi.T,X)
            invJ = np.linalg.inv(J)
            detJ = np.linalg.det(J)

            dNdx = np.dot( dNd_xi , invJ.T )        # <<==== J goes transposed!!!

            F = np.eye(len(dNdx.T)) + np.dot(dNdx.T,u).T
            detF = np.linalg.det(F)
            F11 = F[0,0]
            F12 = F[0,1]
            F13 = F[0,2]
            F21 = F[1,0]
            F22 = F[1,1]
            F23 = F[1,2]
            F31 = F[2,0]
            F32 = F[2,1]
            F33 = F[2,2]
            Finv11= (F22*F33 - F23*F32)/detF
            Finv12=-(F12*F33 - F13*F32)/detF
            Finv13= (F12*F23 - F13*F22)/detF
            Finv21=-(F21*F33 - F23*F31)/detF
            Finv22= (F11*F33 - F13*F31)/detF
            Finv23=-(F11*F23 - F13*F21)/detF
            Finv31= (F21*F32 - F22*F31)/detF
            Finv32=-(F11*F32 - F12*F31)/detF
            Finv33= (F11*F22 - F12*F21)/detF

            LogDetF = np.log(detF)

            C1111=Ca + Cb*Finv11**2 + Finv11**2*(Ca - Cb*LogDetF)
            C1112=Finv11*Finv21*(Ca + Cb - Cb*LogDetF)
            C1113=Finv11*Finv31*(Ca + Cb - Cb*LogDetF)
            C1121=Finv11*Finv12*(Ca + Cb - Cb*LogDetF)
            C1122=Cb*Finv11*Finv22 + Finv12*Finv21*(Ca - Cb*LogDetF)
            C1123=Cb*Finv11*Finv32 + Finv12*Finv31*(Ca - Cb*LogDetF)
            C1131=Finv11*Finv13*(Ca + Cb - Cb*LogDetF)
            C1132=Cb*Finv11*Finv23 + Finv13*Finv21*(Ca - Cb*LogDetF)
            C1133=Cb*Finv11*Finv33 + Finv13*Finv31*(Ca - Cb*LogDetF)
            C1212=Cb*Finv12*Finv21 + Finv11*Finv22*(Ca - Cb*LogDetF)
            C1213=Cb*Finv12*Finv31 + Finv11*Finv32*(Ca - Cb*LogDetF)
            C1221=Ca + Cb*Finv12**2 + Finv12**2*(Ca - Cb*LogDetF)
            C1222=Finv12*Finv22*(Ca + Cb - Cb*LogDetF)
            C1223=Finv12*Finv32*(Ca + Cb - Cb*LogDetF)
            C1231=Finv12*Finv13*(Ca + Cb - Cb*LogDetF)
            C1232=Cb*Finv12*Finv23 + Finv13*Finv22*(Ca - Cb*LogDetF)
            C1233=Cb*Finv12*Finv33 + Finv13*Finv32*(Ca - Cb*LogDetF)
            C1312=Cb*Finv13*Finv21 + Finv11*Finv23*(Ca - Cb*LogDetF)
            C1313=Cb*Finv13*Finv31 + Finv11*Finv33*(Ca - Cb*LogDetF)
            C1322=Cb*Finv13*Finv22 + Finv12*Finv23*(Ca - Cb*LogDetF)
            C1323=Cb*Finv13*Finv32 + Finv12*Finv33*(Ca - Cb*LogDetF)
            C1331=Ca + Cb*Finv13**2 + Finv13**2*(Ca - Cb*LogDetF)
            C1332=Finv13*Finv23*(Ca + Cb - Cb*LogDetF)
            C1333=Finv13*Finv33*(Ca + Cb - Cb*LogDetF)
            C2112=Ca + Cb*Finv21**2 + Finv21**2*(Ca - Cb*LogDetF)
            C2113=Finv21*Finv31*(Ca + Cb - Cb*LogDetF)
            C2122=Finv21*Finv22*(Ca + Cb - Cb*LogDetF)
            C2123=Cb*Finv21*Finv32 + Finv22*Finv31*(Ca - Cb*LogDetF)
            C2132=Finv21*Finv23*(Ca + Cb - Cb*LogDetF)
            C2133=Cb*Finv21*Finv33 + Finv23*Finv31*(Ca - Cb*LogDetF)
            C2213=Cb*Finv22*Finv31 + Finv21*Finv32*(Ca - Cb*LogDetF)
            C2222=Ca + Cb*Finv22**2 + Finv22**2*(Ca - Cb*LogDetF)
            C2223=Finv22*Finv32*(Ca + Cb - Cb*LogDetF)
            C2232=Finv22*Finv23*(Ca + Cb - Cb*LogDetF)
            C2233=Cb*Finv22*Finv33 + Finv23*Finv32*(Ca - Cb*LogDetF)
            C2313=Cb*Finv23*Finv31 + Finv21*Finv33*(Ca - Cb*LogDetF)
            C2323=Cb*Finv23*Finv32 + Finv22*Finv33*(Ca - Cb*LogDetF)
            C2332=Ca + Cb*Finv23**2 + Finv23**2*(Ca - Cb*LogDetF)
            C2333=Finv23*Finv33*(Ca + Cb - Cb*LogDetF)
            C3113=Ca + Cb*Finv31**2 + Finv31**2*(Ca - Cb*LogDetF)
            C3123=Finv31*Finv32*(Ca + Cb - Cb*LogDetF)
            C3133=Finv31*Finv33*(Ca + Cb - Cb*LogDetF)
            C3223=Ca + Cb*Finv32**2 + Finv32**2*(Ca - Cb*LogDetF)
            C3233=Finv32*Finv33*(Ca + Cb - Cb*LogDetF)
            C3333=Ca + Cb*Finv33**2 + Finv33**2*(Ca - Cb*LogDetF)

            C11 = np.array([[C1111,C1112,C1113],
                            [C1121,C1122,C1123],
                            [C1131,C1132,C1133]])
            C12 = np.array([[C1121,C1212,C1213],
                            [C1221,C1222,C1223],
                            [C1231,C1232,C1233]])
            C13 = np.array([[C1131,C1312,C1313],
                            [C1231,C1322,C1323],
                            [C1331,C1332,C1333]])
            C21 = np.array([[C1112,C2112,C2113],
                            [C1212,C2122,C2123],
                            [C1312,C2132,C2133]])
            C22 = np.array([[C1122,C2122,C2213],
                            [C1222,C2222,C2223],
                            [C1322,C2232,C2233]])
            C23 = np.array([[C1132,C2132,C2313],
                            [C1232,C2232,C2323],
                            [C1332,C2332,C2333]])
            C31 = np.array([[C1113,C2113,C3113],
                            [C1213,C2213,C3123],
                            [C1313,C2313,C3133]])
            C32 = np.array([[C1123,C2123,C3123],
                            [C1223,C2223,C3223],
                            [C1323,C2323,C3233]])
            C33 = np.array([[C1133,C2133,C3133],
                            [C1233,C2233,C3233],
                            [C1333,C2333,C3333]])

            PdF = np.array([[C11,C12,C13],
                            [C21,C22,C23],
                            [C31,C32,C33]])

            Kgp = np.tensordot(dNdx,np.tensordot(PdF,dNdx,axes=[3,1]),axes=[1,0]).swapaxes(2,3)
            # Kgp = np.tensordot(dNdx,np.tensordot(PdF.swapaxes(0,1),dNdx.T,axes=1),axes=1).swapaxes(2,3)



            K_el+=Kgp*detJ

        return K_el.reshape(24,24)

    # Residuals at FEAssembly level
    # def get_K(self, Model,NowDiff=False):

    #     return eigen_backend.K_global_alt(Model.X,Model.u,self.hexas,self.DoFs,self.Youngsmodulus,self.Poissonsratio)


    def get_K_Cem(self, Model,NowDiff=False):
        if self.isRigid:
            return None

        if not hasattr(self,"SparseRow"):
            RC = eigen_backend.K_global_cr(Model.X,Model.u_temp,self.hexas,self.DoFs,self.Youngsmodulus,self.Poissonsratio)
            vals = eigen_backend.K_global_cr(Model.X,Model.u_temp,self.hexas,self.DoFs,self.Youngsmodulus,self.Poissonsratio)
            self.SparseRow = RC.rows
            self.SparseCol = RC.cols

        vals = eigen_backend.K_global_val(Model.X,Model.u_temp,self.hexas,self.DoFs,self.Youngsmodulus,self.Poissonsratio)
        Ksparse = sparse.csr_matrix((vals,(self.SparseRow,self.SparseCol)),shape=Model.K.shape)

        Model.K += Ksparse


    def get_K(self, Model,NowDiff=False):
        if self.isRigid:
            return None

        computeRC = False
        if not hasattr(self,"SparseRow"):
            self.SparseRow = np.empty(len(self.hexas)*576,dtype=int  )
            self.SparseCol = np.empty(len(self.hexas)*576,dtype=int  )
            computeRC = True

        V = np.empty(len(self.hexas)*576,dtype=float)
        for i,hexa in enumerate(self.hexas):
            # V[i*576:(i+1)*576] = self.K_el(hexa, Model.u_temp).ravel()
            V[i*576:(i+1)*576] = self.K_el_classic(hexa, Model.u_temp).ravel()
            if computeRC:
                dofs = self.DoFs[hexa].ravel()
                self.SparseRow[i*576:(i+1)*576] = np.repeat(dofs,len(dofs))
                self.SparseCol[i*576:(i+1)*576] = np.  tile(dofs,len(dofs))
        Kbody = sparse.csr_matrix((V,(self.SparseRow,self.SparseCol)),shape=Model.K.shape)

        if NowDiff:
            return Kbody
        else:
            Model.K += Kbody

    def get_fint(self, Model,temp=False):
        if self.isRigid:
            return None

        u = Model.u_temp if temp else Model.u
        for hexa in self.hexas:
            fint_el = self.fint_el(hexa,u)
            dofs = self.DoFs[hexa].ravel()
            Model.fint[np.ix_(dofs)] += fint_el


    # functions for BFGS
    def compute_m(self,u):
        if self.isRigid:
            return 0

        m = 0
        for hexa in self.hexas:
            m += self.m_el(hexa,u)

        return m
    
    def compute_f(self,u, Model):
        force=np.zeros(Model.fint.shape)
        if self.isRigid:
            return 0

        for hexa in self.hexas:
            fint_el = self.fint_el(hexa,u)
            dofs = self.DoFs[hexa].ravel()
            force[dofs] += fint_el
        return force


    # VISUAL
    def plot(self,ax = None, u = None,plotWhat = "surf",nodes="all",plotcoords = False, sed = None,ref=10):

        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        if plotcoords: 
            from PyClasses.Utilities import plot_coords
            plot_coords(ax)
        if plotWhat in ["nodes","all"]:
            if u is not None:
                X = self.X + u[:3*len(self.X)].reshape(-1,3)+ np.array([0,0,-1e-2])
            else:
                X = self.X
            for idx,x in enumerate(X):
                if nodes=="all" or idx in nodes:
                    ax.scatter(x[0],x[1],x[2], c="blue", s = 0.2)
        if plotWhat in ["surf","all"]:
            if u is None:
                u = np.zeros_like(self.X)
            surfObj = self.surf.plot(ax,u, sed=sed,ref=ref)

        return surfObj


def K_hexa_glob(self,Model,hexa):
    dofs = self.DoFs[hexa].ravel()
    vals = self.K_el(hexa,Model.u_temp).ravel()     # ,-- "original"
    r = np.repeat(dofs,len(dofs))
    c = np.  tile(dofs,len(dofs))
    return vals, r ,c

