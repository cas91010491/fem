import numpy as np
from numpy import sign
from numpy.linalg import norm, inv
from pdb import set_trace



def HH(vect,i=0):
    """
    returns HouseHolder transformation matrix
    """
    ei = np.zeros_like(vect)
    ei[i] = 1
    sigma = sign(vect[0])
    u = vect + sigma*norm(vect)*ei
    Q = np.eye(len(vect))-2*np.outer(u,u)/np.inner(u,u)
    return Q

def NNLS(E,f,x,w,z,P,Z):
    m2,n = E.shape
    P=[]
    Z=range(n)
    while True:
        x=np.zeros_like(f)
        w = E.T@(f-E@x)
        if len(Z)==0 or np.all(w,w<0):
            return x
        P.append(Z.pop(np.argmax(w)))   # moves the index with maximum w from Z to P
        Ep = np.zeros_like(E)
        for j in range(n):
            Ep[j] = E[j] if j in P else np.zeros(m2)
            z = np.zeros_like(x)
            z[P] = (inv(Ep.T@Ep)@Ep.T@f)[P]
        if np.all(z[P],z[P]>0):
            x = z.copy()
            continue
        

    



C = np.array([[0.4087,0.1593]])
E = np.array([[0.4302,0.3516],[0.6246,0.3384]])
d = 0.1376
f = np.array([[0.6593],[0.9666]])

K = HH(C,0)
C1 = (C@K)[0]
E1,E2 = (E@K).T
y1 = d/C1
fbar = f-E1*y1
y1 = 0

set_trace()