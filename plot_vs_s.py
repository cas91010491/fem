import numpy as np
from matplotlib import pyplot as plt
from numpy.linalg import norm
from math import sqrt
from pdb import set_trace


def ms(u,X,C,EA):
    
    x=X+u.reshape(-1,2)
    ms=0
    for (i,j) in C:
        L0= norm(X[i]-X[j])
        L = norm(x[i]-x[j])
        ms+=1/2*EA*L0*((L-L0)/L0)**2
    return ms

def dms(u,X,C,EA,dofs,fr = None):
    x=X+u.reshape(-1,2)
    dms=np.zeros_like(x.ravel())
    for (i,j) in C:
        dofij = np.append(dofs[i],dofs[j])
        L0= norm(X[i]-X[j])
        L = norm(x[i]-x[j])
        dLdu = np.array([ x[i,0]-x[j,0] , x[i,1]-x[j,1] , -x[i,0]+x[j,0], -x[i,1]+x[j,1] ])/L
        dms[np.ix_(dofij)]+=EA*((L-L0)/L0)*dLdu

    # di = np.delete(np.arange(ndofs), fr)       # free DOFs 'b'
    # dms[np.ix_(di)]=0.0

    return dms if fr is None else dms[np.ix_(fr)]

def plotMS(u,X,C, xp,yc, ax= None, undef = False):
    u=u.reshape(-1,2)
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)

    x = X+u
    for (i,j) in C:
        xa, xb = x[i],x[j]
        ax.plot([xa[0],xb[0]],[xa[1],xb[1]], c="blue")
        if undef:
            Xa, Xb = X[i],X[j]
            ax.plot([Xa[0],Xb[0]],[Xa[1],Xb[1]],lw=0.5, c="grey")

    ax.plot([-2,3],[yc,yc],c="red")
    ax.scatter(xp,u.ravel()[3], marker="X")
    ax.set_xlim([-2,3])
    ax.set_ylim([-1.5,0.5])


    if ax is None: plt.show()


def minimize_u(u0,fr):

    alpha_init = 1
    c, r = 0.5, 0.5
    nreset = 10

    u = np.array(u0)
    gn = min(yc-u[3],0.0)
    mN = 1/2*kn*gn**2
    fN = kn*gn
    mnew = ms(u,X,C,EA) + mN
    alpha = alpha_init
    mold = 10**100
    fnew,hnew, cnt = np.zeros(ndofs), 0.0, 0

    while abs(mnew-mold)>tol:
        mold=mnew
        fold = fnew
        hold = hnew

        fnew[np.ix_(fr)] = dms(u,X,C,EA,dofs,fr)
        gn = min(yc-u[3],0.0)
        fN = -kn*gn
        fnew[3] += fN

        if cnt%nreset==0:
            hnew=-fnew
        else:
            beta = (fnew@fnew)/(fold@fold)
            # beta = (fnew@(fnew-fold))/(fold@fold)
            # beta = (fnew@(fnew-fold))/(hold@(fnew-fold))
            # beta = (fnew@fnew)/(hold@(fnew-fold))
            hnew = -fnew+max(0,beta)*hold
        
        alpha = (1/r)*alpha_init
        mx = 10**100
        xniter= 0 
        ux = np.array(u)
        while mx>mnew+c*alpha*(hnew@fnew):
            alpha = r*alpha
            ux= u + (alpha*hnew)
            
            gn = min(yc-ux[3],0.0)
            mN = 1/2*kn*gn**2
            mx = ms(ux,X,C,EA) + mN
            xniter +=1

        mnew = mx
        u = ux
        cnt += 1

    # fprint = dms(u,X,C,EA,dofs,fr)
    # fprint += fN
    # print("\t\tdone optimizing u =",u,"\tniter:",cnt,"\txniter:",xniter,"\tFPRINT =", norm(fprint) )

    return u



X = np.array([[-1.0,-1.0],[0.0,0.0],[1.0,-1.0]])
dofs = np.array([[0,1],[2,3],[4,5]])
C = [[0,1],[1,2]]
EA = 1.0

xp, yc = 0.05, -0.1        # hook(x)  and wall(y) position


ndofs = len(dofs.ravel())

kn = 100
mu= 100000
tol = 1e-20

fr = [3]
u = np.zeros(6)
S = np.linspace(0.0,0.06,1001)
GN, GT2 , dFT, MUFN, ABSFT, OBJ= [], [], [], [], [], []
for xs in S:
    u[2] = xs
    u = minimize_u(u,fr)
    gn = yc - u[3]
    Ft, Fn = dms(u,X,C,EA,dofs)[2:4]
    Fn = -kn*gn
    gt = u[2]-xp
    Tau = -gt/abs(gt)
    dft = abs(Ft)- Ft*Tau
    muFn = mu*Fn
    pen = 1e5
    obj = gt**2 + 0.5*pen*(Ft*Tau - abs(Ft))**2 + 0.5*pen*(mu*Fn-abs(Ft))**2*(Ft>mu*Fn)

    GN.append(100*gn)
    GT2.append(gt**2)
    dFT.append(dft)
    MUFN.append(muFn)
    ABSFT. append(abs(Ft))
    OBJ.append(obj)

    # plotMS(u,X,C,xp,yc, ax= None, undef = True)
    # plt.show()
    

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax1. plot(S,GN   ,label="100x gn")
ax1. plot(S,GT2  ,label="gt2")
ax2. plot(S,dFT  ,label="τ·Ft-|Ft|")
ax3. plot(S,MUFN ,label="μFn")
ax3. plot(S,ABSFT,label="|Ft|")

ax1.set_xlabel("s"), ax1.legend()
ax2.set_xlabel("s"), ax2.legend()
ax3.set_xlabel("s"), ax3.legend()

fi = plt.figure()
ax = fi.add_subplot(111)
ax. plot(S,OBJ,label="Obj")

plt.show()


