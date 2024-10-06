from pdb import set_trace
from re import X
import numpy as np
from numpy.linalg import norm
from math import sqrt

def m1(u, g, a, kn, kt, mu, xp, yc, cub=None, showGN=False):
    a1, a2 = a
    g1, g2 = g

    # set_trace()


    mint = a1*(sqrt((1+u[0])**2+(1+u[1])**2)-sqrt(2))**2 + a2*(sqrt((1-u[0])**2+(1+u[1])**2)-sqrt(2))**2 - g1*u[0] - g2*u[1]
    gn = min(yc-u[1],0)
    
    if showGN: print(gn)

    if cub is None:
        mN = 0.5*kn*gn**2
    else:
        mN = kn/2*(gn**2-gn*cub+(cub**2)/3) if gn<cub else kn/(6*cub)*gn**3
    gt = abs(xp-u[0]) if gn<0 else 0
    gt_c = mu*(kn/kt)*abs(gn)
    elastic = kt*gt<=mu*abs(kn*gn)
    # if elastic: set_trace()
    mT = 0.5*kt*gt**2 if elastic else 0.5*kt*gt_c**2+mu*kn*abs(gn)*(gt-gt_c)
    # return mT, elastic 
    return mint + mN + mT, elastic 

def dm1(u, g, a, kn, kt, mu, xp, yc, cub=None, showGN=False):
    fint = np.array([2*a[0]*(u[0] + 1)*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((u[0] + 1)**2 + (u[1] + 1)**2) + 2*a[1]*(u[0] - 1)*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((1 - u[0])**2 + (u[1] + 1)**2) - g[0],
            2*a[0]*(u[1] + 1)*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((u[0] + 1)**2 + (u[1] + 1)**2) + 2*a[1]*(u[1] + 1)*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((1 - u[0])**2 + (u[1] + 1)**2) - g[1]])
    gn = min(yc-u[1],0)

    dgndu = np.array([0,1])
    if cub is None:
        fN = kn*gn*dgndu
    else:
        fN = kn/2*(2*gn-cub)*dgndu if gn<cub else kn/(2*cub)*gn**2*dgndu

    gt = abs(xp-u[0]) if gn<0 else 0
    gt_c = mu*(kn/kt)*abs(gn)
    elastic = kt*gt<=mu*abs(kn*gn)

    dgtdu = np.array([1,0])
    fT = kt*gt*dgtdu if elastic else mu*kn*gn*dgtdu

    return fint + fN + fT

"""

    def d2m1(u,g,a):
        return [[2*(a[0]*(u[0] + 1)**2/((u[0] + 1)**2 + (u[1] + 1)**2) - a[0]*(u[0] + 1)**2*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/((u[0] + 1)**2 + (u[1] + 1)**2)**(3/2) + a[0]*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((u[0] + 1)**2 + (u[1] + 1)**2) + a[1]*(1 - u[0])*(u[0] - 1)*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/((1 - u[0])**2 + (u[1] + 1)**2)**(3/2) + a[1]*(u[0] - 1)**2/((1 - u[0])**2 + (u[1] + 1)**2) + a[1]*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((1 - u[0])**2 + (u[1] + 1)**2)), 2*(u[1] + 1)*(a[0]*(u[0] + 1)/((u[0] + 1)**2 + (u[1] + 1)**2) - a[0]*(u[0] + 1)*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/((u[0] + 1)**2 + (u[1] + 1)**2)**(3/2) + a[1]*(1 - u[0])*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/((1 - u[0])**2 + (u[1] + 1)**2)**(3/2) + a[1]*(u[0] - 1)/((1 - u[0])**2 + (u[1] + 1)**2))], [2*(u[1] + 1)*(a[0]*(u[0] + 1)/((u[0] + 1)**2 + (u[1] + 1)**2) - a[0]*(u[0] + 1)*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/((u[0] + 1)**2 + (u[1] + 1)**2)**(3/2) + a[1]*(u[0] - 1)/((1 - u[0])**2 + (u[1] + 1)**2) - a[1]*(u[0] - 1)*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/((1 - u[0])**2 + (u[1] + 1)**2)**(3/2)), 2*(a[0]*(u[1] + 1)**2/((u[0] + 1)**2 + (u[1] + 1)**2) - a[0]*(u[1] + 1)**2*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/((u[0] + 1)**2 + (u[1] + 1)**2)**(3/2) + a[0]*(sqrt((u[0] + 1)**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((u[0] + 1)**2 + (u[1] + 1)**2) + a[1]*(u[1] + 1)**2/((1 - u[0])**2 + (u[1] + 1)**2) - a[1]*(u[1] + 1)**2*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/((1 - u[0])**2 + (u[1] + 1)**2)**(3/2) + a[1]*(sqrt((1 - u[0])**2 + (u[1] + 1)**2) - sqrt(2))/sqrt((1 - u[0])**2 + (u[1] + 1)**2))]]


    def mCN(u,kn,yc):
        gn = u[1]-yc
        return 0.5*kn*gn**2 if gn>0 else 0

    def fCN(u,kn,yc):
        gn = u[1]-yc
        dmdgn = -kn*gn if gn>0 else 0
        dgndu = np.array([0,1])
        return dmdgn*dgndu


    def mCT(u,kt,xp):
        gt = u[0]-xp
        return 0.5*kt*gt**2

    def fCT(u,kt,xp):
        gt = u[0]-xp
        dmdgt = -kt*gt
        dgtdu = np.array([1,0])
        return dmdgt*dgtdu


"""



a = np.ones(2)
g = np.zeros(2)
u0 = np.array([-0.2,-0.015])
kn = 100
kt = 100
mu = 0.2
cub = None
xp= -1.0
yc = -0.01

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# U1 = [-0.1,0.1]
U1 = [-0.1,0.1]
n1 = 20
eps = 5.0e-10
U2 = [-0.1,0.1]
n2 = 5000

x = np.zeros((n1,n2),dtype = float)
y = np.zeros((n1,n2),dtype = float)
z = np.zeros((n1,n2),dtype = float)
c = np.zeros((n1,n2,4),dtype = float)


minM = 10**100

for i in range(n1):
    for j in range(n2):
        u1 = ((n1-i)*U1[0] + i*(U1[1]))/n1
        u2 = ((n2-j)*U2[0] + j*(U2[1]))/n2

        x[i,j] = u1
        y[i,j] = u2
        current_m,elas = m1(np.array([u1,u2]),g,a,kn,kt,mu,xp,yc,cub=cub,showGN=False)
        minM = min(minM,current_m)
        z[i,j]=current_m
        c[i,j] = (0,0,1,0.5) if elas else (1,0,0,0.5)
        
ax.plot_surface(x,y,z,facecolors=c)
ax.set_xlabel("u1")
ax.set_ylabel("u2")
ax.set_zlabel("m")
# plt.show()





# a = np.ones(2)
# g = np.ones(2)

alpha_init = 0.8
c = 0.2
r = 0.7
nreset = 200

u = np.array(u0)
mnew = m1(u, g, a, kn, kt, mu, xp, yc )[0]
alpha = alpha_init
mold = 10**100
fnew = 0
hnew = 0
cnt = 0

while mnew<mold:
    mold=mnew
    fold = fnew
    hold = hnew

    fnew=dm1(u,g,a,kn,kt,mu,xp,yc,cub=cub,showGN=False)
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
    while mx>mnew+c*alpha*(hnew@fnew):
        alpha = r*alpha
        ux = u + alpha*hnew
        mx = m1(ux,g,a,kn,kt,mu,xp,yc,cub=cub, )[0]
        xniter +=1

    mnew = mx
    u = ux
    cnt += 1

m = mold
u = u - alpha*hnew

print("with u0=",u0,"\tcubic=",cub, ":")
print("")
print("u:", u)
print("m ( CG ):", m)
print("m (grid):", minM)
print("gn:", u[1]-yc)
print("niter", cnt)
print("xniter", xniter)

ax.scatter(u[0],u[1],m,s=5,c='r')
plt.show()