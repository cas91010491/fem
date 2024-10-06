import torch
from pdb import set_trace

def jacobian(y, x, create_graph=False):                                                            
    jac = []                                                                                       
    flat_y = y.reshape(-1)                                                                         
    grad_y = torch.zeros_like(flat_y)                                                                 
    for i in range(len(flat_y)):                                                                   
        grad_y[i] = 1.                                                                             
        grad_x, = torch.autograd.grad(flat_y, x, grad_y, retain_graph=True, create_graph=create_graph)
        jac.append(grad_x.reshape(x.shape))                                                        
        grad_y[i] = 0.                                                                                
    return torch.stack(jac).reshape(y.shape + x.shape)        

def hessian(y, x):                                                                                 
    return jacobian(jacobian(y, x, create_graph=True), x)                                          



xac = torch.tensor([0.0, 1.0], dtype = torch.float64, requires_grad=True)
x1b = torch.tensor([0.0, 0.0], dtype = torch.float64, requires_grad=True)
x2b = torch.tensor([2.0, 2.0], dtype = torch.float64, requires_grad=True)
xs  = torch.tensor([2.0, 0.0], dtype = torch.float64, requires_grad=True)

L = torch.norm(x2b-x1b)         # length of 'patch'

tau = (x2b-x1b)/L               # unit vector along patch, from 1 to 2
Lf1 = torch.dot(xs-x1b,tau)      # 'left' distance of patch (until xc)
t = L1/L                        # parametric coordinate

xc = x1b+(t*L)*tau              # projection point 

gt = torch.norm(xc-xac)

set_trace()

