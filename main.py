'''
This script finds a numerical approximation for the solution to
 the following 1D heat equation with Neumann conditions:
   u_t = u_xx  for  x \in (0,1),  t \in (0,T),
   u(x,t=0) = x,
   u'(x=0,t) = 0, u'(x=1,t) = 0,
 using a backward Euler scheme.
''' 
def U_exact2(x,t):
    M = np.size(x)
    u_ex = 0.5*np.ones(M)  
            
    for s in range(1,1000,2):
        npi= s*np.pi
        c_n = -(4.0)/npi/npi                            
        u_ex = u_ex + c_n*np.cos(npi*x)*np.exp(-1*npi*npi*t)
        
    return u_ex
