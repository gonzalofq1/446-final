from scipy import sparse
from timesteppers import StateVector, CrankNicolson, RK22
import finite
import numpy as np
import timesteppers
import scipy.sparse.linalg as spla
import timesteppers


class SchrodingerBCNonLinear:

    def __init__(self, c, spatial_order, domain,g):
        self.c = c
        self.X = timesteppers.StateVector([c])
        self.t = 0
        self.iter = 0
       
        x = domain.grids[0]
        y= domain.grids[1]
                
        dx = finite.DifferenceUniformGrid(1, spatial_order, x)
        dy = finite.DifferenceUniformGrid(1, spatial_order, y)
        dx2 = finite.DifferenceUniformGrid(2, spatial_order, x)
        dy2 = finite.DifferenceUniformGrid(2, spatial_order, y)


        diffx = self.Diffusion(c, dx,dx2,0)
        diffy = self.Diffusion(c,dy,dy2,1)
        rec = self.Reaction(c,g,x,y)
        
        self.ts_x = timesteppers.CrankNicolson(diffx, 0)
        self.ts_y = timesteppers.CrankNicolson(diffy, 1)
        self.ts_rec = timesteppers.RK22(rec)

    class Diffusion:

        def __init__(self, c, dx,dx2,axis):
            self.X = timesteppers.StateVector([c], axis)
            N = c.shape[axis]
            M = sparse.eye(N, N,dtype=complex)
            M=M.tocsr()
            M[0,0]=0
            M[-1,-1]=0
            M.eliminate_zeros()
            self.M = M

            L = dx2.matrix.astype(complex)
            L=L.tocsr()
            L = L*(-1j)/2
            L[0,:] = 0
            L[-1,:] = 0
            L[0,0] = 1
            L[-1,-1] = 1
            L.eliminate_zeros()
            self.L = L

    
    class Reaction:

        def __init__(self,c,g,x,y):
            self.X = timesteppers.StateVector([c])

            N = len(c)
            I = sparse.eye(N, N,dtype=complex)
            Z = sparse.csr_matrix((N, N),dtype=complex)

            M00 = I
            M01 = Z
            M10 = Z
            M11 = I
            self.M = sparse.bmat([[M00, M01],
                                  [M10, M11]])

            L00 = Z
            L01 = Z
            L10 = Z
            L11 = Z
            self.L = sparse.bmat([[L00, L01],
                                  [L10, L11]])


            def f(X):
                rec = X.data*abs(X.data)**2
                return -1j*rec*g

            self.F =f

            def BC(X):
                X.data[0,:]=0
                X.data[-1,:]=0
                X.data[:,-1]=0
                X.data[:,0]=0
            self.BC =BC

    
    def step(self, dt):
        self.t = self.t+dt
        self.iter = self.iter+1
    
     
        self.ts_y.step(dt/2)
        self.ts_x.step(dt/2)
        self.ts_rec.step(dt/2)
        self.ts_rec.step(dt/2)
        self.ts_x.step(dt/2)
        self.ts_y.step(dt/2)
        
        
class SchrodingerBCLinearSlit:

    def __init__(self, c, spatial_order, domain,g):
        self.c = c
        self.X = timesteppers.StateVector([c])
        self.t = 0
        self.iter = 0
       
        x = domain.grids[0]
        y= domain.grids[1]
                
        dx = finite.DifferenceUniformGrid(1, spatial_order, x)
        dy = finite.DifferenceUniformGrid(1, spatial_order, y)
        dx2 = finite.DifferenceUniformGrid(2, spatial_order, x)
        dy2 = finite.DifferenceUniformGrid(2, spatial_order, y)


        diffx = self.Diffusion(c, dx,dx2,0)
        diffy = self.Diffusion(c,dy,dy2,1)
        rec = self.Reaction(c,g)
        
        self.ts_x = timesteppers.CrankNicolson(diffx, 0)
        self.ts_y = timesteppers.CrankNicolson(diffy, 1)
        self.ts_rec = timesteppers.RK22(rec)

    class Diffusion:

        def __init__(self, c, dx,dx2,axis):
            self.X = timesteppers.StateVector([c], axis)
            N = c.shape[axis]
            M = sparse.eye(N, N,dtype=complex)
            M=M.tocsr()
            M[0,0]=0
            M[-1,-1]=0
            M.eliminate_zeros()
            self.M = M

            L = dx2.matrix.astype(complex)
            L=L.tocsr()
            L = L*(-1j)/2
            L[0,:] = 0
            L[-1,:] = 0
            L[0,0] = 1
            L[-1,-1] = 1
            L.eliminate_zeros()
            self.L = L
            
    class Reaction:

        def __init__(self,c,g):
            self.X = timesteppers.StateVector([c])

            N = len(c)
            I = sparse.eye(N, N,dtype=complex)
            Z = sparse.csr_matrix((N, N),dtype=complex)

            M00 = I
            M01 = Z
            M10 = Z
            M11 = I
            self.M = sparse.bmat([[M00, M01],
                                    [M10, M11]])

            L00 = Z
            L01 = Z
            L10 = Z
            L11 = Z
            self.L = sparse.bmat([[L00, L01],
                                    [L10, L11]])


            def f(X):
                rec = X.data*abs(X.data)**2
                return -1j*rec*g

            self.F =f

            def BC(X):
                m = int(N/2)
                X.data[0:m-7,m-10:m+10]=0
                X.data[m+7:-1,m-10:m+10]=0
            self.BC =BC

    
    def step(self, dt):
        self.t = self.t+dt
        self.iter = self.iter+1
    
       
        self.ts_y.step(dt/2)
        self.ts_x.step(dt/2)
        self.ts_rec.step(dt/2)
        self.ts_rec.step(dt/2)
        self.ts_x.step(dt/2)
        self.ts_y.step(dt/2)
           
