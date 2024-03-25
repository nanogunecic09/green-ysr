import numpy as np
from scipy import constants as const
import sys
sys.path.append('../src/')
import green as gr
import time
import functions as fun
import sys
import matplotlib.pyplot as plt

class lattice():
    def __init__(self,N,coords,direction=(1,0),U=0,J=0.04,spiral = 0,mode=0,m=18.7,pf = 0.21,delta_s = 0.78e-3,gamma_s=20e-6) -> None:
        self.N = N
        self.direction = direction
        self.mode = mode
        self.m=m
        self.pf=pf
        self.delta_s = delta_s
        self.gamma_s = gamma_s
        self.U = np.zeros(self.N)+U
        self.c = const.physical_constants['Hartree energy'][0]/const.e
        self.J = np.zeros(self.N)+J
        self.coords = coords
        # create angles
        self.angles = []
        a = 0
        for i in range(0,N):
            a += spiral
            self.angles.append(a)
        self.E = np.linspace(-5*self.delta_s,5*self.delta_s,100)
        self.sim = gr.lattice(self.N , self.J , self.angles , self.coords,U=self.U ,m=self.m,pf=self.pf,delta=self.delta_s/self.c,mode=self.mode)
        self.par = { # to save parameters
            'N' : N,
            'direction' : direction,
            'coords' : coords,
            'mode' : mode,
            'm' : m,
            'pf' : pf,
            'delta_s' : delta_s,
            'gamma_s' : gamma_s,
            'U' : U,
            'J' : J,
            'E' : [1.5*self.delta_s,100],
            'angles' : spiral,
        }


    def show_lattice(self):
        f,ax = plt.subplots(1)
        for i in self.coords:
            ax.scatter(i[0],i[1],color='C0')

    def map_coord_gen(self,spac,resolution,size): # spac is the point spacing of the grid, resolution is the number of points in one line, size is the length in units of spacing 
        self.resolution = resolution
        A = np.arange(-spac*size/2,spac*size/2+(spac*size)/resolution,(spac*size)/resolution)
        Gx = np.meshgrid(A,A)[0]
        Gy = np.meshgrid(A,A)[1]
        self.map_coords = [Gx,Gy]

    def show_lattice_map(self):
        f,ax = plt.subplots(1)
        for i in range(len(self.map_coords[0])):
            ax.scatter(self.map_coords[0][i],self.map_coords[1][i],color='C1')    
        for i in self.coords:
            ax.scatter(i[0],i[1],color='C0')



    def didv(self,coord):
        spec = []
        if self.U[0] != 0:
            for k in self.E/self.c:
                spec.append(np.sign(k)*self.sim.ElecDOS(coord, k + self.gamma_s*1j*np.sign(k)/self.c))
        else:
            for k in self.E/self.c:
                spec.append(np.sign(k)*(self.sim.DOS(coord, k + self.gamma_s*1j*np.sign(k)/self.c)))
        
        return np.array(spec)
    
    def didv_conv(self,coord,Delta_t,Gamma_t):
        dos = self.didv(coord)
        self.V = np.linspace(-self.delta_s*5,self.delta_s*5,100)
        conv_dos = fun.dynesConvolute(self.V,self.E,dos,Delta_t,1.3,Gamma_t)
        return np.array(conv_dos)

    def didv_map_calc(self):
        #timecalc
        t0 = time.time()
        self.didv([0,0])
        t1 = time.time()
        total_time = (t1-t0)*self.resolution*self.resolution
        print('Simulation time = {}'.format(np.round(total_time/60,2)))
        ####
        self.didv_map = np.zeros((self.resolution,self.resolution,self.E.shape[0]))
        for i in range(self.resolution):
            for j in range(self.resolution):
                self.didv_map[i,j,:] = self.didv([self.map_coords[0][i,j],self.map_coords[1][i,j]])

    