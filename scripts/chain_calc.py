import numpy as np
from scipy import constants as const
import sys
sys.path.append('../src/')
import green as gr
import time
import functions as fun
import sys

class chain():
    def __init__(self,N,pitch_x,direction=(1,0),U=0,alpha=0.04,spiral = 0,mode=0,m=18.7,pf = 0.21,delta_s = 1e-3,gamma_s=20e-6,E_px=100,E_range=5,V_range=2.8) -> None:
        self.N = N
        self.direction = direction
        self.pitch = pitch_x
        self.mode = mode
        self.m=m
        self.V_range = V_range
        self.pf=pf
        self.delta_s = delta_s
        self.gamma_s = gamma_s
        self.U = np.zeros(self.N)+U
        self.c = const.physical_constants['Hartree energy'][0]/const.e
        self.alpha = np.zeros(self.N)+alpha
        self.coords = self.coord_gen()
        self.E_px = E_px
        # create angles
        self.angles = []
        a = 0
        for i in range(0,N):
            a += spiral
            self.angles.append(a)
        self.E = np.linspace(-E_range*self.delta_s,E_range*self.delta_s,self.E_px)

        self.sim = gr.lattice(self.N , self.alpha , self.angles , self.coords,U=self.U ,m=self.m,pf=self.pf,delta=self.delta_s/self.c,mode=self.mode)
        self.par = { # to save parameters
            'N' : N,
            'direction' : direction,
            'pitch_x' : pitch_x,
            'mode' : mode,
            'm' : m,
            'pf' : pf,
            'delta_s' : delta_s,
            'gamma_s' : gamma_s,
            'U' : U,
            'alpha' : alpha,
            'E' : [self.E[-1],E_px],
            'angles' : spiral,
        }


    def coord_gen(self):
        n = 0
        coords = []
        for i in range(0,self.N):
            coords.append([n*self.direction[0],n*self.direction[1]])
            n+=self.pitch
        return coords


    def didv(self,coord):
        spec = []
        if self.U[0] != 0:
            for k in self.E/self.c:
                spec.append(np.sign(k)*self.sim.ElecDOS(coord, k + self.gamma_s*1j*np.sign(k)/self.c))
        else:
            for k in self.E/self.c:
                spec.append(np.sign(k)*(self.sim.DOS(coord, k + self.gamma_s*1j*np.sign(k)/self.c)))
        
        return np.array(spec/spec[0])

    def didv_conv(self,coord,Delta_t,Gamma_t):
        dos = self.didv(coord)
        self.V = np.linspace(-self.delta_s*5,self.delta_s*5,self.E_px)
        conv_dos = fun.dynesConvolute(self.V,self.E,dos,Delta_t,1.3,Gamma_t)
        return np.array(conv_dos)

    def linescan(self,density):
        x = self.pitch
        y = self.pitch*self.direction[1]
        self.LSx = np.linspace(self.coords[0][0]-x,self.coords[-1][0]+x,(self.N)*density)
        self.LSy = np.linspace(self.coords[0][1]-y,self.coords[-1][1]+y,(self.N)*density)
        
        #timecalc
        t0 = time.time()
        self.didv((self.LSx[0],self.LSy[0]))
        t1 = time.time()
        total_time = (t1-t0)*len(self.LSx)
        print('Simulation time = {}'.format(np.round(total_time/60,2)))
        ####

        self.LS = []

        for i in range(0,len(self.LSx)):
            self.LS.append(self.didv((self.LSx[i],self.LSy[i])))
        self.LS = np.array(self.LS)

        
    def LSconvolute(self,Delta_t,Gamma_t):
        self.LSC = np.zeros(self.LS.shape)
        self.V = np.linspace(-self.delta_s*self.V_range,self.delta_s*self.V_range,self.E_px)
        for i in range(self.LS.shape[0]):
            self.LSC[i,:] = (fun.dynesConvolute(self.V,self.E,self.LS[i,:],Delta_t,1.3,Gamma_t))
    
