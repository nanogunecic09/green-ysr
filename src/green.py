import numpy as np
import scipy.special as scisp
from mpmath import struveh
from scipy import constants as const

class lattice():

    def __init__(self,N,alpha,theta,r__,U=0,m=20.956,pf=0.274,delta=0.0000287,gamma=40e-6,mode=1):
        self.N = N
        self.s0 = np.array([[1,0],[0,1]])
        self.s1 = np.array([[0,1],[1,0]])
        self.s3 = np.array([[1,0],[0,-1]])
        self.theta_ = theta # theta of the spins
        self.alpha_ = alpha # alpha of the atoms in order
        self.U_ = U # potential scattering
        self.r__ = r__ #atomic positions
        self.m = m # mass electron
        self.pf = pf # fermi momentum
        self.delta = delta # delta superconductor
        self.gamma = gamma
        self.mode = mode
    def G0(self,r1,r2,E):
        delta = self.delta
        w = np.sqrt(delta**2-E**2)
        m=self.m
        pf=self.pf
        BCS = E/w*np.kron(self.s0,self.s0) + delta/w*np.kron(self.s1,self.s0)
        xi=np.kron(self.s3,self.s0)
        if self.mode==0:
            x = r1[0]-r2[0]
            y = r1[1]-r2[1]
            x1 = np.sqrt(x**2+y**2)
            ## mode=0 is for the spherical Fermi surface ##
            ## for this mode x1 is the radial distance and x2 is the angle ##
            u=complex(x1*pf , x1*m*w/(pf))
            a=scisp.jv(0,u) #bessel function 0th order u=argument
            b=complex(struveh(0,u)) #gives the struve function
            G0=m/2*(complex(np.real(a+complex(0,1)*b),0))*BCS + m/2*(complex(np.imag(a+complex(0,1)*b),0)) * xi
     
        if self.mode==1:
            x1=r1[0]-r2[0]
            x2=r1[1]-r2[1]
            delta = self.delta
            m=self.m
            pf=self.pf
            if np.abs(x1)<0.00001 and np.abs(x2)<0.00001:
                    G1=(2*m/np.pi)*BCS
                    G2=0.0
            elif np.abs(x2)<0.00001:
                G1=(np.exp(-(m*w/pf)*np.abs(x1))*( (1/np.abs(x1))*pf*np.sin(pf*(np.abs(x1)))+(pf**2)*np.cos(pf*np.abs(x1)) ) )*(m/(np.pi*pf**2))*BCS
                G2=(np.exp(-(m*w/pf)*np.abs(x1))*( (1/np.abs(x1))*pf*np.cos(pf*(np.abs(x1)))+(pf**2)*np.sin(pf*np.abs(x1)) ) -pf/np.abs(x1))*(m/(np.pi*pf**2))*xi
            elif np.abs(x1)<0.00001:
                G1=(np.exp(-(m*w/pf)*np.abs(x2))*( (1/np.abs(x2))*pf*np.sin(pf*(np.abs(x2)))+(pf**2)*np.cos(pf*np.abs(x2)) ) )*(m/(np.pi*pf**2))*BCS
                G2=(np.exp(-(m*w/pf)*np.abs(x2))*( (1/np.abs(x2))*pf*np.cos(pf*(np.abs(x2)))+(pf**2)*np.sin(pf*np.abs(x2)) ) -pf/np.abs(x2))*(m/(np.pi*pf**2))*xi
            elif x1+x2>=0.00001 and x1-x2>=0.00001:
                G1=(np.exp(-(m*w/pf)*(x1+x2))*(1/x1+1/x2)*pf*np.sin(pf*(x1+x2))+np.exp(-(m*w/pf)*(x1-x2))*(1/x1-1/x2)*pf*np.sin(pf*(x1-x2)) )*(m/(2*np.pi*pf**2))*BCS
                G2=(np.exp(-(m*w/pf)*(x1+x2))*(1/x1+1/x2)*pf*np.cos(pf*(x1+x2))+np.exp(-(m*w/pf)*(x1-x2))*(1/x1-1/x2)*pf*np.cos(pf*(x1-x2)) -2*pf/x1)*(m/(2*np.pi*pf**2))*xi
            elif x1+x2>0.00001 and x1-x2<0.00001:
                G1=(np.exp(-(m*w/pf)*(x1+x2))*(1/x1+1/x2)*pf*np.sin(pf*(x1+x2))+np.exp(-(m*w/pf)*(x2-x1))*(1/x2-1/x1)*pf*np.sin(pf*(x2-x1)) )*(m/(2*np.pi*pf**2))*BCS
                G2=(np.exp(-(m*w/pf)*(x1+x2))*(1/x1+1/x2)*pf*np.cos(pf*(x1+x2))+np.exp(-(m*w/pf)*(x2-x1))*(1/x2-1/x1)*pf*np.cos(pf*(x2-x1)) -2*pf/x2)*(m/(2*np.pi*pf**2))*xi
            elif x1+x2<=0.00001 and x1-x2<=0.00001:
                G1=(np.exp((m*w/pf)*(x1+x2))*(-1/x1-1/x2)*pf*np.sin(pf*(-x1-x2))+np.exp(-(m*w/pf)*(-x1+x2))*(-1/x1+1/x2)*pf*np.sin(pf*(-x1+x2)) )*(m/(2*np.pi*pf**2))*BCS
                G2=(np.exp((m*w/pf)*(x1+x2))*(-1/x1-1/x2)*pf*np.cos(pf*(-x1-x2))+np.exp(-(m*w/pf)*(-x1+x2))*(-1/x1+1/x2)*pf*np.cos(pf*(-x1+x2)) +2*pf/x1)*(m/(2*np.pi*pf**2))*xi
            elif x1+x2<0.00001 and x1-x2>0.00001:
                G1=(np.exp((m*w/pf)*(x1+x2))*(-1/x1-1/x2)*pf*np.sin(pf*(-x1-x2))+np.exp(-(m*w/pf)*(-x2+x1))*(-1/x2+1/x1)*pf*np.sin(pf*(-x2+x1)) )*(m/(2*np.pi*pf**2))*BCS
                G2=(np.exp((m*w/pf)*(x1+x2))*(-1/x1-1/x2)*pf*np.cos(pf*(-x1-x2))+np.exp(-(m*w/pf)*(-x2+x1))*(-1/x2+1/x1)*pf*np.cos(pf*(-x2+x1)) +2*pf/x2)*(m/(2*np.pi*pf**2))*xi
            G0 = G1+G2
        return G0
    
    def V(self,theta,alpha,U):
        return alpha * np.cos(theta)*np.kron(self.s0,self.s3) + alpha*np.sin(theta)*np.kron(self.s0,self.s1) + U*np.kron(self.s3,self.s0)
    
    def M(self,E):
        n=self.N
        M = np.zeros((4*n,4*n),dtype=complex)
        for i in range(0,4*n,4):
            for j in range(0,4*n,4):
                M[i:i+4,j:j+4] = np.dot(self.G0(self.r__[i//4],self.r__[j//4],E),self.V(self.theta_[j//4],self.alpha_[j//4],self.U_[j//4]))
        return M

    def G0_(self,r_,E):
        G0_ = np.zeros((4*self.N,4),dtype=complex)
        n=0
        for i in range(0,4*self.N,4):
            G0_[i:i+4,0:4] = self.G0(self.r__[n],r_,E)
            n+=1
        return G0_
    
    def G(self,r_,E):
        MM = np.linalg.inv(np.identity(self.N*4,dtype=complex)-self.M(E))
        GG = np.dot(MM,self.G0_(r_,E))
        G = np.zeros((4,4),dtype=complex)
        G += self.G0((0,0),(0,0),E)
        n=0
        for i in range(0,4*self.N,4):
            G += np.dot(np.dot( self.G0(r_,self.r__[n],E) , self.V(self.theta_[n] ,self.alpha_[n],self.U_[n])), GG[i:i+4,0:4])
            n+=1
        return G

    def DOS(self,r_,E):
        return np.imag(np.trace(self.G(r_,E)))

    def ElecDOS(self,r_,E):
        return np.imag(np.trace(np.dot(self.G(r_,E),np.diag((1,1,0,0)))))

    def HoleDOS(self,r_,E):
        return np.imag(np.trace(np.dot(self.G(r_,E),np.diag((0,0,1,1)))))

    def rhoss(self,L,E):
        self.alphaL_ = self.alpha_*L
        rhoss = 0
        for n in range(0,self.N):
            G = -np.imag(self.G(self.r__[n],E))
            rhoss += self.alphaL_[n]*( (G[0,0]-G[1,1]) * np.cos(self.theta_[n]) + (G[1,0]-G[0,1]) * np.sin(self.theta_[n]) )
        return rhoss
    def spectra(self,r_):
        c = const.physical_constants['Hartree energy'][0]/const.e
        E = np.linspace(-4*self.delta ,4*self.delta,300)/c
        spectra=[]
        for i in E:
            spectra.append(np.sign(i)*self.DOS( r_, i + self.gamma*1j*np.sign(i)/c))
        return spectra

    def rho4(self,r_,E):
        return np.imag(self.G(r_,E))