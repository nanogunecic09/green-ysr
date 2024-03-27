## Reference

### *class* green_ysr.lattice(**kwargs)

>The atomic lattice to simulate.

**kwargs**:

type : string  
>type of lattice, 'single atom', '1D' or '2D'

N : int  
>Number of atoms in the system  

coords : list  
>in the case of a '2D' lattice, a list of the atomic coordinates must be provided

pitch_x : float  
>for a equi-spaced '1D' chain the pitch_x is the spacing of the atoms in the chain, needed to generate the coordinates

direction : tuple  
>direction of elongation of the chain in cristallographic notation, e.g. (1,0) is horizontal, (1,1) is 45 degrees

alpha : float  
>Adimensional exchange coupling parameter

theta : float  
>angle describing the orientation in rad of the classical spins

U : float  
>adimensional potential scattering parameter  

spiral : float  
>if inserted is the pace of the spin spiral along the chain in radiants  

m : float  
>effective mass multiplied by pi in atomic unit  

pf : float  
>fermi momentum in atomic units  

delta : float  
>superconducting order parameter in atomic units  

gamma : float  
>dynes parameter in atomic units  

mode : int  
>fermi surface dimensionality, 1 for circular, 2 for squared  

E_px : int  
>point of energy to calculate  

E_range : tuple  
>extremes of the energy, set 0,2 for calculate half of the spectrum  

V_range : tuple  
>extremes of the voltage, used when doing the convolution  

spin_texture : list  
>list of the spin orientation of the atoms in the structure, if not provided will be ferromagnetic  

**Attributes:**

E : array  
>energy array  

V : array  
>voltage array, for convolution. The V range has to be always lower than E for the convolution.  

sim : class  
>is the solver of the green function  

par : dict  
>save all parameters  

coord_gen()  
>for generating chain coordinates fiven the pitch  
>**return** array of coordinates

show_lattice()  
>show the lattice coordinates in a plot

map_coord_gen()  
>generate the coordinate for a didv map

**Parameters :**  

spac : float  
>spacing in atomic units of the lattice  

resolution : int  
>number of pixel in one line of the didv map  

size : int  
>size of one side of the didv map in units of spac (spacing)  

show_lattice_map()
>show the lattice and the didv map coordinates in a plot

didv()  
>calculate the didv in one point  
>**Parameters :**  
>coord : tuple  
>coordinates of the measurement point  
>**return** array with spectra (normalized to normal conductance)  

didv_map_calc()  
>caculates the didv map and store it in self.didv_map

didv_conv()  
>calculate the convolution with the superconducting tip at the temperature specified in self.T    
>**Parameters :**  
>coord : tuple  
>coordinates of the measurement point  

Delta_t : float  
>pairing amplitude in the tip electrode in meV  

Gamma_t : float  
>dynes parameter in meV  
>**return** array with tip convoluted spectra (normalized to normal conductance)  

linescan()  
>compute a line of spectra (saved in self.LS) along the 1D chain direction, including initial and final portion of superconductor (pitch_x)  
>**Parameters :**  
>density : int  
>number of points per atom in the linescan

LSconvolute()  
>compute the convoluted linescan (saved in self.LSC) at the temperature specified in self.T  
>**Parameters :**  
>coord : tuple  
>coordinates of the measurement point

Delta_t : float  
>pairing amplitude of the tip electrode in meV

Gamma_t : float  
>Dynes parameter of the tip

explorer()  
>opens a slider plot to explore self.didv_map in constant energy cuts

update_energy()  
>vital to explorer()

update_cscale()  
>vital to explorer()




### *class* green_ysr.green(**kwargs)

>To calculate the green function.

**Parameters :**

N : int  
>Number of atoms in the system  

alpha : float  
>Adimensional exchange coupling parameter  

theta : float  
>angle describing the orientation in rad of the classical spins  

r__ : list   
>set of atomic coordinates  

U : float  
>adimensional potential scattering parameter  

m : float  
>effective mass multiplied by pi in atomic unit  

pf : float  
>fermi momentum in atomic units  

delta : float  
>superconducting order parameter in atomic units  

gamma : float  
>dynes parameter in atomic units  

mode : int  
>fermi surface dimensionality, 1 for circular, 2 for squared  

**Attributes:**

s0-1-3 : numpy array  
>pauli matrices

G0()  
>analitical form of the unperturbed green function

V()  
>Potential matrix of each impurity
>**return** a 4x4 matrix

G0_()  
>green function of the Dyson series
>**return** a 4x4N matrix

G()  
>full green function
>**return** a 4Nx4N matrix

DOS()  
>compute the spectrum of the green function

ElecDOS()  
>Compute the electron part of the spectral function

ElecDOS()  
>Compute the hole part of the spectral function


### **Additional functions**

fdd()  
>Fermi dirac distribution

dynesdos()  
>Dynes dos of a superconductor

dynesConvolute()  
>numerical convolution of a conducance array

save_obj()  
>save an arbitrary object

load_obj()  
>load an arbitrary object given its path without extension

sim_save()  
>to save a simulation, including important parameters in the name

plot_LS()  
>to show the result of a line of spectra calculated with lattice()

set_size_cm()  
>to change the size of a figure axis

