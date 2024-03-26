# Green_YSR documentation

This is a basic tutorial with few examples on how to calculate LDOS using this green function model.

The aim of this package is the calculation of the local density of states (LDOS) of magnetic impurities interacting with superconductors. The model is a standard Green function model able to deal with anisotropic Fermi surfaces. See [this paper](https://doi.org/10.1103/PhysRevB.105.245403) for details on the model.

## Requirements

Clone the following repo: [Repo Green_YSR](https://github.com/nanogunecic09/Green_YSR)

I suggest using Anaconda python, the calculations are based on standard numpy, scipy, math pkgs.

The examples in the test folder are written with .ipynb files. VS code, among others, provides an add-on to visualize and work with notebooks.

## Quick tutorial

The simplest calculation is the spectra of a single atom. You can find an example in test/single_atom.ipynb.

To load the simulator use:
    
    sim = ch.chain(1,0,alpha=0.05,gamma_s=50e-6,U=0,E_px=500,E_range=2)

The initialization parameters are in order:
1. Numbe of atoms, in this case 1.
2. Spacing between the atoms, in this case is not used.
3. alpha: is the adimensional exchange coupling.
4. gamma_s : dynes parameter
5. U : asimentional potential scattering parameter
6. E_px : energy points
7. E_range : energy range (E_range*Delta)

To calculate the LDOS in the point (0,0) and plot it:

    plt.plot(sim.E*1e3,sim.didv([0,0]))



The rest of the arguments of ch.chain are implicit and can be modified if needed. The material parameters that determine the period of the YSR wavefunction and the coherence length are: 'm' the effective mass and 'pf0' the Fermi momentum.



## About