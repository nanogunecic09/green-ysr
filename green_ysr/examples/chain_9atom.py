import matplotlib.pyplot as plt
import green_ysr as gr
import scipy.ndimage
import numpy as np

chain_sim = gr.chain(9,6.3,J=0.03,direction=(1,2),mode=1,U=0,spiral=0,m=18.7,pf = 0.28) # load the simulator
chain_sim.linescan(3)# compute the linescan with 3 points per atom
chain_sim.LSconvolute(0.65e-3,40e-6) # convolute linescan with a supercondcting tip [delta,dynes]
gr.sim_save(chain_sim) # save the calculation in the output

# loading saved simulation
chain_sim = gr.load_obj('green_ysr/out/sim_type1D_N9_direction(1, 2)_pitch_x6.3_mode1_U0_alpha0.03_angles0')


# account for non-local tunneling effects
sigma = [1, 0]
y = scipy.ndimage.filters.gaussian_filter(np.array(chain_sim.LSC), sigma, mode='constant')

# plot the calculation
f, ax = plt.subplots()
ax.imshow(y,extent=[chain_sim.V[0]/chain_sim.par['delta_s'],chain_sim.V[-1]/chain_sim.par['delta_s'],0,7],aspect='auto')
ax.set_xlabel('Bias (mV)')
ax.set_ylabel('Distance (nm)')
ax.tick_params(axis='both',direction='in')
gr.set_size_cm(8,5)
ax.set_xlim()
