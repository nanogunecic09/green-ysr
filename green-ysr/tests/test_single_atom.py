import green_ysr as gr
import numpy as np
import matplotlib.pyplot as plt


sim = gr.chain(1,0,alpha=0.05,gamma_s=50e-6,U=0,E_px=500,E_range=2)

plt.plot(sim.E*1e3,sim.didv([0,0]))

plt.xlabel('Energy (meV)')
plt.ylabel('LDOS (Gn)')
plt.tick_params(axis='both',direction='in')
gr.set_size_cm(8,5)
plt.tight_layout()