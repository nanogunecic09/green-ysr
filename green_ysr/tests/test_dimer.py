import matplotlib.pyplot as plt
import green_ysr as gr


sim = gr.chain(2,8,gamma_s=50e-6,U=0,E_px=500,E_range=2)
plt.plot(sim.E*1e3,sim.didv([0,0]))

plt.xlabel('Energy (meV)')
plt.ylabel('LDOS (Gn)')
