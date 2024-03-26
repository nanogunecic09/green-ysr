import sys
import matplotlib.pyplot as plt
import green_ysr as gr
import numpy as np

map = []
Js = np.linspace(0,0.4,20)

for J in Js:
    sim = gr.chain(1,8,J=J)
    map.append(sim.didv([0,0]))



plt.imshow(map,aspect='auto',extent=[sim.E.min()*1e3,sim.E.max()*1e3,0,0.1])
plt.xlabel('Energy (meV)')
plt.ylabel('LDOS (Gn)')