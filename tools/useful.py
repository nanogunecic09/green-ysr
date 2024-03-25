import pickle
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL) 
def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)

def sim_save(sim):
    keys = ['N','direction','pitch_x','mode','U','J','angles']
    fname = '../out/sim'
    for i in keys:
        fname  =  fname +'_'+ i  + '{}'.format(sim.par[i])
    save_obj(sim,fname)

def plot_LS(chain_sim):
    plt.figure()
    plt.imshow(chain_sim.LS,extent=[chain_sim.E[0]/chain_sim.par['delta_s'],chain_sim.E[-1]/chain_sim.par['delta_s'],0,1],aspect='auto')

def set_size_cm(w,h, ax=None):
    """ w, h: width, height in cm """
    cm = 1/2.54
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w*cm)/(r-l)
    figh = float(h*cm)/(t-b)
    ax.figure.set_size_inches(figw, figh)

import matplotlib
def spines(ax):
    plt.setp(ax.spines.values(), lw=0.3)
    ax.tick_params(width=0.3)
    new_rc_params = {'text.usetex': False,
    "svg.fonttype": 'none'
    }
    matplotlib.rcParams.update(new_rc_params)
    matplotlib.rcParams['axes.unicode_minus']=False


class grid():

    def __init__(self,bias,didv_map):
        self.conductance = didv_map
        self.bias = bias

    
    def explorer(self):
        self.figure = plt.figure(figsize=(6,6))
        self.axMap = self.figure.add_subplot(1,1,1)
        self.figure.subplots_adjust(bottom=0.35)
        self.ax1 = self.figure.add_axes([0.20, 0.10, 0.65, 0.03])
        self.ax2 = self.figure.add_axes([0.20, 0.15, 0.65, 0.03])
        self.ax3 = self.figure.add_axes([0.20, 0.20, 0.65, 0.03])
        self.energyCut_slider = Slider(self.ax1,'Energy cut',self.bias.min()*1e3,self.bias.max()*1e3,valinit=0, valstep=(0.01))
        self.smin_slider = Slider(self.ax2, 'Min', self.conductance.min(), self.conductance.max(), valinit =self.conductance.min())
        self.smax_slider = Slider(self.ax3, 'Max', self.conductance.min(), self.conductance.max(), valinit =self.conductance.max()*0.5)            
        self.energyCut_slider.on_changed(self.update_energy)
        self.smin_slider.on_changed(self.update_cscale)
        self.smax_slider.on_changed(self.update_cscale)
        self.conductance_cut = self.conductance[:,:,0]
        self.im1 = self.axMap.imshow(np.flipud(self.conductance_cut),interpolation='nearest')
        #axis labels
        self.axMap.set_xlabel('x (nm)')
        self.axMap.set_ylabel('y (nm)')

    def update_energy(self,val):
        self.cutIdx = (abs(self.bias-val*1e-3)).argmin()
        self.im1.set_data(np.flipud(self.conductance[:,:,self.cutIdx]))
        self.im1.set_clim(self.conductance[:,:,self.cutIdx].min(),self.conductance[:,:,self.cutIdx].max())
        # self.label.set_text('{} mV'.format(np.round(val,2)))
        self.figure.canvas.draw()

    def update_cscale(self,val):
        self.im1.set_clim([self.smin_slider.val,self.smax_slider.val])
        self.figure.canvas.draw()