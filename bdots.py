import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import scipy.integrate

class ScopeChannel:
    def __init__(self, shot, scope, channel):
        fn="/media/nest_files/magpie_scopes/scope"+scope+"_"+shot
        self.time=np.loadtxt(fn+"time")
        self.data=np.loadtxt(fn+"_"+channel)[1:]
        
class Bdot_pair:
    def __init__(self, shot, scope="1", bdot1='A1', bdot2='A2'):
        self.shot=shot
        #bdot signals 1 and 2
        self.bd1=ScopeChannel(shot, scope, bdot1)
        self.bd2=ScopeChannel(shot, scope, bdot2)
    def truncate(self, threshold=1.0, window=1000, cal=[1,1]):
        #find the start of the current pulse with a  high threshold
        sig1=self.bd1.data
        start=np.nonzero(abs(sig1)>threshold)[0][0]
        #back off a bit so we can see the zero signal
        self.start=start-100
        #reverse the array to find the end of the current pulse with a high threshold
        #end=np.nonzero(abs(sig1[::-1])>threshold)[0][0]
        #back off a bit so we can see the zero signal
        #end=end-100
        #self.end=sig1.size-end #find the index in the non-reversed array
        self.time=self.bd1.time[self.start:self.start+window]
        self.bd1_tr=self.bd1.data[self.start:self.start+window]*cal[0]
        self.bd2_tr=self.bd2.data[self.start:self.start+window]*cal[1]
        self.add()
        self.subtract()
    def add(self):
        self.estat=(self.bd1_tr+self.bd2_tr)/2.0      
    def subtract(self):
        self.dBdt=(self.bd1_tr-self.bd2_tr)/2.0
    def integrate(self):
        self.B=scipy.integrate.cumtrapz(self.dBdt,self.time)/1e9
        self.time_B=self.time[:-1]
    def plot(self, data, ax=None, flip=1, bdname=None):
        if ax is None:
            fig, ax=plt.subplots()
        if bdname is not None:
            b1=bdname[0:2]
            b2=bdname[0]+bdname[2]
        if data is "raw":
            t=self.bd1.time
            d1=self.bd1.data
            d2=self.bd2.data
            l1=b1+' raw'
            l2=b2+' raw'
        if data is "tr":
            t=self.time
            d1=self.bd1_tr
            d2=self.bd2_tr
            l1=b1+' truncated'
            l2=b2+' truncated'
        if data is "sum_diff":
            t=self.time
            d1=self.estat
            d2=self.dBdt
            l1=bdname+' Electrostatic'
            l2=bdname+' dB/dt'
        if data is "B":
            t=self.time_B
            d1=self.B
            d2=None
            l1=bdname+' Magnetic Field'
        ax.plot(t, flip*d1, label=l1, lw=4)
        if d2 is not None:
            ax.plot(t, flip*d2, label=l2, lw=4)
        ax.legend()
        
class Bdots:
    def __init__(self, shot, pairs, attenuations, diameters, scope="1", threshold=1.0):
        self.shot=shot
        self.bd={}
        for k, v in  pairs.items():
            bd1=v+"1"
            bd2=v+"2"
            area=(1e-3*diameters[k]/2.0)**2*np.pi
            calibration=[attenuations[bd1]/area, attenuations[bd2]/area]
            self.bd[k]=Bdot_pair(shot, scope, bdot1=bd1, bdot2=bd2)
            self.bd[k].truncate(threshold=threshold,cal=calibration)
            self.bd[k].integrate()
    def plot(self, name, data, ax=None, flip=1):
        self.bd[name].plot(data, ax, flip, bdname=name)
    def plot_raw(self, name):
        self.bd[name].plot_raw()
    def plot_estat_dBdt(self, name):
        self.bd[name].plot_estat_dBdt()
    def plot_B(self, name):
        self.bd[name].plot_B()
            