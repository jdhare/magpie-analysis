import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import colormaps as cmaps
from scipy.ndimage.interpolation import rotate
from scipy.ndimage import zoom
import os
import image_registration as ir



class NeLMap:
    def __init__(self, filename, scale, multiply_by=1, flip_lr=False):
        d=np.loadtxt(open(filename,"r"),delimiter=",")
        d=d-np.nan_to_num(d).min()
        d=np.nan_to_num(d)
        if flip_lr is True:
            d=np.fliplr(d)
        self.neL=d*multiply_by
        self.scale=scale
    def plot_neL(self, clim=[0,2], ax=None, transpose=False):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.neL/1e18
        if transpose is True:
            d=np.transpose(d)
        ax.imshow(self.neL/1e18)#, cmap='afmhot', interpolation='none', clim=clim)
    def set_origin(self, origin, x_range=11.5, y_range=8.5):
        self.origin=origin
        ymin=origin[0]-y_range*self.scale
        ymax=origin[0]+y_range*self.scale
        xmin=origin[1]-x_range*self.scale
        xmax=origin[1]+x_range*self.scale
        self.neL_crop=self.neL[ymin:ymax, xmin:xmax]
        self.extent=[-x_range,x_range,-y_range,y_range]
    def plot_neL_mm(self,clim=[0,2], ax=None, transpose=False, cmap=cmaps.cmaps['inferno']):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.neL_crop/1e18
        ex=self.extent
        if transpose is True:
            d=np.transpose(d)
            ex=ex[2:4]+ex[0:2]
        return ax.imshow(d, cmap=cmap, interpolation='none', clim=clim, extent=ex, aspect=1)
    def create_lineout(self, axis, centre=None,average_over_px=20, mm_range=10):
        px_range=mm_range*self.scale
        if axis is 1:
            d=np.transpose(self.neL)
            y0=self.origin[1] if centre is None else centre
            x0=self.origin[0]
        if axis is 0:
            d=self.neL
            y0=self.origin[0] if centre is None else centre
            x0=self.origin[1]
        section=d[y0-average_over_px:y0+average_over_px, x0-px_range:x0+px_range]
        self.lo=np.mean(section, axis=0)
        self.mm=np.linspace(px_range, -px_range, 2*px_range)/self.scale #flip range to match images
    def plot_lineout(self, ax=None, label=''):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.plot(self.mm, self.lo/1e18, label=label, lw=4)
        
class OpticalFrames:
    def __init__(self, start, IF):
        self.load_images()
        self.normalise()
        self.start=start
        self.IF=IF
        self.frame_times=np.arange(start, start+12*IF, IF)
    def load_images(self):
        shot=os.path.split(os.getcwd())[-1][0:8] #automatically grab the shot number
        b=[]
        s=[]
        for i in range(1,13):
            if i<10:
                st="0"+str(i)
            else:
                st=str(i) 
            bk_fn=shot+" Background_0"+st+".png"
            bk_im=plt.imread(bk_fn) #read background image
            #bk_im=np.asarray(np.sum(bk_im,2), dtype=float)
            b.append(bk_im)#np.asarray(np.sum(bk_im,2), dtype=float)) #convert to grrayscale
            sh_fn=shot+" Shot_0"+st+".png" 
            sh_im=plt.imread(sh_fn)
            s.append(sh_im)
           
        self.shot=shot
        self.b=b
        self.s=s
    def normalise(self):
        norms=[b_im[100:-100,100:-100].sum() for b_im in self.b]
        n_max=max(norms)
        nn=[n/n_max for n in norms]
        self.s_n=[s_im[100:-100,100:-100]/n for s_im, n in zip(self.s, nn)]
    def logarithm(self, lv_min=-4, lv_max=0.2):
        self.s_l=[np.log(s_im) for s_im in self.s_n]
        self.s_nl=[(np.clip(s_im, a_min=lv_min, a_max=lv_max)-lv_min)/(lv_max-lv_min) for s_im in self.s_l]
    def rotate(self, angle_deg=0):
        self.s_r=[rotate(s_im, angle_deg)for s_im in self.s_nl]
    def crop(self, origin, xcrop=400, ycrop=400):
        x0=origin[1]
        y0=origin[0]
        self.origin=[ycrop,xcrop]
        self.s_c=[s_im[y0-ycrop:y0+ycrop,x0-xcrop:x0+xcrop] for s_im in self.s_r]
    def plot(self, array, frame=1, clim=None, ax=None):
        fin=frame-1
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.imshow(array[fin], cmap='afmhot', clim=clim)
        ax.axis('off')
        ax.set_title('t='+str(self.frame_times[fin])+' ns', fontsize=22)
    def plot_norm(self, frame=1, clim=None, ax=None):
        self.plot(self.s_n, frame=frame, clim=clim, ax=ax)
    def plot_log(self, frame=1, clim=None, ax=None):
        self.plot(self.s_nl, frame=frame, clim=clim, ax=ax)
    def plot_rot(self, frame=1, clim=None, ax=None):
        self.plot(self.s_r, frame=frame, clim=clim, ax=ax)
    def plot_crop(self, frame=1, clim=None, ax=None):
        self.plot(self.s_c, frame=frame, clim=clim, ax=ax)
    def plot_sequence(self, array=None, frames=list(range(1,13)), clim=None, figsize=None):
        xframes=round(len(frames)/2)
        if array is None:
            array=self.s_c
        if figsize is None:
            figsize=(xframes*4,16)
        fig, ax=plt.subplots(2,xframes, figsize=figsize)
        ax=ax.flatten()
        for fin, f in enumerate(frames):
            fn=f-1 #shift to 0 indexed arrays
            a=ax[fin]
            a.imshow(array[fn], cmap='afmhot', clim=clim)
            a.axis('off')
            a.set_title('t='+str(self.frame_times[fn])+' ns', fontsize=22)
        fig.suptitle("Optical Framing images from "+self.shot, fontsize=32)
        fig.tight_layout(w_pad=0, h_pad=0)
        self.fig=fig
    def save_sequence(self, filename=None):
        if filename is None:
            filename=self.shot+" frame sequence"
        self.fig.savefig(filename+".png")        
    def create_lineout(self, axis=0, frame=1,centre=None,average_over_px=20, mm_range=10, scale=29.1, ax=None):
        px_range=mm_range*scale
        fn=frame-1 #shift to 0 indexed arrays
        if axis is 1:
            d=np.transpose(self.s_c[fn])
            y0=self.origin[1] if centre is None else centre
            x0=self.origin[0]
        if axis is 0:
            d=self.s_c[fn]
            y0=self.origin[0] if centre is None else centre
            x0=self.origin[1]
        section=d[y0-average_over_px:y0+average_over_px, x0-px_range:x0+px_range]
        self.lo=np.mean(section, axis=0)
        self.mm=np.linspace(-px_range, px_range, self.lo.size)/scale
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.plot(self.mm, self.lo, label='t='+str(self.frame_times[fn])+' ns', lw=4)
        
def convert_to_alpha(B0,BT,S0,ST, beta):
    N0=S0/B0
    N1=ST/BT
    diff=N0-N1
    diff=np.nan_to_num(diff)
    beta=beta*np.pi/180
    alpha=(180/np.pi)*0.5*np.arcsin(diff*np.tan(beta)/2.0)
    return alpha

class PolarimetryMap:
    def __init__(self, B0fn, B1fn, S0fn, S1fn):
        self.B0=plt.imread(B0fn)
        self.B1=plt.imread(B1fn)
        self.S0=plt.imread(S0fn)
        self.S1=plt.imread(S1fn)
    def register(self):
        self.BT, self.ST,self.scale, self.angle, (self.t0, self.t1)=ir.transform_like(self.B0,self.B1, self.S1)
    def convert_to_alpha(self, beta=3.0):
        self.N0=self.S0/self.B0
        self.N1=self.ST/self.BT
        diff=self.N0-self.N1
        self.diff=np.nan_to_num(diff)
        beta=beta*np.pi/180
        self.alpha=(180/np.pi)*0.5*np.arcsin(self.diff*np.tan(beta)/2.0)
    def set_origin(self, origin, x_range=11.5, y_range=8.5):
        self.origin=origin
        ymin=origin[0]-y_range*self.scale
        ymax=origin[0]+y_range*self.scale
        xmin=origin[1]-x_range*self.scale
        xmax=origin[1]+x_range*self.scale
        self.alpha_crop=self.alpha[ymin:ymax, xmin:xmax]
        self.extent=[-x_range,x_range,-y_range,y_range]
    def plot_alpha(self,clim=[-2,2], ax=None, transpose=False, cmap=plt.cm.seismic):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.alpha
        if transpose is True:
            d=np.transpose(d)
        return ax.imshow(d, cmap=cmap, interpolation='none', clim=clim, aspect=1)
    def plot_alpha_mm(self,clim=[-2,2], ax=None, transpose=False, cmap=plt.cm.seismic):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.alpha_crop
        ex=self.extent
        if transpose is True:
            d=np.transpose(d)
            ex=ex[2:4]+ex[0:2]
        return ax.imshow(d, cmap=cmap, interpolation='none', clim=clim, extent=ex, aspect=1)
    
class FaradayMap:
    def __init__(self, polmap, I0, ne):
        I0=plt.imread(I0)
        self.I0s=np.sum(I0,2)
        I1=np.loadtxt(ne, delimiter=',')
        self.I1=np.nan_to_num(I1)
        self.pm=polmap
    def scale_and_crop(self):
        B0=self.pm.B0
        scale=B0.shape[0]/self.I0s.shape[0]
        I0z=zoom(self.I0s, scale)
        crop=(I0z.shape[1]-B0.shape[1])/2
        I0zc=I0z[:,crop:-crop]
        self.I0zcn=np.fliplr(I0zc/I0zc.max())
        
        I1z=zoom(self.I1, scale)
        I1zc=I1z[:,crop:-crop]
        self.I1zc=np.fliplr(I1z[:,crop:-crop])
    def register(self):
        self.I0T, self.I1T, self.scale, self.angle, (self.t0, self.t1)=ir.transform_like(self.pm.B0,self.I0zcn, self.I1zc)
        self.B=5.99e18*self.pm.alpha/self.I1T