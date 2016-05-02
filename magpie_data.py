import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import colormaps as cmaps
from scipy.ndimage.interpolation import rotate
from scipy.ndimage import zoom
import os
import image_registration as ir
from skimage.measure import profile_line
import imreg_dft as ird
import images2gif as ig

class DataMap:
    def __init__(self, flip_lr, rot_angle, multiply_by, scale):
        if flip_lr is True:
            self.d=np.fliplr(self.d)
        if rot_angle is not None:
            self.d=rotate(self.d, rot_angle)
        self.data=self.d*multiply_by
        self.scale=scale
    def plot_data_px(self, clim=None, multiply_by=1, ax=None):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.data*multiply_by
        return ax.imshow(d, clim=clim, cmap=self.cmap)
    def set_origin(self, origin, extent):
        self.origin=origin
        ymin=origin[0]-extent[1]*self.scale
        ymax=origin[0]-extent[0]*self.scale
        xmin=origin[1]+extent[2]*self.scale
        xmax=origin[1]+extent[3]*self.scale
        self.origin_crop=(-extent[0]*self.scale,-extent[2]*self.scale)
        self.data_c=self.data[ymin:ymax, xmin:xmax]
        self.extent=extent[2:4]+extent[0:2]
    def plot_data_mm(self, clim=None, multiply_by=1, ax=None):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.data_c*multiply_by
        return ax.imshow(d, cmap=self.cmap, interpolation='none', clim=clim, extent=self.extent, aspect=1)
    def plot_contours_px(self, levels=None, multiply_by=1, ax=None, color='k'):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.data*multiply_by
        return ax.contour(d, levels, origin='image', hold='on',colors=color)
    def plot_contourf_mm(self, levels=None, multiply_by=1, ax=None):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        d=self.data_c*multiply_by
        return ax.contourf(d, levels=levels, cmap=self.cmap,extent=self.extent,origin='image',aspect=1)
    def create_lineout(self, start=(0,0), end=(0,0), lineout_width=20):
        '''
        start and end are in mm on the grid defined by the origin you just set
        '''
        #find coordinates in pixels
        start_px=self.mm_to_px(start)
        end_px=self.mm_to_px(end)
        print(start_px,end_px)
        #use scikit image to do a nice lineout on the cropped array
        self.lo=profile_line(self.data_c, start_px,end_px,linewidth=lineout_width)
        #set up a mm scale centred on 0
        px_range=self.lo.size/2
        self.mm=np.linspace(-px_range, px_range, 2*px_range)/self.scale #flip range to match images
    def plot_lineout(self, ax=None, label='', multiply_by=1):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.plot(self.mm, self.lo*multiply_by, label=label, lw=4)        
    def mm_to_px(self,mm):
        scale=self.scale
        px_origin=self.origin_crop
        return (int(mm[0]*scale+px_origin[0]),int(mm[1]*scale+px_origin[1]))
    
class NeLMap2(DataMap):
    def __init__(self, filename, scale, multiply_by=1, flip_lr=False, rot_angle=None):
        d=np.loadtxt(open(filename,"r"),delimiter=",")
        d=d-np.nan_to_num(d).min()
        d=np.nan_to_num(d)
        if flip_lr is True:
            d=np.fliplr(d)
        if rot_angle is not None:
            d=rotate(d, rot_angle)
        self.data=d*multiply_by
        self.scale=scale
        self.cmap=cmaps.cmaps['inferno']
        
class PolarimetryMap2(DataMap):
    def __init__(self, R0fn, R1fn, B0fn, B1fn, S0fn, S1fn, rot_angle=None):
        self.R0=plt.imread(R0fn)
        self.R1=np.fliplr(plt.imread(R1fn))
        self.B0=plt.imread(B0fn)
        self.B1=np.fliplr(plt.imread(B1fn))
        self.S0=plt.imread(S0fn)
        self.S1=np.fliplr(plt.imread(S1fn))
        if rot_angle is not None:
            self.R0=rotate(self.R0, rot_angle)
            self.R1=rotate(self.R1, rot_angle)
            self.B0=rotate(self.B0, rot_angle)
            self.B1=rotate(self.B1, rot_angle)
            self.S0=rotate(self.S0, rot_angle)
            self.S1=rotate(self.S1, rot_angle)
        #normalise registration images
        R0s=self.R0.sum()
        R1s=self.R1.sum()
        self.R0=self.R0*R1s/R0s
        self.cmap='seismic'
    def register(self):
        self.result=ird.similarity(self.R0, self.R1, numiter=3)
        self.BT=ird.transform_img_dict(self.B1, self.result)
        self.ST=ird.transform_img_dict(self.S1, self.result)
    def convert_to_alpha(self, beta=3.0):
        self.N0=self.S0/self.B0
        self.N1=self.ST/self.BT
        diff=self.N0-self.N1
        self.diff=np.nan_to_num(diff)
        beta=beta*np.pi/180
        self.data=(180/np.pi)*0.5*np.arcsin(self.diff*np.tan(beta)/2.0)
        
class InterferogramOntoAlpha(DataMap):
    def __init__(self, polmap, I0, I1):
        I0=plt.imread(I0)
        self.I0s=np.sum(I0,2)
        self.pm=polmap
        #scale and flip registration to polarisation data
        R0=self.pm.R0
        scale=R0.shape[0]/self.I0s.shape[0]
        I0z=zoom(self.I0s, scale)
        crop=(I0z.shape[1]-R0.shape[1])/2
        I0zc=I0z[:,crop:-crop]
        self.I0zcn=np.flipud(I0zc/I0zc.max())
        #do the same to the inteferogram
        I1=plt.imread(I1)
        I1s=np.sum(I1,2)
        I1z=zoom(I1s, scale)
        I1zc=I1z[:,crop:-crop]
        self.I1zcf=np.flipud(I1zc)
        self.cmap='gray'
    def register(self):
        self.transform=ird.similarity(self.pm.R0, self.I0zcn, numiter=3)
        self.data=ird.transform_img_dict(self.I1zcf, self.transform)
    def plot_overlay_px(self, clim=None, ax=None, transparency=0.8):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.imshow(self.pm.data, cmap='RdBu', clim=clim)
        ax.imshow(self.data, cmap='gray', alpha=transparency)


class FaradayMap2(DataMap):
    def __init__(self, polmap, I0, ne):
        I0=plt.imread(I0)
        self.I0s=np.sum(I0,2)
        I1=np.loadtxt(ne, delimiter=',')
        I1=I1-np.nan_to_num(I1).min()
        self.I1=np.nan_to_num(I1)
        self.pm=polmap
        #scale and flip to data
        B0=self.pm.B0
        scale=B0.shape[0]/self.I0s.shape[0]

        I0z=zoom(self.I0s, scale)
        crop=(I0z.shape[1]-B0.shape[1])/2
        I0zc=I0z[:,crop:-crop]
        self.I0zcn=np.flipud(I0zc/I0zc.max())
        
        I1z=zoom(self.I1, scale)
        self.I1zc=np.flipud(I1z[:,crop:-crop])
        
        self.cmap='seismic'
    def register(self):
        self.transform=ird.similarity(self.pm.R0, self.I0zcn, numiter=3)
        self.IT=ird.transform_img_dict(self.I1zc, self.transform)
        self.data=5.99e18*self.pm.data/self.IT
        
class NeLMap:
    def __init__(self, filename, scale, multiply_by=1, flip_lr=False, rot_angle=None):
        d=np.loadtxt(open(filename,"r"),delimiter=",")
        d=d-np.nan_to_num(d).min()
        d=np.nan_to_num(d)
        if flip_lr is True:
            d=np.fliplr(d)
        if rot_angle is not None:
            d=rotate(d, rot_angle)
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
        y0=y_range*self.scale
        x0=x_range*self.scale
        ymin=origin[0]-y_range*self.scale
        ymax=origin[0]+y_range*self.scale
        xmin=origin[1]-x_range*self.scale
        xmax=origin[1]+x_range*self.scale
        self.origin_crop=[y_range*self.scale,x_range*self.scale]
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
    def create_lineout(self, start=(0,0), end=(0,0), lineout_width=20):
        '''
        start and end are in mm on the grid defined by the origin you just set
        '''
        #find coordinates in pixels
        start_px=mm_to_px(start)
        end_px=mm_to_px(stop)
        #use scikit image to do a nice lineout on the cropped array
        self.lo=profile_line(self.neL_crop, start_px,end_px,linewidth=lineout_width)
        #set up a mm scale centred on 0
        px_range=self.lo.size/2
        self.mm=np.linspace(px_range, -px_range, 2*px_range)/self.scale #flip range to match images
    def plot_lineout(self, ax=None, label=''):
        if ax is None:
            fig, ax=plt.subplots(figsize=(12,8))
        ax.plot(self.mm, self.lo/1e18, label=label, lw=4)        
    def mm_to_px(mm):
        scale=self.scale
        px_origin=self.origin_crop
        return (int(mm[0]*scale+px_origin[0]),int(mm[1]*scale+px_origin[1]))
            
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
    def save_gif(self, filename, clim):
        hot_im=[]
        for im in self.s_c:
            ax.imshow(im, cmap='afmhot', clim=clim)
            plt.axis('off')
            fig.tight_layout()
            fig.canvas.draw()
            w,h=fig.canvas.get_width_height()
            buf=np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8)
            buf.shape=(h,w,3)
            hot_im.append(buf)
        ig.writeGif(filename+'.gif',hot_im, duration=0.2)
        


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
        self.I1zc=np.flipud(I1z[:,crop:-crop])
    def register(self):
        self.I0T, self.I1T, self.scale, self.angle, (self.t0, self.t1)=ir.transform_like(self.pm.B0,self.I0zcn, self.I1zc)
        self.B=5.99e18*self.pm.alpha/self.I1T
        
