from spectral_density_functions import *
import csv
from ipywidgets import interact, fixed
import matplotlib.pyplot as plt
from lmfit import Model
from lmfit.models import VoigtModel
import os


class Fibre:
    def __init__(self, wavelength, background, shot, bkgd_ferr, shot_ferr, theta):
        self.lamb=wavelength*1e-9
        self.bkgd=background
        self.shot=shot
        self.shot_ferr=shot_ferr
        self.bkgd_ferr=bkgd_ferr
        self.theta=theta
    def voigt_response(self, sigma=None, gamma=None):
        '''
        Fit the background with a Voigt profile to determine the response
        of the spectrometer

        If you have a good, clear signal, set sigma and gamma to None (done by default)

        If your signal is poor, set sigma and gamma using a fit to a good signal, and then
        only the position of the central wavelength will be altered.
        '''
        vm=VoigtModel()
        par_v=vm.guess(self.bkgd, x=self.lamb)
        par_v['center'].set(value=532e-9, vary=True)
        err=(self.bkgd_ferr*self.shot)
        if sigma is not None: #if a width is provided, fix it.
            par_v['sigma'].set(value=sigma, vary=False)
        if gamma is not None: #if a width is provided, fix it.
            par_v['gamma'].set(value=gamma, vary=False, expr='')
        elif gamma is None: #vary gamma for better fit - this is not done by default
            par_v['gamma'].set(value=par_v['sigma'].value,vary=True, expr='')

        ##Fit the Voigt Model to the data
        vm_fit=vm.fit(self.bkgd,par_v,x=self.lamb)
        self.vm_fit=vm_fit
        #now crop the data so that the response is symmetric for the convolution to work
        l0=vm_fit.best_values['center']
        self.sigma=vm_fit.best_values['sigma']
        self.l0=l0
        l0_i=find_nearest(self.lamb, l0)
        l_size=self.lamb.size
        take_l=min(l0_i,l_size-l0_i) #trim the shortest distance from the central wavelength
        low_i=l0_i-take_l
        high_i=l0_i+take_l
        self.lamb=self.lamb[low_i:high_i]
        self.bkgd=self.bkgd[low_i:high_i]
        self.shot=self.shot[low_i:high_i]
        self.shot_ferr=self.shot_ferr[low_i:high_i]
        self.bkgd_ferr=self.bkgd_ferr[low_i:high_i]
        #the response is taken from the model so it is nice and smooth
        self.response=vm_fit.best_fit[low_i:high_i]
        self.shift=self.lamb-l0 #this is useful for plotting data
    def fit_fibre(self, pp):
        '''
        Fit the shot data. This is complicated!
        This examines the dictionary provided, determines which are dependent and independent variables
        It then chooses the correct model, and sets up lmfit.
        '''
        pp_valid={}
        iv_dict={} #dictionary for independent variables

        '''Choose the model'''
        if pp['model']=='two stream':
            valid_keys=['model','n_e','T_e','V_fe','Aj','Zj','Fj','T_i1','V_fi1','T_i2','V_fi2','stry','amp', 'offset']
            for k in valid_keys:
                pp_valid[k]=pp[k]
            skw_func=Skw_two_stream_stray_light_convolve
        elif pp['model']=='multi species':
            valid_keys=['model','n_e','T_e','V_fe','Aj','Zj','Fj','T_i1','V_fi1','stry','amp', 'offset']
            for k in valid_keys:
                pp_valid[k]=pp[k]
            skw_func=Skw_multi_species_stray_light_convolve
        elif pp['model']=='nLTE':
            valid_keys=['model','n_e','T_e','V_fe','A','T_i','V_fi','stry','amp', 'offset', 'shift']
            for k in valid_keys:
                pp_valid[k]=pp[k]
            skw_func=Skw_nLTE_stray_light_convolve
            if pp['A'][0][0]==12:
                iv_dict['Z_Te_table']=np.genfromtxt('zb_Spk_C_NLTE.txt', delimiter='       ', skip_header=4)
            elif pp['A'][0][0]==27:
                iv_dict['Z_Te_table']=np.genfromtxt('zb_Al.dat', delimiter=' ', skip_header=2)
            elif pp['A'][0][0]==183:
                iv_dict['Z_Te_table']=np.genfromtxt('zb_W.dat', delimiter=' ')
        elif pp['model']=='Collisional nLTE':
            print('Not a suitable model in most cases, use nLTE instead')
            valid_keys=['model','n_e','T_e','V_fe','Aj','T_i1','V_fi1','stry','amp', 'offset', 'shift']
            for k in valid_keys:
                pp_valid[k]=pp[k]
            skw_func=Skw_nu_nLTE_stray_light_convolve
            if pp['Aj'][0][0]==12:
                iv_dict['Z_Te_table']=np.genfromtxt('zb_Spk_C_NLTE.txt', delimiter='       ', skip_header=4)
            elif pp['Aj'][0][0]==27:
                iv_dict['Z_Te_table']=np.genfromtxt('zb_Al.dat', delimiter=' ', skip_header=2)


        #sort into dictionaries based on another dictionary
        for k,v in pp_valid.items():
            if v[1] is True: #independent variable
                iv_dict[k]=pp_valid[k][0]#get first element of list
        iv_dict['lambda_range']=self.lamb
        iv_dict['lambda_in']=self.l0
        iv_dict['response']=self.response
        iv_dict['theta']=self.theta
        skw=Model(skw_func, independent_vars=list(iv_dict.keys())) #create our model with our set variables.
        #our best guesses at what the fitting parameters should be
        for k,v in pp_valid.items():
            if v[1] is False: #dependent variable
                try:
                    skw.set_param_hint(k, value = v[0], min=v[2]) #if a minimum is provided, use it
                except IndexError:
                    skw.set_param_hint(k, value = v[0])

        '''now do the fitting'''
        skw_res=skw.fit(self.shot,verbose=False, **iv_dict)
        self.skw_res=skw_res
        ##Calculating alpha. This is really ugly
        try:
            te=skw_res.best_values['T_e']
        except KeyError: #if this is an independent variable it isn't in the fit
            te=pp_valid['T_e'][0]
        try:
            ti=skw_res.best_values['T_i']
        except KeyError: #if this is an independent variable it isn't in the fit
            ti=pp_valid['T_i'][0]
        ne=pp_valid['n_e'][0]
        lambda_De=7.43*(te/ne)**0.5 #in m
        k=4*np.pi*np.sin(self.theta/2.0)/self.l0
        alpha=np.abs(1/(k*lambda_De))

        #create a dictionary which contains the fit parameters
        params=skw_res.best_values.copy()
        for k,v in iv_dict.items():
                params[k]=v
        params['alpha']=alpha
        params['model']=pp_valid['model']
        #calculate the expected TS intensity
        if pp['model']=='nLTE':
            params['Z']=Z_nLTE(te, iv_dict['Z_Te_table'])
            Z=params['Z']
            TS_norm=Z*ne*alpha**4/((1+alpha**2)*(1+alpha**2+alpha**2*Z*te/ti))
            params['TS_norm']=TS_norm
        if pp['model']=='Collisional nLTE':
            params['Z']=Z_nLTE(te, iv_dict['Z_Te_table'])
            Z=params['Z']
            TS_norm=Z*ne*alpha**4/((1+alpha**2)*(1+alpha**2+alpha**2*Z*te/ti))
            params['TS_norm']=TS_norm
        remove_keys=['lambda_range','lambda_in','response','theta'] #remove pointless keys
        [params.pop(k, None) for k in remove_keys]
        self.fit_parameters=params
    def export_data(self, filename):
        data=list(zip(self.shift*1e10, self.bkgd,self.response, self.shot, self.skw_res.best_fit))
        headings=('Wavelength shift', 'Background', 'Response', 'Shot', 'Fit')
        units=('Angstroms', 'a.u.', 'a.u.', 'a.u.', 'a.u.')
        with open(filename+'.dat', 'w',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headings)
            writer.writerow(units)
            for f in data:
                writer.writerow(f)

class TS_Analysis:
    def __init__(self, folder, shot, backgrounds, calibration=False, skip_footer=0):
        '''
        Loads .asc files.
        shot should be the path to a single.asc
        backgrounds is a list of paths to multiple backgrounds - useful if you want to sum several together
        calibration is if your .asc files don't have the wavelength calibration embdeed,
        you can borrow it from another. .asc.
        skip_footer is if the end of your file is corrupted I guess? But maybe don't use it.
        '''
        ts_dir=os.getcwd()
        os.chdir(folder)
        self.s_name=os.path.basename(shot)[:8]
        data=np.genfromtxt(open(shot,"rb"),delimiter="\t", skip_footer=skip_footer)
        self.shot=np.rot90(data[:,1:-1]) #shot data
        self.x_axis=data[:,0]#x_axis, either pixel number or wavelength scale
        #Create an empty array to load multiple backgrounds into
        self.background=np.zeros((self.shot.shape))
        #load multiple backgrounds and sum them together
        for b in backgrounds:
            data=np.genfromtxt(open(b,"rb"),delimiter="\t", skip_footer=skip_footer)
            bkgd=np.rot90(data[:,1:-1])
            self.background+=bkgd
        #If a file is provided with a seperate calibration, use that calibration
        if calibration:
            data=np.genfromtxt(open(calibration,"rb"),delimiter="\t")
            self.x_axis=data[:,0]
        os.chdir(ts_dir)

    def plot_fibre_edges(self, spacing=19.0, offset=1.0):
        '''
        Plots the intensity of the signal on each fibre.
        And also little red dots at where the fiber edges are currently meant to be.
        '''
        pe=self.shot[:,950:1100].sum(1)
        self.fig_edge, self.ax_edge=plt.subplots(figsize=(15,4))
        self.plot_edge=self.ax_edge.plot(pe)
        fe=np.arange(offset, pe.size-1, spacing)
        fe=np.append(fe,pe.size-1)
        fe=np.around(fe)
        fe=fe.astype(int)#need integers to index arrays
        #can't have the first bin start before the image
        if fe[0]<0:
            fe[0]=0
        self.fibre_edges=fe
        self.plot_dots=self.ax_edge.plot(fe, pe[fe], 'r.')
    def find_fibre_edges(self):
        interact(self.plot_fibre_edges, spacing=(10,30,0.1),offset=(0,512,1))
    def split_into_fibres(self, discard_rows=3, numA=14):
        '''Splits the images into 1D arrays for each fibre'''
        fe=np.array(self.fibre_edges)
        self.N_fibres=fe.size
        ##Shot Fibres
        shot_fibres=np.zeros((fe.size-1, self.shot.shape[1]))
        #take the data in each fibre, discard some rows where there is cross talk
        #and then sum together in the y direction to improve signal to noise.
        for i in np.arange(fe.size-1):
            shot_fibres[i]=self.shot[fe[i]+discard_rows:fe[i+1]-discard_rows,:].sum(0)
        self.shot_fibres=shot_fibres
        self.shot_frac_err=1/np.sqrt(shot_fibres)#assuming poisson stats, sigma=sqrt(y), so sigma/y=1/sqrt(y)
        ##Background Fibres
        bkgd_fibres=np.zeros((fe.size-1, self.background.shape[1]))
        for i in np.arange(fe.size-1):
            bkgd_fibres[i]=self.background[fe[i]:fe[i+1],:].sum(0)
        self.bkgd_fibres=bkgd_fibres
        self.bkgd_frac_err=1/np.sqrt(bkgd_fibres)#assuming poisson stats, sigma=sqrt(y), so sigma/y=1/sqrt(y)
    def zero_fibres(self, lower=500, upper=1500):
        self.shot_fibres_z=np.zeros((self.N_fibres,upper-lower))
        self.bkgd_fibres_z=np.zeros((self.N_fibres,upper-lower))
        self.S_T=np.zeros(self.N_fibres)
        for fin, f in enumerate(self.shot_fibres):
            #remove offset due to dark counts
            mean1=f[0:300].mean()
            mean2=f[-300:].mean()
            mean=(mean1+mean2)/2
            f=f-mean #zero the fibres
            self.S_T[fin]=np.trapz(f, x=self.x_axis) #calculate the total scattered amplitude
            f=f[lower:upper]
            self.shot_fibres_z[fin]=f
        for fin, f in enumerate(self.bkgd_fibres):
            mean1=f[0:300].mean()
            mean2=f[-300:].mean()
            mean=(mean1+mean2)/2
            f=f-mean #zero the fibres
            f=f[lower:upper]
            self.bkgd_fibres_z[fin]=f
        self.shot_frac_err=self.shot_frac_err[:,lower:upper]
        self.bkgd_frac_err=self.bkgd_frac_err[:,lower:upper]
        self.x_axis=self.x_axis[lower:upper]
    def pair_fibres(self, angle_a, angle_b):
        fibre_angles=angle_a+angle_b
        fibres=[]
        l=self.x_axis
        params=list(zip(self.bkgd_fibres_z,self.shot_fibres_z, self.bkgd_frac_err, self.shot_frac_err, fibre_angles))
        self.fibres=[Fibre(l,bkgd,shot,bkgd_err,shot_err, angle) for bkgd,shot,bkgd_err,shot_err, angle in params]
        self.fibres_a=self.fibres[:len(angle_a)]
        self.fibres_b=self.fibres[len(angle_a):]
    def copy_background(self, good, bad):
        self.fibres[bad].bkgd=self.fibres[good].bkgd.copy() #copy a good background to overwrite a bad one
    def select_fibre(self, Fnum, Fset):
        if Fset=='A':
            f=self.fibres_a[Fnum-1]
        elif Fset=='B':
            f=self.fibres_b[Fnum-1]
        return f
    def plot_data(self, Fnum, Fset, sr=8, tm=1.0):
        '''Probably the prettiest plot you've ever seen'''
        f=self.select_fibre(Fnum,Fset)
        text_mul=tm
        fig, ax=plt.subplots(figsize=(16,10))
        bk_norm=0.5*f.shot.max()/f.bkgd.max()
        plot_data=ax.plot(f.lamb,bk_norm*f.bkgd, label='Background', lw=2, marker='o', c='0.5')
        plot_data=ax.plot(f.lamb,f.shot,label='Data', marker='o',lw=2, c='b')
        #plotting region
        ax.set_ylim(bottom=0.0)
        ax.set_xlim([531e-9,533e-9])
        ax.set_xlabel(r'Wavelength (nm)',fontsize=20*text_mul)
        ax.set_ylabel('Intensity (a.u.)',fontsize=20*text_mul)
        ax.tick_params(labelsize=20*text_mul, pad=5, length=10, width=2)
        title_str=self.s_name+': Thomson Scattering for fibre '+str(Fnum)+Fset+r', $\theta=$'+str(f.theta)+r'$^{\circ}$'
        ax.set_title(title_str,fontsize=20*text_mul)
        ax.legend(fontsize=18*text_mul)
        plt.tight_layout()
        self.fig=fig
        self.ax=ax
    def pretty_plot(self, Fnum, Fset ,sr=8, tm=1.0):
        '''Probably the prettiest plot you've ever seen'''
        f=self.select_fibre(Fnum,Fset)
        try:
            shift=f.fit_parameters['shift']
        except KeyError:
            shift=0
        response=np.interp(f.lamb+shift, f.lamb, f.response)
        text_mul=tm

        bk_norm=0.5*f.shot.max()/f.bkgd.max()
        fig, ax=plt.subplots(figsize=(16,10))
        plot_data=ax.plot(f.shift*1e10,bk_norm*f.bkgd, label='Background', lw=2, marker='o', color='gray')
        plot_data=ax.plot(f.shift*1e10,bk_norm*response, label='Response', lw=2, ls='--', color='black')
        plot_data=ax.plot(f.shift*1e10,f.shot,label='Data', marker='o',lw=2, color='blue')
        plot_fit=ax.plot(f.shift*1e10,f.skw_res.best_fit, label='Best Fit', lw=3, ls='--', color='red')
        #plotting region
        ax.set_ylim(bottom=0.0)
        ax.set_xlim([-sr,sr])
        ax.set_xticks(np.arange(-sr,sr+1,2))
        ax.set_xlabel(r'Wavelength shift, $(\AA)$',fontsize=20*text_mul)
        ax.set_ylabel('Intensity (a.u.)',fontsize=20*text_mul)
        ax.tick_params(labelsize=20*text_mul, pad=5, length=10, width=2)
        kms=r' $km\,s^{-1}$'
        if f.fit_parameters['model']=='multi species':
            string_list=[
                    r'$F\,= $'+str(f.fit_parameters['Fj']),
                    r'$A\,= $'+str(f.fit_parameters['Aj']),
                    r'$Z\,= $'+str(f.fit_parameters['Zj']),
                    r'$n_e= $'+str_to_n(f.fit_parameters['n_e']/1e17,2)+r'$\times$10$^{17} cm^{-3}$',
                    r'$T_e= $'+str_to_n(f.fit_parameters['T_e'],2)+' $eV$',
                    r'$T_i= $'+str_to_n(f.fit_parameters['T_i1'],2)+' $eV$',
                    r'$V_{fi}= $'+str_to_n(f.fit_parameters['V_fi1']/1e3,2)+kms,
                    r'$V_{fe}= $'+str_to_n(f.fit_parameters['V_fe']/1e3,2)+kms,
                    r'$\alpha\,= $'+str_to_n(f.fit_parameters['alpha'],2),
                    ]
        if f.fit_parameters['model']=='two stream':
            string_list=[
                r'$F\,= $'+str(f.fit_parameters['Fj']),
                r'$A\,= $'+str(f.fit_parameters['Aj']),
                r'$Z\,= $'+str(f.fit_parameters['Zj']),
                r'$n_e= $'+str_to_n(f.fit_parameters['n_e']/1e17,2)+r'$\times$10$^{17} cm^{-3}$',
                r'$T_e= $'+str_to_n(f.fit_parameters['T_e'],2)+' $eV$',
                r'$T_{i,1}= $'+str_to_n(f.fit_parameters['T_i1'],2)+' $eV$',
                r'$T_{i,2}= $'+str_to_n(f.fit_parameters['T_i2'],2)+' $eV$',
                r'$V_{fi,1}= $'+str_to_n(f.fit_parameters['V_fi1']/1e3,2)+kms,
                r'$V_{fi,2}= $'+str_to_n(f.fit_parameters['V_fi2']/1e3,2)+kms,
                r'$V_{fe}= $'+str_to_n(f.fit_parameters['V_fe']/1e3,2)+kms,
                r'$\alpha\,= $'+str_to_n(f.fit_parameters['alpha'],2),
                ]
        if f.fit_parameters['model']=='nLTE':
            string_list=[
                    r'$A\,= $'+str(f.fit_parameters['A'][0]),
                    r'$Z\,= $'+str_to_n(f.fit_parameters['Z'],2),
                    r'$n_e= $'+str_to_n(f.fit_parameters['n_e']/1e17,2)+r'$\times$10$^{17} cm^{-3}$',
                    r'$T_e= $'+str_to_n(f.fit_parameters['T_e'],2)+' $eV$',
                    r'$T_i= $'+str_to_n(f.fit_parameters['T_i'],2)+' $eV$',
                    r'$V_{fi}= $'+str_to_n(f.fit_parameters['V_fi']/1e3,2)+kms,
                    r'$V_{fe}= $'+str_to_n(f.fit_parameters['V_fe']/1e3,2)+kms,
                    r'$\alpha\,= $'+str_to_n(f.fit_parameters['alpha'],2),
                    ]
        if f.fit_parameters['model']=='Collisional nLTE':
            string_list=[
                    r'$A\,= $'+str(f.fit_parameters['Aj'][0]),
                    r'$Z\,= $'+str_to_n(f.fit_parameters['Z'],2),
                    r'$n_e= $'+str_to_n(f.fit_parameters['n_e']/1e17,2)+r'$\times$10$^{17} cm^{-3}$',
                    r'$T_e= $'+str_to_n(f.fit_parameters['T_e'],2)+' $eV$',
                    r'$T_i= $'+str_to_n(f.fit_parameters['T_i1'],2)+' $eV$',
                    r'$V_{fi}= $'+str_to_n(f.fit_parameters['V_fi1']/1e3,2)+kms,
                    r'$V_{fe}= $'+str_to_n(f.fit_parameters['V_fe']/1e3,2)+kms,
                    r'$\alpha\,= $'+str_to_n(f.fit_parameters['alpha'],2),
                    ]

        text_str=''
        for st in string_list:
            text_str=text_str+st+'\n'
        text_str=text_str[:-1]

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='gray', alpha=0.2)

        # place a text box in upper left in axes coords
        ax.text(0.02, 0.96, text_str, transform=ax.transAxes, fontsize=20*text_mul,
            verticalalignment='top', bbox=props)

        title_str=self.s_name+': Fit of Thomson Scattering for fibre '+str(Fnum)+Fset+r', $\theta=$'+str(f.theta)+r'$^{\circ}$'
        ax.set_title(title_str,fontsize=20*text_mul)
        ax.legend(fontsize=18*text_mul)
        plt.tight_layout()
        self.fig=fig
        self.ax=ax
    def export_data(self, Fnum, Fset):
        f=self.select_fibre(Fnum,Fset)
        filename=self.s_name+' fit dat files/'+self.s_name+'_'+str(Fnum)+Fset+'_data_and_fit'
        f.export_data(filename)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def ZTe_finder(n_e, ZTe_experimental, Z_guess, element='Al'):
    if element=='Al':
        Z_Te_table=np.genfromtxt('zb_Al.dat', delimiter=' ', skip_header=2)
        ni=np.array([1e17,5e17,1e18,5e18,1e19])
    if element=='C':
        Z_Te_table=np.genfromtxt('zb_C.dat')
        ni=np.array([1e19])#always choose lowest for now
    ni_guess=n_e/float(Z_guess)
    ind=find_nearest(ni, ni_guess)+1
    ZTe_list=Z_Te_table[:,0]*Z_Te_table[:,ind]
    index=np.where(ZTe_list>=ZTe_experimental)[0][0]
    Te=Z_Te_table[index,0]
    Z=Z_Te_table[index,ind]
    return Z, Te

def Z_finder(n_e, Te_experimental, Z_guess=4, element='Al'):
    if element=='Al':
        Z_Te_table=np.genfromtxt('zb_Al.dat', delimiter=' ', skip_header=2)
        ni=np.array([1e17,5e17,1e18,5e18,1e19])
    if element=='C':
        Z_Te_table=np.genfromtxt('zb_C.dat')
        ni=np.array([1e19])#always choose lowest for now
    ni_guess=n_e/float(Z_guess)
    ind=find_nearest(ni, ni_guess)+1
    Te_list=Z_Te_table[:,0]
    index=np.where(Te_list>=Te_experimental)[0][0]
    Z=Z_Te_table[index,ind]
    return Z
