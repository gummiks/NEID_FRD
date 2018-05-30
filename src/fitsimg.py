from __future__ import print_function
import glob
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch, SqrtStretch, AsinhStretch, HistEqStretch,ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import phothelp
import gkfit
import numpy as np
import gklib as gk
import butter
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import filepath
import matplotlib.gridspec as gridspec
import radial_data



class FitsImage(object):
    """
    A helper class when reading fits files. Depends on pyFits.

    """
    DIRLOC = ''
    
    def __init__(self,filename):
        self.filename = filename
        self.hdulist = fits.open(self.filename)
        #print self.hdulist.info()
        self.header = self.hdulist[0].header
        data = self.hdulist[0].data
        self.data = data.astype(float)
    
    def resize(self, shape):
        """
        A function to resize the image data to a shape, preserving the range.
        
        INPUT:
            shape - tuple
        
        EXAMPLE:
            f = FitsImage("fits_10_blue_red_ir_1500/0001_266.6147_27.72068_dss2blue_N333.fits")
            f.resize((100,100))
            f.plot()
            f.data.shape
        """
        self.data = resize(self.data,shape,preserve_range=True)

    def get_radial_profile(self,rmax=None,plot=False,z=20.,return_hwzm=False,ax=None,xcen=None,ycen=None,annulus_width=1,subtract_min=False):
        """
        Plot radial profile
        """
        if not hasattr(self,"radial"):
            print("Calculating radial data")
            self.radial = radial_data.radial_data(self.data,rmax=rmax,x=xcen,y=ycen,annulus_width=annulus_width)
        else:
            print("Warning: using stored radial data")
        if subtract_min:
            print("Subtracting min value from azimuthal average", np.min(self.radial.mean))
            self.radial.mean -= np.min(self.radial.mean)
        self.radial_hwhm = gkfit.calc_hwhm(self.radial.r,self.radial.mean)
        self.radial_hwzm = gkfit.calc_hwzm(self.radial.r,self.radial.mean,z=z)
        z = float(z)

        print("HWZM",self.radial_hwzm)
        
        if plot:
            if ax==None:
                self.fig, self.ax = plt.subplots()
            else:
                self.ax = ax
            self.ax.plot(self.radial.r,self.radial.mean)
            ymin, ymax = self.ax.get_ylim()
            self.ax.vlines(self.radial_hwhm,ymin,ymax,label="HWHM="+gk.num2str(self.radial_hwhm),color="orange",linestyle="--",lw=1)
            self.ax.vlines(self.radial_hwzm,ymin,ymax,label="HWZM(z="+gk.num2str(z)+")="+gk.num2str(self.radial_hwzm),color="red",linestyle="--",lw=1)
            self.ax.legend(loc="upper right",fontsize=14)
            self.ax.grid(lw=0.5,alpha=0.3)


        if return_hwzm:
            return self.radial.r,self.radial.mean,self.radial_hwzm
        else:
            return self.radial.r,self.radial.mean

    def get_ee_profile(self,rmax=None,plot=False,ax=None,xcen=None,ycen=None,annulus_width=1):
        if not hasattr(self,"radial"):
            self.radial = radial_data.radial_data(self.data,rmax=rmax,x=xcen,y=ycen,annulus_width=1)
        self.radial_hwhm = gkfit.calc_hwhm(self.radial.r,self.radial.mean)

        if plot:
            if ax==None:
                self.fig, self.ax = plt.subplots()
            else:
                self.ax = ax
            self.ax.plot(self.radial.r,self.radial.sum)
            ymin, ymax = self.ax.get_ylim()
            #self.ax.vlines(self.radial_hwhm,ymin,ymax,label="HWHM="+gk.num2str(self.radial_hwhm),color="orange",linestyle="--",lw=1)
            #self.ax.vlines(self.radial_hwzm,ymin,ymax,label="HWZM(z="+gk.num2str(z)+")="+gk.num2str(self.radial_hwzm),color="red",linestyle="--",lw=1)
            self.ax.legend(loc="upper left",fontsize=14)
            self.ax.grid(lw=0.5,alpha=0.3)





    def get_centroid(self,plot_cross=False,ax=None,plot_lines=False):
        """
        Find centroid using Howell centroiding.

        See phothelp for the method
        """
        self.xcenter, self.ycenter = phothelp.howell_center(self.data)
        if plot_cross:
            self.plot(ax=ax)
            self.ax.scatter(self.xcenter,self.ycenter,marker="+",s=50,color="green")
        if plot_lines:
            self.plot(ax=ax)
            self.ax.hlines(int(self.ycenter),0,self.data.shape[1],color='#1f77b4',lw=1)
            self.ax.vlines(int(self.xcenter),0,self.data.shape[0],color='#1f77b4',lw=1)
        return self.xcenter, self.ycenter

    def get_centroid_line_cut(self,line="X",plot=False,ax=None,return_FWHM=False,use_butter_filter=False,butter_cutoff_freq=0.03,fwzm_z = 10):
        """
        A convenience method to get a horizontal (line=="x") or a vertical cut (line=="y") at the centroid of the image.
        This is useful for fitting those curves.

        INPUT:
         line="X"  - get the horizontal cut at the centroid
         line="Y"  - get the vertical cut at the centroid
         plot=True - plot the  cut on an axis *ax*
         butter=True - use a butter low pass filter to help with return_FWHM

        OUTPUT:
         cut - the horizontal / vertical cut that goes through the centroid of the image

        NOTES:
        """
        x, y = self.get_centroid()
        x, y = int(x), int(y)
        if line=="Y": # vertical cut
            cut = self.data[:,x]
        elif line=="X": # horizontal cut
            cut = self.data[y,:]
        else:
            raise ValueError("Line must be equal to 'X' (a single row) or 'Y' (a single column)")

        
        self.cut_x = np.arange(0,len(cut))
        self.cut_y = cut

        if return_FWHM:
            if use_butter_filter:
                self.fwhm = gkfit.calc_fwzm(self.cut_x,butter.low_pass_butter(self.cut_y,cutoff_freq=butter_cutoff_freq),z=fwzm_z)
            else:
                self.fwhm = gkfit.calc_fwzm(self.cut_x,self.cut_y,z=fwzm_z)

        if plot:
            if ax==None:
                self.fig, self.ax = plt.subplots()
            else:
                self.ax = ax
            self.ax.plot(self.cut_x,self.cut_y,label=line+" cut")
            self.ax.minorticks_on()
            self.ax.set_xlabel("X")
            self.ax.set_ylabel("Counts")
            if return_FWHM:
                height = (self.cut_y.max()-self.cut_y.min())/fwzm_z
                self.ax.hlines(height,x-self.fwhm/2.,x+self.fwhm/2.,color="orange",lw=1,label=line+" FWHM="+gk.num2str(self.fwhm))
            self.ax.legend(loc="upper right",fontsize=12)

        if return_FWHM:
            return cut, self.fwhm
        else:
            return cut

    def crop(self,x,y,w,h):
        """
        Crop to a box centered at (x,y), of size w x h
        """
        x, y = int(x),int(y)
        self.data = self.data[y-h/2:y+h/2,x-w/2:x+w/2]
        
    def cropcenter(self,w,h,points=False):
        """
        Returns an image array around the center of an image array.
        """
        shape = self.data.shape
        x,y = (shape[0] / 2, shape[1] / 2)
        self.crop(x,y,w,h)

    def cropcentroid(self,w,h):
        x,y = self.get_centroid()
        self.crop(x,y,w,h)

        
    def plot(self,stretch="hist",cmap="gray",origin="lower",ax=None,colorbar=False,title=""):
        if ax == None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = ax
        if stretch=="hist":
            norm = ImageNormalize(stretch=HistEqStretch(self.data))
            self.im = self.ax.imshow(self.data,cmap=cmap,origin=origin,norm=norm)
        else:
            self.im = self.ax.imshow(self.data,cmap=cmap,origin=origin)
        self.ax.set_xlim(0,self.data.shape[1]) # cols
        self.ax.set_ylim(0,self.data.shape[0]) # rows
        self.ax.set_title(title,y=1.02)
        self.ax.set_xlabel("X pixels")
        self.ax.set_ylabel("Y pixels")
        if colorbar:
            self.fig.colorbar(self.im)
        
    def hist(self,num=1000,figsize=(20,8)):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=1,figsize=figsize)
        self.axes.hist(self.data.flat,num)
        
    def plot_hist(self,cmap='gray',figsize=(20,8),vmin=None,vmax=None,scaled=False,num=1000,titlepre=""):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=2,figsize=figsize)
        im = self.axes[0].imshow(self.data,cmap=cmap,vmin=vmin,vmax=vmax)
        #self.axes[0].set_title(titlepre + "DATE:" + self.header["DATE-OBS"] + " EXP: " + str(self.header["EXPOSURE"]) + " SHAPE: " +str(self.data.shape))
        self.fig.colorbar(im)
        self.axes[1].hist(self.data.flat,num,log=True)


    def plot_3D(self,cut_levels=None,stride=3,surface_lw=0.1,surface_alpha=0.96,ax=None,azim=None,elev=None):
        """
        Plot a 3D layout of the image

        INPUT:
            cut_levels - a list of levels to cut. This is in pixel coordinates.
            stride     - the stride of the grid to plot on the surface. Lower is finer grid
            surface_lw - line width of the surface
            surface_alpha - alpha of surface
            ax - optional axis object
            azim - set the azimuth camera viewing angle
            elev - set the elevation of the camera viewing angle

        OUTPUT:

        EXAMPLE:
            fimg= fitsimg.FitsImage("/Users/gks/Dropbox/mypylib/notebooks/PAPERS/DIFFUSER_PAPER/data/on_sky_psf/16Cyg_semrock_diffuser_rot.0065_out_crop_400_400.fits")
            fimg.cropcentroid(w=200,h=200)
            fimg.plot_3D(cut_levels=[100],stride=3,surface_lw=0.1)
        """
        # Define figure
        if ax==None:
            self.fig = plt.figure(figsize=(20,12))
            self.ax  = plt.axes(projection="3d")
        else:
            self.ax = ax
        
        # Define the data
        X = range(0,self.data.shape[0])
        Y = range(0,self.data.shape[1])
        Z = self.data
        xx, yy = np.meshgrid(X,Y)

        # Plot
        print("Surface_linewidth",surface_lw)
        cmap = sns.cubehelix_palette(as_cmap=True, start=0.45,dark=0.3, light=1., reverse=True,rot=-0.75)
        surf = self.ax.plot_surface(xx,yy,Z,cmap=cmap,linewidth=surface_lw,alpha=surface_alpha,rstride=stride,cstride=stride,edgecolors='k')

        # Plot cuts
        if cut_levels!=None:
            cmap_green = sns.cubehelix_palette(as_cmap=True, start=0.6,dark=0.48, light=1., reverse=True,rot=-0.9)
            cset = self.ax.contour(xx, yy, Z, zdir='x', offset=0, cmap=cmap_green,levels=cut_levels)
            cset = self.ax.contour(xx, yy, Z, zdir='y', offset=200, cmap=cmap_green,levels=cut_levels)

        self.ax.set_frame_on(True)
        for i in self.ax.get_xgridlines():
            i.set_alpha(0.1)
            i.set_lw(0.5)
            
        self.ax.set_xlabel("X pixel",fontsize=25)
        self.ax.set_ylabel("Y pixel",fontsize=25)
        #ax.set_xticks([0,50,100,150,200])
        #ax.set_yticks([0,50,100,150,200])
        self.ax.set_zlabel("Counts",fontsize=25,labelpad=40)
        #ax.set_zticks([0,10000,20000,30000,40000])
            
        for i in self.ax.get_zticklabels():
            i.set_horizontalalignment("left")
        self.ax.view_init(elev=elev,azim=azim)
        #self.fig.tight_layout()

    def plot_image_3D_overview(self,save=True,zlim=(0,45000.),titlestr_prefix="",stretch=None,cut_levels=[90,100,110],**plot_3D_kwargs):
        self.fig = plt.figure(figsize=(20,10))
        gs = gridspec.GridSpec(1, 3)
        ax = plt.subplot(gs[0, 0])
        bx = plt.subplot(gs[0, 1:],projection="3d")

        self.plot(stretch=stretch,ax=ax)
        self.plot_3D(cut_levels=cut_levels,ax=bx,**plot_3D_kwargs)

        bx.set_zlim(*zlim)
        titlestr = titlestr_prefix
        tt = astropy.time.Time(self.header["DATE-OBS"],format="isot").iso
        titlestr += "\nDate: " + tt
        self.fig.tight_layout()
        self.fig.suptitle(titlestr,fontsize=30,x=0.05,horizontalalignment="left")
        self.fig.subplots_adjust(right=0.81,top=0.94)
        if save:
            fp = filepath.FilePath(self.filename)
            self.fig.savefig(fp.basename+".png",dpi=120)
            print("Saved to",fp.basename+".png")
        
    def savefits(self,filename="",verbose=True,overwrite=True):
        if filename=="":
            fp = filepath.FilePath(self.filename)
            fp.add_suffix("_out")
            filename=fp._fullpath
        astropy.io.fits.writeto(filename,data=self.data,header=self.header,overwrite=overwrite)
        if verbose: print("Saved to",filename)

    def close(self):
        self.hdulist.close()

class FitsImageFLI(FitsImage):
    """
    A class that works with the FLI specifically, and takes out the overscan
    """

    def __init__(self,filename):
        super(FitsImageFLI,self).__init__(filename)
        self.data = self.data[:,8:self.data.shape[1]-8]
