from __future__ import print_function
import photutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gklib as gk


class CircularBackgroundSubractor(object):
    """
    A class to calculate and subtract the background in a FITS image after masking out a circle of radius r at
    position x and y.

    NOTE:
    For the HPF FRD measurements, have r ~330
    """
    
    def __init__(self,data,x,y,r,box_size=(205,205)):
        self.data = data
        self.x    = x
        self.y    = y
        self.r    = r
        self.box_size = box_size
        
    def subtract_background(self,subtract_min_value=True,plot_background=False,ax=None):
        """
        A function to subtract the background for the FRD tests

        INPUT:
        subtract_min_value - subtracts the .min value (no negative numbers)
        plot_background    - if True, plots the estimated background with the meshes it used.
        """
        self.aper = photutils.CircularAperture(positions=(self.x,self.y),r=self.r)

        # get the mask from the aperture
        # hack
        mask = self.aper.to_mask()[0].to_image(self.data.shape)

        # use sigma clipping
        sigma_clip = photutils.SigmaClip(sigma=3., iters=10)

        bkg_estimator = photutils.MedianBackground()
        self.bkg = photutils.Background2D(self.data, self.box_size, filter_size=(3, 3),
                            sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,mask=mask,edge_method="crop")

        self.data_background = self.data - self.bkg.background

        if subtract_min_value:
            self.subtracted_value      = self.data_background.min()
            print("Subtracted min value:",self.subtracted_value)
            self.data_background       -= self.subtracted_value
        
        if plot_background:
            if ax==None:
                self.fig, self.ax = plt.subplots()
            else:
                self.ax = ax
            im = self.ax.imshow(self.bkg.background, origin='lower', cmap='Greys_r')
            self.ax.set_xlim(0,self.bkg.background.shape[1])
            self.ax.set_ylim(0,self.bkg.background.shape[0])
            self.ax.set_title("Estimated Background",y=1.02)
            self.ax.set_xlabel("X pixels")
            self.ax.set_ylabel("Y pixels")
            if ax==None:
                # Only do this if ax is not supplied
                # Don't know how to deal with this without passing the figure explicitly
                self.fig.colorbar(im)
            self.bkg.plot_meshes(outlines=True, color='#1f77b4',ax=self.ax)
            
        return self.data_background


def howell_center(postage_stamp):
    """
    Howell centroiding, from Howell's Handbook of CCD astronomy

    INPUT:
     postage_stamp - A 2d numpy array to do the centroiding

    OUTPUT:
     x and y center of the numpy array

    NOTES:
    Many thanks to Thomas Beatty and the MINERVAphot.py pipeline for this method
    see here: https://github.com/TGBeatty/MINERVAphot/blob/master/MINERVAphot.py
    """
    xpixels = np.arange(postage_stamp.shape[1])
    ypixels = np.arange(postage_stamp.shape[0])
    I = np.sum(postage_stamp, axis=0)
    J = np.sum(postage_stamp, axis=1)

    Isub = I-np.sum(I)/I.size
    Isub[Isub<0] = 0
    Jsub = J-np.sum(J)/J.size
    Jsub[Jsub<0] = 0
    xc = np.sum(Isub*xpixels)/np.sum(Isub)
    yc = np.sum(Jsub*ypixels)/np.sum(Jsub)
    return xc, yc

def get_encircled_energy_and_rad_at_EE(data,x,y,radii,get_rad_at_EE=0.9,plot=False,ax=None):
    """
    A function to calculate the encircled energy at a given position, summing up the flux in apertures of size *radii*,
    normalizing the EE to the last flux value.
    
    INPUT:
    data - a two dimensional np.array
    x - x centroid
    y - y centroid
    radii - an array of radii to calculate the EE
    get_rad_at_EE - the EE of which to calculate the radius value
    plot - plot a EE vs radii plot
    
    OUTPUT:
    df
        a dataframe with two columns:
        - radii
        - EE
    rad_at_EE
        - the radius when EE is a given input value

    NOTE:
    Assumes that the data is background subtracted

    EXAMPLE:
    radii = np.arange(1,450)
    df_ee, r_at_EE90 = phothelp.get_encircled_energy(fimg.data,fimg.xcenter,fimg.ycenter,radii,plot=True)
    """
    apertures = [photutils.CircularAperture((x,y), r=r) for r in radii]
    phot_table = photutils.aperture_photometry(data, apertures)
    df = phot_table[phot_table.colnames[3:]].to_pandas()
    EE = df.loc[0]/df.loc[0][-1]
    
    df = pd.DataFrame(zip(radii,EE),columns=["radii","EE"])

    r_at_EE = df[df["EE"] > get_rad_at_EE]["radii"].values[0]

    if plot:
        if ax==None:
            fig, ax = plt.subplots()
        ax.plot(radii,EE.values)
        ax.set_xlabel("Radii")
        ax.set_ylabel("Encircled Energy")
        ax.set_title("EE"+gk.num2str(get_rad_at_EE*100,1)+"% at r="+str(r_at_EE))
        
        ax.vlines(r_at_EE,0,get_rad_at_EE,color="orange",linestyle="--")
        ax.hlines(get_rad_at_EE,0,r_at_EE,color="orange",linestyle="--")
        ax.minorticks_on()
        
    return df, r_at_EE
