#!/usr/bin/python
# -*- coding: utf-8 -*-
# File: fitsimagefli.py
# Created: 2016-03-09 by gks 
"""
Description: 
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class FitsImageFLI(object):
    """
    A helper class when reading fits files from the FLI camera. Depends on pyFits.
    """
    DIRLOC = ''
    
    def __init__(self,filename):
        self.filename = filename
        
        self.hdulist = fits.open(self.filename)
        #print self.hdulist.info()
        self.header = self.hdulist[0].header
        data = self.hdulist[0].data
        self.data = data.astype(float)
        #print self.data.shape
        self.exptime = self.header["EXPTIME"]
        self.dateobs = self.header["DATE-OBS"]
        #self.xpixelsz = self.header["XPIXELSZ"] # Pixel x size in microns
        #self.ypixelsz = self.header["YPIXELSZ"] # Pixel y size in microns
        #self.pltscale = self.header["PLTSCALE"] # "/mm
        #self.pltsizex = self.header["PLTSIZEX"] # plate x dimension in mm
        #self.pltsizey = self.header["PLTSIZEY"] # plate y dimension in mm
    
    #def check_size(self):
    #    """Print out """
    #    print "Platescale: ", self.pltscale, '"/mm'
    #    print "Shape: ", self.data.shape[0], self.data.shape[1]
    #    print "X pixel size: ", self.xpixelsz, "micron"
    #    print "Y pixel size: ", self.ypixelsz, "micron"
    #    print "X-arcsec/pixel:", self.pltscale * self.xpixelsz / 1000.0, '"/pixel'
    #    print "Y-arcsec/pixel:", self.pltscale * self.ypixelsz / 1000.0, '"/pixel'
    #    print "X-arcsec:", self.pltscale * self.xpixelsz / 1000.0 * self.data.shape[0], '"'
    #    print "Y-arcsec:", self.pltscale * self.ypixelsz / 1000.0 * self.data.shape[1], '"'
    #    print "Platesize X: ", self.pltsizex, "mm"
    #    print "Platesize Y: ", self.pltsizey, "mm"
    #    print "\n"
    
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
        
    def cropcenter(self,w,h,points=False):
        """
        Returns an image array around the center of an image array.
        """
        shape = self.data.shape
        center = (shape[0] / 2, shape[1] / 2)
        hl = center[1]-h/2
        hr = center[1]+h/2
        wl = center[0]-w/2
        wr = center[0]+w/2
        self.data = self.data[hl:hr,wl:wr]
        #if points == False:
        #    return self.data_crop
        #else:
        #    return self.data_crop, (wl,wl,wr,wr), (hl,hr,hl,hr)
        
    def cropcord(self,cord,w,h,points=False):
        """
        Returns an image array around the center of an image array.
        
        cord is a tuple
        """
        shape = self.data.shape
        center = cord
        hl = center[1]-h/2
        hr = center[1]+h/2
        wl = center[0]-w/2
        wr = center[0]+w/2
        self.data = self.data[hl:hr,wl:wr]
        
    def plot(self,cmap='gray',figsize=(20,8),vmin=None,vmax=None,scaled=False,colorbar=True,numticks=5):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=1,figsize=figsize)
        if scaled == False:
            #cmap: hot, gray, spectral, hot_r, gray_r, spectral_r
            im = self.axes.imshow(self.data, cmap=cmap,vmin=vmin,vmax=vmax)
            self.axes.set_title("Date:" + self.dateobs + " EXP: " + str(self.exptime))
            if colorbar==True:
                self.fig.colorbar(im)
        else:
            im = self.axes.imshow(self.data/self.exptime, cmap=cmap,vmin=vmin,vmax=vmax)
            self.axes.set_title("SCALED: Date:" + self.dateobs + " EXP: " + str(self.exptime))
            if colorbar==True:
                self.fig.colorbar(im)
        # 5 bins in 
        plt.locator_params(nbins=numticks)
        
    def hist(self,num=1000,figsize=(20,8)):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=1,figsize=figsize)
        self.axes.hist(list(self.data.flat.base),num)
        
    def plot_hist(self,cmap='gray',figsize=(20,8),vmin=None,vmax=None,scaled=False,num=1000,titlepre=""):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=2,figsize=figsize)
        im = self.axes[0].imshow(self.data,cmap=cmap,vmin=vmin,vmax=vmax)
        self.axes[0].set_title(titlepre + "DATE:" + self.dateobs + " EXP: " + str(self.exptime) + " SHAPE: " +str(self.data.shape))
        self.fig.colorbar(im)
        self.axes[1].hist(list(self.data.flat),num,log=True)

    def get_row(self,row):
        """
        Get row data 
        """
        self.data_row = self.data[row,:]
        return self.data_row
        
    def plot_row(self,row,scaled=False,figsize=(20,8)):
        self.fig, self.axes = plt.subplots(nrows=1,ncols=1,figsize=figsize)
        self.data_row = self.data[row,:]
        if scaled == False:
            self.axes.plot(self.data_row)
        else:
            self.axes.plot(self.data_row/np.max(self.data_row)) 
            self.axes.set_ylim(0,1.05)

    def plot_3d(self):
        """

        See here: http://matplotlib.org/examples/mplot3d/contourf3d_demo2.html
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        X = range(0,self.data.shape[0])
        Y = range(0,self.data.shape[1])
        Z = self.data

        xx, yy = np.meshgrid(X,Y)
        surf = ax.plot_surface(xx,yy,Z,cmap=cm.coolwarm,antialiased=True,linewidth=0.1,alpha=0.9,rstride=80,cstride=80)
        #ax.set_xlabel("X")
        #ax.set_ylabel("Y")
        #ax.set_zlabel("Counts")
        #plt.tight_layout()

        # Extra countours
        #cset = ax.contourf(X, Y, Z, zdir='z', offset=30000, cmap=cm.coolwarm)
        #cset = ax.contourf(xx, yy, Z, zdir='x', offset=0, cmap=cm.coolwarm)
        #cset = ax.contourf(xx, yy, Z, zdir='y', offset=max(Y), cmap=cm.coolwarm)

        fig.colorbar(surf,shrink=0.5,aspect=8)
        fig.tight_layout()

        
    def close(self):
        self.hdulist.close()
