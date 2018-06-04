import utils
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import butter
import photutils
import glob
import fitsimg
#import fitsimagefli
import phothelp
from astropy.modeling import models, fitting
from frdimg_help import add_y_in_input, resample_df_mean, resample_df_mean, get_EE_from_rad_in_pix
import filepath

class FRDImg(object):
    """
    A class to analyze FRDImages
    """
    
    def __init__(self,filename):
        self.filename   = filename
        self.basename   = self.filename.split(os.sep)[-1]
        self.fimg       = fitsimg.FitsImageFLI(self.filename)
        print("Trimming problematic first row !")
        self.fimg.data = self.fimg.data[1:1027,0:1056]
        
    def run(self,
            plot_suffix,
            plot_folder,
            get_rad_at_EE = 0.96,
            MAXRAD_FACTOR = 0.54,
            fwzm_z=20,
            use_azimuthal_averaging=True,
            force_max_radii_to_390=False):
        """
        Run the analysis
        
        INPUT:
            plot_suffix = "f01"
            plot_folder = "results/"
            get_rad_at_EE = 0.96
            FACTOR = 0.54
            fwzm_z=20

        OUTPUT:
            Saves files to folders  
        """
        utils.make_dir(plot_folder)
        fig, ax = plt.subplots(ncols=3,nrows=2,figsize=(20,15))
        
        ax_img         = ax.flat[0]
        ax_img_bkg     = ax.flat[1]
        ax_img_bkg_sub = ax.flat[2]
        ax_cut         = ax.flat[3]
        ax_ee          = ax.flat[4]
        ax_img_bkg_sub2= ax.flat[5]

        # -------------
        # Plot #1: Normal image
        x,y = self.fimg.get_centroid(plot_lines=True,ax=ax_img)
        self.fimg.ax.set_title("Bias + Dark Subtracted (Hist stretch)")
        #cbar_ax = fig.add_axes([0.38, 0.59, 0.015, 0.3],label="Counts")
        #fig.colorbar(self.fimg.im, cax=cbar_ax)
        
        # -------------
        # Plot #2: Background estimation
        x,y = self.fimg.get_centroid()
        print("Centroids",x,y)
        CBS = phothelp.CircularBackgroundSubractor(self.fimg.data,x=x,y=y,r=400)
        data_back_sub = CBS.subtract_background(plot_background=True,ax=ax_img_bkg)
        self.subtracted_value = CBS.subtracted_value
        
        # -------------
        # Plot #3: Background subtracted image
        self.fimg.data = data_back_sub
        self.fimg.plot(ax=ax_img_bkg_sub,colorbar=False,title="Bias + Dark + Background subtracted (linear stretch)",stretch="linear")
        # add colorbar
        #cbar_ax = fig.add_axes([0.38, 0.09, 0.015, 0.3],label="Counts")
        #fig.colorbar(self.fimg.im, cax=cbar_ax)
        # Get centroid again after taking out the background
        x,y = self.fimg.get_centroid()        
        
        
        # -------------
        # Plot #4: Cuts
        if use_azimuthal_averaging == True:
            #self.max_radii_for_EE = 
            print("Using azimuthal averaging")
            _r,_mean,_r_hwzm = self.fimg.get_radial_profile(rmax=500,plot=True,z=fwzm_z,return_hwzm=True,ax=ax_cut,xcen=x,ycen=y,subtract_min=True) # 500 ~2/1024
            self.max_radii_for_EE = _r_hwzm*2.*MAXRAD_FACTOR
            if force_max_radii_to_390 == True:
                print("FORCING MAX RADII TO 390")
                self.max_radii_for_EE = 490.
        else:
            print("Not using azimuthal averaging, Using mean(X+Y) cutouts")
            cut_x, self.fwhm_x = self.fimg.get_centroid_line_cut(line="X",return_FWHM=True,plot=True,ax=ax_cut,use_butter_filter=True,butter_cutoff_freq=0.03,fwzm_z=fwzm_z)
            cut_y, self.fwhm_y = self.fimg.get_centroid_line_cut(line="Y",return_FWHM=True,plot=True,ax=ax_cut,use_butter_filter=True,butter_cutoff_freq=0.03,fwzm_z=fwzm_z)        
            #return cut_x, cut_y
            ax_cut.plot(butter.low_pass_butter(cut_x,cutoff_freq=0.03))
            ax_cut.plot(butter.low_pass_butter(cut_y,cutoff_freq=0.03))
            self.max_radii_for_EE = ((self.fwhm_x + self.fwhm_y)/2.)*MAXRAD_FACTOR
        
        ax_cut.set_ylim(0,64000.)

        # -------------
        # Plot #4: Encircled energy plot
        

        radii = np.arange(1,self.max_radii_for_EE)
        df, self.r_ee = phothelp.get_encircled_energy_and_rad_at_EE(self.fimg.data,
                                                                    self.fimg.xcenter,
                                                                    self.fimg.ycenter,
                                                                    radii,
                                                                    plot=True,
                                                                    ax=ax_ee,
                                                                    get_rad_at_EE=get_rad_at_EE)
        ax_ee.set_xlim(0,500)
        ax_ee.set_xlabel("Radii (pixels)")
        
        # Overplot aperture in figure $ with EE
        aper = photutils.CircularAperture((x,y),r=self.r_ee)
        aper.plot(ax=ax_img_bkg_sub,color="green")
        ax_img_bkg_sub.annotate("EE"+str(100*get_rad_at_EE)+"%",xy=(x,y+self.r_ee),color="white",fontsize=8)
        # overplot max aperture 
        aper = photutils.CircularAperture((x,y),r=self.max_radii_for_EE)
        aper.plot(ax=ax_img_bkg_sub,color="red")
        ax_img_bkg_sub.annotate("EE100% @ r = "+str(MAXRAD_FACTOR*2)+" HWHM",xy=(x,y+self.max_radii_for_EE),color="white",fontsize=8)
        
        # Hist plot
        self.fimg.plot(ax=ax_img_bkg_sub2,colorbar=False,title="Bias + Dark + Background subtracted (Hist stretch)",stretch="hist")
        aper = photutils.CircularAperture((x,y),r=self.r_ee)
        aper.plot(ax=ax_img_bkg_sub2,color="green")
        ax_img_bkg_sub2.annotate("EE"+str(100*get_rad_at_EE)+"%",xy=(x,y+self.r_ee),color="white",fontsize=8)
        # overplot max aperture 
        aper = photutils.CircularAperture((x,y),r=self.max_radii_for_EE)
        aper.plot(ax=ax_img_bkg_sub2,color="red")
        ax_img_bkg_sub2.annotate("EE100% @ r = "+str(MAXRAD_FACTOR*2)+" HWHM",xy=(x,y+self.max_radii_for_EE),color="white",fontsize=8)

        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        fig.suptitle("FRD Overview Plot, file: "+self.basename,fontsize=20)
        fig.savefig(plot_folder + self.basename+"_"+plot_suffix+".png")
        print("Saved file:",plot_folder + self.basename+"_"+plot_suffix+".png")
        df.to_csv(plot_folder + self.basename+"_"+plot_suffix+".csv")
        print("Saved file:",plot_folder + self.basename+"_"+plot_suffix+".csv")
        



        
    
class AnalyzeFRDImages(object):
    """
    Analyze a set of FRDImages.
    
    NOTES:
    See *analyze_all_frames()*
    """
    # CONSTANTS
    PIXEL_SCALE = 0.013 # mm
    F_LENGTH_IN = 20.0  # mm
    
    def __init__(self,
                 fitsfiles,
                 plot_suffix,
                 plot_folder,
                 setup_csv_file,
                 FOLDER_CSV_SAVE,
                 motorized=True,
                 get_rad_at_EE=0.96,
                 MAXRAD_FACTOR=0.54,
                 fwzm_z=20,
                 force_max_radii_to_390=False,
                 use_azimuthal_averaging=False):
        self.fitsfiles        = fitsfiles
        self.plot_suffix      = plot_suffix
        self.plot_folder      = plot_folder
        self.setup_csv_file   = setup_csv_file
        self.FOLDER_CSV_SAVE  = FOLDER_CSV_SAVE
        self.MOTORIZED        = True # using motorized iris or not ?
        self.f_number_names   = pd.Series([filepath.FilePath(filename).basename.split("_")[2] for filename in self.fitsfiles]).unique()
        self.NUM_F_NUMBERS    = len(self.f_number_names) # Number of different F-numbers measured
        self.get_rad_at_EE    = get_rad_at_EE # What EE to get
        self.MAXRAD_FACTOR    = MAXRAD_FACTOR # The factor to multiply the average of FWHM_X and FWHM_Y to set as the maximum radius
        self.BACKUP_FILE_NAME = self.FOLDER_CSV_SAVE + self.plot_suffix + "_radii_fwhm_subtracted.npz"
        self.fwzm_z           = fwzm_z
        self.force_max_radii_to_390=force_max_radii_to_390
        self.use_azimuthal_averaging = use_azimuthal_averaging
        utils.make_dir(self.FOLDER_CSV_SAVE)
        
    def _get_df(self,df_config,frat,startnum=0):
        """
        Get data frames, and finds all of the records for a given F input, and 
        calculates the input and output F-ratios.
        
        INPUT: 
            df_config - dataframe with the measurements
            frat      - fiber input measurement set
            
        OUTPUT:
        
        NOTES:
            When using the motorized setup, then need to use the line model to convert to mm.
            
        """
        # Filenames are f_psm_01_d_01.fits
        # Select only measurements for a given input F-number
        df = df_config[df_config.filename.isin(["f_psm_"+frat+"_d_"+str(i).zfill(2)+".fits" for i in range(1,12)])]
        df = df.reset_index(drop=True)
        df = df[startnum:]
        
        #if method=="FIRST":
        df["x_out_fiber_dist_delta"] = df["x_out_fiber_dist"].values[0] - df["x_out_fiber_dist"].values
        df["y_out_fiber_dist_delta"] = df["y_out_fiber_dist"].values - df["y_out_fiber_dist"].values[0]
        df["f_ratio_out"] = df["x_out_fiber_dist_delta"].values/(2.*df["y_out_fiber_dist_delta"].values)
        #elif method=="SEQUENTIAL"
        
        if self.MOTORIZED==True:
            # motorized setup: convert from percent to mm
            df["iris"] = self._percent2mm(df["iris_diam_perc"])
        df["f_ratio_in"]  = self.F_LENGTH_IN/df["iris"]
        return df
    
    def load_backup():
        backup_file = np.load(self.BACKUP_FILE_NAME)
        a, b, c = backup_file.items()
        radii_at_ee, fwhm, subtracted_values = a[1], b[1], c[1]
        return radii_at_ee, fwhm, subtracted_values

    def _calc_f_out(self,radii_at_EE,subtracted_values,save=True):
        """
        Calculate the F output F number for a given EE.
        
        INPUT:
            radii_at_EE - the radius (in pix) used to calculate the 
        
        setup_csv_file = "/Users/gks/Dropbox/mypylib/notebooks/HPF_FRD/frd_files_final_run_f01.csv"
        """
        filepath = utils.FilePath(self.setup_csv_file)
        filepath.add_prefix("df_config_")

        df_config = pd.read_csv(self.setup_csv_file)[0:len(self.fitsfiles)]
        df_config["radii_at_EE"] = radii_at_EE
        df_config["sub_values"]  = subtracted_values
        df_config["y_out_fiber_dist"] = df_config["radii_at_EE"]*self.PIXEL_SCALE
        df_config = pd.concat([self._get_df(df_config,frat=str(i).zfill(2)) for i in range(1,self.NUM_F_NUMBERS+1)],ignore_index=True)
        if save:
            df_config.to_csv(self.FOLDER_CSV_SAVE+filepath.basename)
            print("Saved to:",self.FOLDER_CSV_SAVE+filepath.basename)
        return df_config
    
    def _percent2mm(self,percent):
        """
        Converts between iris size in percent to mm
        """
        m_line = models.Linear1D(slope=0.333201810375,intercept=2.02718927701)
        print("LINE",m_line(percent))
        return m_line(percent)

    def analyze_all_frames(self):
        """
        Analyze all the frames, and saves as a .csv file

        OUTPUT:
        """
        plt.switch_backend("agg")
        self.radii_at_EE       = np.zeros(len(self.fitsfiles))
        self.fwhm              = np.zeros(len(self.fitsfiles))
        self.subtracted_values = np.zeros(len(self.fitsfiles))
        for i,fitsfile in enumerate(self.fitsfiles):
            frd = FRDImg(fitsfile)
            print("Analyzing file #",i,fitsfile)
            frd.run(get_rad_at_EE=self.get_rad_at_EE,
                    MAXRAD_FACTOR=self.MAXRAD_FACTOR,
                    plot_suffix=self.plot_suffix,
                    plot_folder=self.plot_folder,
                    fwzm_z=self.fwzm_z,
                    force_max_radii_to_390=self.force_max_radii_to_390,
                    use_azimuthal_averaging=self.use_azimuthal_averaging);
            print("Radii at EE:",frd.r_ee)
            self.radii_at_EE[i] = frd.r_ee
            self.fwhm[i] = frd.max_radii_for_EE
            self.subtracted_values[i] = frd.subtracted_value
            np.savez(self.BACKUP_FILE_NAME,self.radii_at_EE,self.fwhm,self.subtracted_values)
            print("Saved backup to",self.BACKUP_FILE_NAME)
            print("FWHM",frd.max_radii_for_EE)
            print("")
        self.df_config = self._calc_f_out(self.radii_at_EE,self.subtracted_values,save=True)
        return self.df_config
    
    def _get_EE_in_input_cone(self,df_config,suffix):
        """
        Get encircled energy in the given input cone.
        """
        df_config["y_for_EE_in_input_cone"] = np.zeros(len(df_config))
        df_config["EE_in_input_cone"] = np.zeros(len(df_config))

        for i, filename in enumerate(df_config["filename"].values):
            print(i)
            print(filename)
            if (np.isfinite(df_config["f_ratio_out"].values[i])):
                filepath = utils.FilePath(filename)
                filepath.basename 
                #df_ee = pd.read_csv(self.plot_folder + filepath.barename + "_out.fits_" + suffix + ".csv")
                df_ee = pd.read_csv(self.plot_folder + filepath.barename + ".fits_" + suffix + ".csv")

                
                df_config["y_for_EE_in_input_cone"].ix[i] = df_config["y_in_input"].values[i]/self.PIXEL_SCALE
                print("Rad_in_pix",df_config["y_for_EE_in_input_cone"].ix[i])
                df_config["EE_in_input_cone"].ix[i]       = get_EE_from_rad_in_pix(df_ee,df_config["y_for_EE_in_input_cone"].ix[i])
            else:
                df_config["y_for_EE_in_input_cone"].ix[i] = np.nan
                df_config["EE_in_input_cone"].ix[i]       = np.nan
                print("Skip!")
        return df_config