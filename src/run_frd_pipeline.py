from __future__ import print_function
import pandas as pd
import numpy as np
import datetime
import matplotlib.dates as mdates
import utils
import phothelp
from frdimg import FRDImg, AnalyzeFRDImages
from frdimg_help import add_y_in_input, add_y_in_xcone, resample_df_mean, resample_df_mean, get_EE_from_rad_in_pix
import os
import frd_plot
import glob
import sys
from astropy.io import fits

# ---------- Graphics ------------
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns; sns.set()
sns.set_context("poster",font_scale=1.2,rc={"font":"helvetica"});
sns.set_style("white"); #sns.set_style("ticks")
cp = sns.color_palette("colorblind") #sns.palplot(current_palette)
rcParams["savefig.dpi"] = 100
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['font.weight'] = "normal"
rcParams["axes.formatter.useoffset"] = False
rcParams['xtick.major.width']=1
rcParams['xtick.major.size']=4
rcParams['xtick.minor.width']=0.5
rcParams['xtick.minor.size']=2
rcParams['xtick.direction'] = "in"
rcParams['ytick.direction'] = "in"
rcParams['ytick.major.width']=1
rcParams['ytick.major.size']=4
rcParams['ytick.minor.width']=0.5
rcParams['ytick.minor.size']=2
# ---------- Graphics ------------
sys.path.append("../src/")
try:
    sys.path.remove('/Users/gks/Dropbox/mypylib')
    sys.path.remove('/Users/gks/Dropbox/mypylib/GIT')
except Exception as e:
    print('Already removed')
   

   
#BASEFOLDER         = "/Users/gks/Dropbox/mypylib/data/NEID_FRD/20180602_HPF_fiberC_200um/"   

def fn_analyze_FRD_data(BASEFOLDER = "C:\\Users\\szk381\\Google Drive\\PSU-file_storage\\NEID\\FRD_data\\20180602_science6_polished_50um\\",
                        FOLDER_CSV_SETUP = None,FOLDER_CSV_SAVE = None, PLOT_FOLDER = None, MASTER_PLOT_FOLDER = None, TITLE = None,
                        MAXRAD_FACTOR      = 0.52, FWZM = 200.,FIBER_NAMES = None, soft_bg_est = True, ADD_F_NUMBER_OUTCONE_VALUE = 3.65):
    print(BASEFOLDER)                   
    if FOLDER_CSV_SETUP == None:
        FOLDER_CSV_SETUP = os.path.join(BASEFOLDER,"ANALYSIS","CSV_SETUP","")
    if FOLDER_CSV_SAVE == None:
        FOLDER_CSV_SAVE = os.path.join(BASEFOLDER,"ANALYSIS","CSV_RESULTS","")    
    if PLOT_FOLDER == None:
        PLOT_FOLDER = os.path.join(BASEFOLDER,"ANALYSIS","PLOTS","")
    if MASTER_PLOT_FOLDER == None:
        MASTER_PLOT_FOLDER = os.path.join(BASEFOLDER,"ANALYSIS","MASTER_PLOTS","")
    if TITLE == None:
        TITLE = BASEFOLDER.split(os.sep)[-2]    
        
    if FIBER_NAMES == None:
        contents = np.array(os.listdir(BASEFOLDER))
        fiber_dir = np.where(np.isin(contents,['ANALYSIS', 'background', 'bias',  'psm_images','desktop.ini']) == False)[0]
        FIBER_NAMES = contents[fiber_dir]
        print(FIBER_NAMES)
        
    try:
        bias = fits.open(os.path.join(BASEFOLDER,'bias','median_bias.fits'))[0].data
        print('Found bias frame')
    except:
        bias = None
        print('Bias frame not found')
        

    print('\n\nRUNNING FRD PIPELINE FOR {}\n\n'.format(TITLE))
        
    # Run for all fibers
    for FIBER_NAME in FIBER_NAMES: 
    
        # Files
        setup_csv_files  = sorted(glob.glob(FOLDER_CSV_SETUP+"*.csv"))
        fitsfiles_f01    = glob.glob(BASEFOLDER+FIBER_NAME+"/*.fits")
        
        if len(fitsfiles_f01) <= 2:
            next
    
        plt.switch_backend("agg")
    
        ##########################
        # MAIN analysis step
        AFRDImg = AnalyzeFRDImages(fitsfiles = sorted(fitsfiles_f01),
                                biasframe = bias,
                                plot_suffix = FIBER_NAME,
                                plot_folder = PLOT_FOLDER,
                                setup_csv_file = setup_csv_files[0],
                                FOLDER_CSV_SAVE = FOLDER_CSV_SAVE,
                                motorized = True,
                                MAXRAD_FACTOR = MAXRAD_FACTOR,
                                fwzm_z = FWZM,
                                use_azimuthal_averaging = True,
                                soft_bg_est = soft_bg_est)
        AFRDImg.analyze_all_frames()
        ##########################
    
    
        ##########################
        # Find  .csv files
        list_of_df_config = sorted(glob.glob(FOLDER_CSV_SAVE+"*.csv"))
    
        df_config_f01 = pd.read_csv(list_of_df_config[0])
    
        print("Calculating input/output cone values")
        # Add column for input cone 
        df_config_f01 = add_y_in_input(df_config_f01)
        df_config_f01 = AFRDImg._get_EE_in_input_cone(df_config_f01,suffix=FIBER_NAME)
        
        print("Calculating 3.65 cone values")
        # Add column for 3.65 ( ADD_F_NUMBER_OUTCONE_VALUE = 3.65 ) input cone
        df_config_f01 = add_y_in_xcone(df_config_f01,f_ratio_of_output_cone=ADD_F_NUMBER_OUTCONE_VALUE)
        df_config_f01 = AFRDImg._get_EE_in_xcone(df_config_f01,suffix=FIBER_NAME,f_ratio_of_output_cone=ADD_F_NUMBER_OUTCONE_VALUE)
        ##########################
       
        # Plot main plot
        frd_plot.plot_final_panel(df_config_f01,fibername=FIBER_NAME,
                                  title=TITLE,outfolder=MASTER_PLOT_FOLDER,
                                  f_ratio_of_output_cone=ADD_F_NUMBER_OUTCONE_VALUE)
        
        print('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    return df_config_f01
        
    
if __name__=="__main__":
    home = r'C:\Users\szk381\Google Drive\PSU-file_storage\NEID\FRD_data'
    #home = '/Users/szk381/Google Drive/PSU-file_storage/NEID/FRD_data'	
    folder_name = "20180718_science1_puck_polished_final_50um"
    folder_name = '20180720_neid_test_fiberhead'
    folder_name = '20180725_science6_polished_200um_input'
    folder_name = '20180803_science4_retest'
    #folder_name = '20180723_science1_puck_polished_final_brass_polished_final_50um'
    folder_name = '20180806_science4_thorlabs_635'
    fn_analyze_FRD_data(BASEFOLDER = os.path.join(home,folder_name,""),FIBER_NAMES = ['HR2_Intsph - Copy'])
        		
