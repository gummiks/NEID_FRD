from __future__ import print_function
import pandas as pd
import numpy as np
import utils

def get_EE_from_rad_in_pix(df,radius_in_pix):
    print("Now getting EE from rad in pix",np.interp(radius_in_pix,df.radii.values,df.EE.values))
    #print("Now getting EE from rad in pix",df[df["radii"] > radius_in_pix]["EE"].values[0])
    return np.interp(radius_in_pix,df.radii.values,df.EE.values)
    #return df[df["radii"] > radius_in_pix]["EE"].values[0]

    #r_at_EE = df[df["EE"] > get_rad_at_EE]["radii"].values[0]
    #print(r_at_EE)
    # better to interpolate

def add_y_in_input(df_config):
    df_config["y_in_input"] = (df_config["f_ratio_out"]/df_config["f_ratio_in"])*df_config["y_out_fiber_dist"]
    return df_config

def add_y_in_xcone(df_config,f_ratio_of_output_cone=3.65):
    """
    Similar to add_y_in_input, but can do for arbitrary cone
    """
    df_config["y_in_"+str(f_ratio_of_output_cone)+"_cone"] = (df_config["f_ratio_out"]/f_ratio_of_output_cone)*df_config["y_out_fiber_dist"]
    return df_config

def resample_df_mean(df,f_ratio_of_output_cone=3.65):
    eelabel = "EE_in_"+str(f_ratio_of_output_cone)+"_cone"
    chunks = utils.chunker(df,4)
    dff = pd.DataFrame()
    for chunk in chunks:
        dff = pd.concat([dff,pd.DataFrame(chunk[["f_ratio_in","f_ratio_out","EE_in_input_cone",eelabel]].median(skipna=True)).T],
                        ignore_index=True)
    return dff

