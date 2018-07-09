import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import utils
from frdimg_help import add_y_in_input, resample_df_mean, resample_df_mean, get_EE_from_rad_in_pix

# ---------- Graphics ------------
import seaborn as sns; sns.set()
sns.set_context("poster",font_scale=1.2,rc={"font":"helvetica"});
sns.set_style("white"); #sns.set_style("ticks")
cp = sns.color_palette("colorblind") #sns.palplot(current_palette)
from matplotlib import rcParams
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


def plot_final_panel(df,fibername,title='FRD main analaysis',outfolder=""):
    """
    Plot main plot and save it to a folder
    
    INPUT:
        df: pandas dataframe with columns: 
            f_ratio_in
            f_ratio_out
            EE_in_input_cone
        title: title of the plot
        outfolder: folder to save to
        outname: output save name
        
    OUTPUT:
        Creates 

    EXAMPLE:
        plot_final_panel(df_config_f01,outfolder=BASEFOLDER)
    """
    # Resample
    df_mean = resample_df_mean(df)
    resample_array = [2.5,3.33,3.5,3.65,4.,5.,6.,7.]
    resample_array = [2.5,2.73,3.,3.33,3.5,3.65,3.75,4.,4.62,5.,6.,6.67,7.]

    df_res = utils.resample_df(df_mean,resample_array,"f_ratio_in")
    
    # Start plot
    fig, axx = plt.subplots(ncols=3,figsize=(18,7))
    ax, bx, cx = axx[0], axx[1], axx[2]
    
    # FRD plot
    ax.plot(df_mean.f_ratio_in,df_mean.f_ratio_out,color=cp[0],label="EE96",alpha=0.8,lw=1,marker="o",markersize=5)
    ax.set_xlabel("Input F/#")
    ax.set_ylabel("Output F/#")
    xx = np.linspace(1.6,8,100)
    ax.plot(xx,xx,color="black",alpha=0.4,lw=1,ls="--")
    ax.legend(loc="upper left",fontsize=14)
    ax.scatter(df.f_ratio_in,df.f_ratio_out,color=cp[0],alpha=0.8,s=10,marker="1")
    
    # EE Plot
    bx.plot(df_res.f_ratio_in,df_res.EE_in_input_cone,label="",color=cp[0],alpha=0.8,lw=1,marker="o",markersize=5)
    bx.set_xlabel("Input F/#")
    bx.set_ylabel("Output EE within input F/#")
    
    for xx in [ax,bx]:
        xx.minorticks_on()
        xx.grid(lw=0.5,alpha=0.5)

    # Plot table
    plot_table(df_res,round=2,ax=cx)
    fig.tight_layout()
    utils.make_dir(outfolder)
    savename = os.path.join(outfolder,title+"_"+fibername+'.png')
    fig.subplots_adjust(top=0.9)
    fig.suptitle(title+"_"+fibername,fontsize=32,y=0.97)
    fig.savefig(savename,dpi=200)
    print("Saved final plot to {}".format(savename))
    csvout = os.path.join(outfolder,title+"_"+fibername+'.csv')
    df_res.to_csv(csvout,index=False)
    print('Saved final csv to {}'.format(csvout))
    
def plot_table(df,round=2,ax=None,fontsize=18,yt_diff=0.07,yt_start=0.9):
    """
    Plot a simple table
    
    INPUT:
        df: dataframe with keys:
            f_ratio_in,
            df.f_ratio_out, 
            df.EE_in_input_cone 
    
    OUTPUT:
        matplotlib plot
    
    EXAMPLE:
        plot_table(df_int_f01)
    """
    if ax is None:
        fig, ax = plt.subplots(dpi=200)
    
    columns = []
    yt = yt_start
    columns = ['F/#in','F/#out','EE in input cone (%)']
    ax.annotate(columns[0],xy=(0.05,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
    ax.annotate(columns[1],xy=(0.3,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
    ax.annotate(columns[2],xy=(0.7,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
    yt -= yt_diff
    for l,m,r in zip(df.f_ratio_in, df.f_ratio_out, df.EE_in_input_cone):
        ax.annotate('{:0.2f}'.format(l),xy=(0.05,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
        ax.annotate('{:0.2f}'.format(m),xy=(0.3,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
        ax.annotate('{:0.1f}'.format(r*100.),xy=(0.7,yt),xycoords='axes fraction',ha='center',fontsize=fontsize)
        yt -= yt_diff
    ax.axis('off')