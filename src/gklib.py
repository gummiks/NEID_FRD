import numpy as np
import matplotlib.dates as mdates
import subprocess
import matplotlib
import matplotlib.pyplot as plt

def num2str(num,precision=3,delim=" "):
    """A function that returns a formatted float number"""
    if type(num)==float or type(num)==np.float64:
        return "%0.*f" % (precision,num)
    else:
        return delim.join("%0.*f" % (precision,i) for i in num)

def num2latex(f):
    float_str = "{0:.4g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        number = r"{0} \times 10^{{{1}}}".format(base, int(exponent))
        with_dollar = "$"+number+"$"
        return with_dollar
    else:
        return float_str

#------------------------------------------------------------------------
#Plotting functions
#------------------------------------------------------------------------
def conf_matplot():
    """Initialize matplotlib things: font size etc"""
    #------------Figure Layout--------------
    #---------------------FONT and other graphics------------------------
    font = {'family'         : 'serif',
            'size'	         : 14}
    matplotlib.rc('font',**font)
    matplotlib.rc('grid',linewidth=1)
    matplotlib.rc('xtick.major',width=2)
    matplotlib.rc('xtick.major',size=7)
    matplotlib.rc('xtick.minor',width=2)
    matplotlib.rc('xtick.minor',size=4)
    matplotlib.rc('ytick.major',width=2)
    matplotlib.rc('ytick.major',size=7)
    matplotlib.rc('ytick.minor',width=2)
    matplotlib.rc('ytick.minor',size=4)
    #-------------------------------------------------
    pass

def gplot(x,y,col="black",lw=2,xlab="",ylab="",title="",filename="",leg="",legloc="",xlim=0,ylim=0,logx=False,logy=False,show=True,invertx=False):
    """
    DESCRIPTION:
         To plot y as a function of x and save to file *filename*.
    
    CALLING SEQUENCE:
           plot_mesa_fig(x,y,"x","y","myfile.pdf",mylegend="a",\
                   xlim=(0,1),ylim=(0,1))
    
    INPUT PARAMETERS:
           x        - array with n elements       
           y        - array with n elements
           col      - color
           lw       - linewidth
           xlab     - r"$x$"
           ylab     - r"$y$"
           title    - The plot title
           filename - Filename e.g. "myfile.pdf"
           leg      - Legend
           legloc   - E.g. "upper right", "center right"
           xlim     - E.g. (0,1)
           ylim     - E.g. (0,1)
           logx     - True/False
           logy     - True/False
           show     - Show window
           invertx  - Invert the x axis
    
    OPTIONAL INPUT:
           See above
    
    OUTPUT PARAMETERS:
           None
    
    NOTES:
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date, Feb   17
                                             March 24, 2015
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,title=title)
    adjustprops = dict(left=0.19,bottom=0.15,right=0.92,top=0.9,wspace=0.,hspace=0.2)
    fig.subplots_adjust(**adjustprops)
    if xlab!="":
        ax.set_xlabel(xlab,size="x-large")
    if ylab!="":
        ax.set_ylabel(ylab,size="x-large")
    ax.minorticks_on()
    ax.grid()
    if xlim != 0:
        ax.set_xlim(xlim)
    if ylim != 0:
        ax.set_ylim(ylim)
    #PLOT
    ax.plot(x,y,color=col,linewidth=lw,linestyle="-",label=leg)
    if legloc!="":
        ax.legend(loc=legloc,prop={'size':16})
    else:
        ax.legend(loc="upper center",prop={'size':16})
    if logx == True:
        ax.set_xscale("log")
    if logy == True:
        ax.set_yscale("log")
    if invertx==True:
        plt.gca().invert_xaxis()
    if filename !="":
        fig.savefig(filename)
        print("Plot saved to file, "+filename)
    if show==False:
        plt.close()
    if show==True:
        fig.show()

def gplot2(x,y1,y2,col="black",lw=2,xlab="",ylab="",title="",filename="",leg1="",leg2="",legloc="",xlim=0,ylim=0,logx=False,logy=False,show=True):
    """
    DESCRIPTION:
         To plot y1 and y2 as a function of x and save to file *filename*.
    
    CALLING SEQUENCE:
    
    INPUT PARAMETERS:
           x        - array with n elements       
           y1       - array with n elements
           y2       - array with n elements
           col      - color
           lw       - linewidth
           xlab     - r"$x$"
           ylab     - r"$y$"
           title    - Title of the plot
           filename - Filename e.g. "myfile.pdf"
           leg1     - Legend1
           leg2     - Legend2
           legloc   - E.g. "upper right", "center right"
           xlim     - E.g. (0,1)
           ylim     - E.g. (0,1)
           logx     - True/False
           logy     - True/False
           show     - Show window
    
    OPTIONAL INPUT:
           See above
    
    OUTPUT PARAMETERS:
           None
    
    NOTES:
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date, Feb   17
                                             March 24, 2015
                                             March 31, 2015
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,title=title)
    adjustprops = dict(left=0.19,bottom=0.15,right=0.92,top=0.9,wspace=0.,hspace=0.2)
    fig.subplots_adjust(**adjustprops)
    if xlab!="":
        ax.set_xlabel(xlab,size="x-large")
    if ylab!="":
        ax.set_ylabel(ylab,size="x-large")
    ax.minorticks_on()
    ax.grid()
    if xlim != 0:
        ax.set_xlim(xlim)
    if ylim != 0:
        ax.set_ylim(ylim)
    #PLOT
    ax.plot(x,y1,color=col,linewidth=lw,linestyle="-",label=leg1,alpha=0.5)
    ax.plot(x,y2,color=col,linewidth=lw,linestyle="--",label=leg2)
    if legloc!="":
        ax.legend(loc=legloc,prop={'size':16})
    else:
        ax.legend(loc="upper center",prop={'size':16})
    if logx == True:
        ax.set_xscale("log")
    if logy == True:
        ax.set_yscale("log")
    if filename !="":
        fig.savefig(filename)
        print("Plot saved to file, "+filename)
    if show==False:
        plt.close()
    if show==True:
        fig.show()

def gplot3(x,y1,y2,y3,col="black",lw=2,xlab="",ylab="",title="",filename="",leg1="",leg2="",leg3="",legloc="",xlim=0,ylim=0,logx=False,logy=False,show=True):
    """
    DESCRIPTION:
         To plot y1 and y2 and y3  as a function of x and save to file *filename*.
    
    CALLING SEQUENCE:
    
    INPUT PARAMETERS:
           x        - array with n elements       
           y1       - array with n elements
           y2       - array with n elements
           y3       - array with n elements
           col      - color
           lw       - linewidth
           xlab     - r"$x$"
           ylab     - r"$y$"
           title    - Title of plot
           filename - Filename e.g. "myfile.pdf"
           leg1     - Legend1
           leg2     - Legend2
           leg2     - Legend3
           legloc   - E.g. "upper right", "center right"
           xlim     - E.g. (0,1)
           ylim     - E.g. (0,1)
           logx     - True/False
           logy     - True/False
           show     - Show window
    
    OPTIONAL INPUT:
           See above
    
    OUTPUT PARAMETERS:
           None
    
    NOTES:
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date, Feb   17
                                             March 24, 2015
                                             March 31, 2015
                                             April 11, 2015
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,title=title)
    adjustprops = dict(left=0.19,bottom=0.15,right=0.92,top=0.9,wspace=0.,hspace=0.2)
    fig.subplots_adjust(**adjustprops)
    if xlab!="":
        ax.set_xlabel(xlab,size="x-large")
    if ylab!="":
        ax.set_ylabel(ylab,size="x-large")
    ax.minorticks_on()
    ax.grid()
    if xlim != 0:
        ax.set_xlim(xlim)
    if ylim != 0:
        ax.set_ylim(ylim)
    #PLOT
    ax.plot(x,y1,color=col,linewidth=lw,linestyle="-", label=leg1,alpha=0.5)
    ax.plot(x,y2,color=col,linewidth=lw,linestyle="--",label=leg2)
    ax.plot(x,y3,color=col,linewidth=lw,linestyle=":", label=leg3)
    if legloc!="":
        ax.legend(loc=legloc,prop={'size':16})
    else:
        ax.legend(loc="upper center",prop={'size':16})
    if logx == True:
        ax.set_xscale("log")
    if logy == True:
        ax.set_yscale("log")
    if filename !="":
        fig.savefig(filename)
        print("Plot saved to file, "+filename)
    if show==False:
        plt.close()
    if show==True:
        fig.show()

def hrd(x,y,title="",filename=""):
    """
    DESCRIPTION:
           A function to plot the HR diagram
    
    CALLING SEQUENCE:
            hrd(x,y)  
    
    INPUT PARAMETERS:
           x - The effective temperature in K
           y - The luminosity in units of log(L/L_sun)
    
    OPTIONAL INPUT:
    
    OUTPUT PARAMETERS:
                 
    NOTES:
           Dependencies: 
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date April 25, 2015
    """
    gplot(x,y,title=title,xlab=r"$T_{\mathrm{eff}}$ [K]",ylab=r"Luminosity $\log( L / L_{\odot} )$",invertx=True,filename=filename)

def plot31(x,y1,y2,y3,xlab="",y1lab="",y2lab="",y3lab="",title="",filename="",show=True):
    """
    DESCRIPTION:
           A plot that plots 3 figure panel, sharing the x axis
    
    CALLING SEQUENCE:
           plot31(x,y1,y2,y3,title="",filename="",show=False)  
    
    INPUT PARAMETERS:
           x  - the x array
           y1 - y1 array
           y2 - y2 array
           y3 - y3 array
           xlab - the x label
           y1lab - the y1 label
           y2lab - the y2 label
           y3lab - the y3 label
    
    OPTIONAL INPUT:
           title=""     - Title of the plot
           filename=""  - Filename to save the figure, with extension
           show=False   - Show plot if True, not if False
    
    OUTPUT PARAMETERS:
           res      
    
    NOTES:
           Dependencies: 
           keywords:
           res = plot31(.. , keyw='')
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date, April 26, 2015
    """

    figprops = dict(figsize=(8.,8./2.118), dpi=256)
    adjustprops = dict(left=0.13,bottom=0.12,right=0.97,top=0.83,wspace=0.,hspace=0.2)

    #fig = plt.figure(**figprops)
    fig = plt.figure()
    fig.subplots_adjust(**adjustprops)

    ax = fig.add_subplot(3,1,1,title=title)
    bx = fig.add_subplot(3,1,2)
    cx = fig.add_subplot(3,1,3)
    #bx = fig.add_subplot(3,1,2,sharex=ax,sharey=ax)
    #cx = fig.add_subplot(3,1,3,sharex=ax,sharey=ax)

    ax.set_ylabel(xlab,size="large")
    bx.set_ylabel(y1lab,size="large")
    cx.set_ylabel(y2lab,size="large")
    cx.set_xlabel(y3lab,size="large")

    ax.grid()
    bx.grid()
    cx.grid()
    ax.minorticks_on()
    bx.minorticks_on()
    cx.minorticks_on()

    plt.setp(ax.get_xticklabels(),visible=False)
    plt.setp(bx.get_xticklabels(),visible=False)

    #PLOT
    ax.plot(x,y1,color="black",linewidth=3,linestyle="-",label="")
    bx.plot(x,y2,color="black",linewidth=3,linestyle="-",label="")
    cx.plot(x,y3,color="black",linewidth=3,linestyle="-",label="")

    if filename!="":
       fig.savefig(filename)
    if show==True:
       plt.show()

def gscatter(x,y,col="black",s=100,xlab="",ylab="",title="",filename="",leg="",legloc="",xlim=0,ylim=0,logx=False,logy=False,show=True,invertx=False):
    """
    DESCRIPTION:
         To make a scatterplot of y as a function of x and save to file *filename*.
    
    CALLING SEQUENCE:
           gk.gscatter(x,y,"x","y","myfile.pdf",mylegend="a",\
                   xlim=(0,1),ylim=(0,1))
    
    INPUT PARAMETERS:
           x        - array with n elements       
           y        - array with n elements
           col      - color
           s        - size of points
           xlab     - r"$x$"
           ylab     - r"$y$"
           title    - The plot title
           filename - Filename e.g. "myfile.pdf"
           leg      - Legend
           legloc   - E.g. "upper right", "center right"
           xlim     - E.g. (0,1)
           ylim     - E.g. (0,1)
           logx     - True/False
           logy     - True/False
           show     - Show window
           invertx  - Invert the x axis
    
    OPTIONAL INPUT:
           See above
    
    OUTPUT PARAMETERS:
           None
    
    NOTES:
    
    MODIFICATION HISTORY:
           Coded by G. K. Stefansson - date, April 26, 2015
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,title=title)
    adjustprops = dict(left=0.19,bottom=0.15,right=0.92,top=0.9,wspace=0.,hspace=0.2)
    fig.subplots_adjust(**adjustprops)
    if xlab!="":
        ax.set_xlabel(xlab,size="x-large")
    if ylab!="":
        ax.set_ylabel(ylab,size="x-large")
    ax.minorticks_on()
    ax.grid()
    if xlim != 0:
        ax.set_xlim(xlim)
    if ylim != 0:
        ax.set_ylim(ylim)
    #PLOT
    ax.scatter(x,y,color=col,s=s,label=leg)
    if legloc!="":
        ax.legend(loc=legloc,prop={'size':16})
    else:
        ax.legend(loc="upper center",prop={'size':16})
    if logx == True:
        ax.set_xscale("log")
    if logy == True:
        ax.set_yscale("log")
    if invertx==True:
        plt.gca().invert_xaxis()
    if filename !="":
        fig.savefig(filename)
        print("Plot saved to file, "+filename)
    if show==False:
        plt.close()
    if show==True:
        fig.show()

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

def plot_param_table(ltable,rtable=None,ax=None,fontsize=8,linespacing=0.05,xmean=0.45):
        """
        Plot a labels and values table on a matplotlib axis
        """
        if ax==None:
            fig, ax = plt.subplots()
            
        ########################################
        ##### TEXT - PLANET
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        yt = 0.975
        if rtable is None:
            for l in ltable:
                ax.annotate(l, xy=(xmean+0.05, yt), xycoords="axes fraction", ha='right', fontsize=fontsize)
                yt -= linespacing
        else:
            for l,r in zip(ltable, rtable):
                ax.annotate(l, xy=(xmean+0.05, yt), xycoords="axes fraction", ha='right', fontsize=fontsize)
                ax.annotate(r, xy=(xmean+0.10, yt), xycoords="axes fraction", fontsize=fontsize)
                yt -= linespacing
        ax.axis("off")
