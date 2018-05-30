import scipy.signal as signal
def low_pass_butter(x,filt_order=3.,cutoff_freq=0.07):
    """
    Implements a butter filter
    
    INPUT:
        x - array to filter
        filt_order - 2, 3 ? 
        cutoff_freq = 0.01
    
    RETURNS:
        Low pass array of same length
        
    EXAMPLE:
        # Butter filter
        N  = 3    # Filter order
        Wn = 0.07 # Cutoff frequency
        B, A = signal.butter(N, Wn, output='ba')
        tempf = signal.filtfilt(B,A, aa.df_all["X(FITS)_T1"])
        plt.plot(aa.df_all["X(FITS)_T1"])
        plt.plot(tempf)
    """
    B, A = signal.butter(filt_order, cutoff_freq, output='ba')
    x_low_pass = signal.filtfilt(B,A,x)
    return x_low_pass
