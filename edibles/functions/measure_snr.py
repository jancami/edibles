def measure_snr(wave, flux, do_plot=False):

    """ call:
      
        snr = measure_snr(wave_array, flux_array)
 
    """
    #split in blocks of 1 Angstrom.
    block_size = 1.0

    xmin = wave[0]
    xmax = xmin + block_size

    SNR =[]
    LAM =[]

    while xmax < wave[-1]:
        idx = np.where( (wave > xmin) & (wave < xmax) )
        wave_block = wave[idx]
        flux_block = flux[idx]
        if (np.nanmean(flux_block) > 0.0):
            #fit polynomial (1st order) to "normalise"
            #c1 = np.polyfit(wave_block,flux_block,1)
            #p1 = np.poly1d(c1)
            #continuum = p1(wave_block)
            #flux_block = flux_block/continuum
            #compute stddev for each block
            sigma_block = np.nanmean(flux_block)/np.nanstd(flux_block)
            SNR.append(sigma_block)
            LAM.append(xmin+(xmax-xmin)/2.0)
        xmin = xmax.copy()
        xmax = xmin + block_size
    if (do_plot == True):
       plt.plot(LAM,SNR)
       plt.plot(LAM,smooth(SNR,20))
    
    maxsnr=np.nanmax(smooth(SNR,20))
    
    return maxsnr, SNR, LAM

   
