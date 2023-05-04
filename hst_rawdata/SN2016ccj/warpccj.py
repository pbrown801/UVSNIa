def countsin_sedout(counts_array):

    print('printing counts within')
    print(counts_array)

    import numpy as np
    import matplotlib.pyplot as plt
    from operator import truediv
#######

    h = 6.6260755e-27
    c = 2.99792458e18
    hc = h*c

    files = ['filters/UVW2_2010.txt', 'filters/UVM2_2010.txt', 'filters/UVW1_2010.txt', 'filters/U_UVOT.txt', 'filters/B_UVOT.txt', 'filters/V_UVOT.txt']

    filter_WL = []
    filter_A = []

    for item in files:

        f = open(item,'r')

        filter_lambda = []
        filter_area = []
        for line in f:
            line = line.rstrip()
            column = line.split()
            wavelen = column[0]
            area = column[1]
            filter_lambda.append(float(wavelen))
            filter_area.append(float(area))

        filter_lambda = np.asarray(filter_lambda,dtype=float)
        filter_area = np.asarray(filter_area,dtype=float)
        
 #       nonzero = np.where(filter_area > 0.0)
        
 #       filter_lambda = filter_lambda[nonzero]
 #       filter_area = filter_area[nonzero]

        filter_WL.append(filter_lambda)
        filter_A.append(filter_area)

        f.close()
###    print(filter_area)
###    plt.plot(filter_lambda,filter_area)
###    plt.show

###    print("something")
###    wait = input("PRESS ENTER TO CONTINUE.")
###    print("something")    
### Starting with the factors associated with Pickles


    factor = [5.77E-16, 7.47E-16, 4.06E-16, 1.53E-16, 1.31E-16, 2.61E-16]
    factor = np.asarray(factor, dtype=float)

    new_WL = [1928,2246,2600,3465,4392,5468]
    flux_top=[]
    new_counts=[12341,12341234,123412341,1234134,12341234,1234132]
    flag = [0,0,0,0,0,0]


    filter_array = np.array([filter_A[0],filter_A[1],filter_A[2],filter_A[3],filter_A[4],filter_A[5]])

    filter_wave = np.array([filter_WL[0],filter_WL[1],filter_WL[2],filter_WL[3],filter_WL[4],filter_WL[5]])

   
    ### Iterating over until we are within 10% of the input counts
    ### Will print the value for the flux for each iteration

    while sum(flag) != 6:
        flux_top=[]
        for count in range(0,len(counts_array)):

            flux_top.append(counts_array[count]*factor[count])

        for x in range(len(flux_top)):
            
            new_spec = np.interp(filter_wave[x], new_WL, flux_top)
            
            new_counts[x] = np.trapz(new_spec*filter_array[x]*filter_wave[x]/hc,filter_wave[x])

            factor = map(truediv, flux_top, new_counts)
            # print(factor)
	    print(new_counts)
            if abs(new_counts[x] - counts_array[x]) <= 0.01*counts_array[x]:
                flag[x] = 1
            else:
                flag[x] = 0
#        print flux_top
    print(new_counts)
#    return(flux_top);
    return(new_spec);





import numpy as np
from countsin_sedout import *

filters = ['UVW2_2010.txt','UVM2_2010.txt','UVW1_2010.txt','U_UVOT.txt','B_UVOT.txt','V_UVOT.txt']
newfilters=[]
for f in range(len(filters)):
	str1='filters/'
	str2=filters[f]
	newstr=str1+str2
	newfilters.append(newstr)

obs_mags = [7.80,7.95,8.01,6.50,6.20,6.00] #arbitrary
vega_wave,vega_flux = np.loadtxt('spectra/vega.dat',dtype = float, usecols = (0,1), unpack=True)

def pivot_wavelength(Filter):
    
    filter_wave,filter_tp = np.loadtxt(Filter, dtype = float, usecols=(0,1), unpack=True)

    numerator = np.trapz(filter_tp*filter_wave,filter_wave)
    denominator = np.trapz(filter_tp/filter_wave,filter_wave)
    
    pivot_lambda = np.sqrt(numerator/denominator)

    return pivot_lambda




def vegaspecphot(wavez,fluxz,Filter):

    h = 6.6260755e-27
    c = 2.99792458e18
    hc = h*c #units of erg*A

    filter_lambda,filter_area = np.loadtxt(Filter,comments='#',usecols=(0,1), unpack=True)

    nonzero = np.where(filter_area > 0.0)
        
    filter_lambda = filter_lambda[nonzero]
    filter_area = filter_area[nonzero]


    ##############   calculate vega zeropoint for every filter from vega spectrum

    in_lambda_range = np.where((vega_wave>=min(filter_lambda))&(vega_wave<=max(filter_lambda)))
    interpolated_flux = np.interp(filter_lambda,vega_wave[in_lambda_range[0]],vega_flux[in_lambda_range[0]])
    zeropoint = round(-2.5*np.log10(np.trapz(filter_area*interpolated_flux*filter_lambda/hc,filter_lambda)),2)

    # Calculated magnitudes

    sp_ea = np.interp(wavez,filter_lambda,filter_area) ### spectrum effective area         
    counts = np.trapz(sp_ea*fluxz*wavez/hc,wavez) ### Integrating under the curve using numpy
    if counts > 0:           
        vegamag = -2.5*np.log10(counts) - zeropoint ### Calculated magnitudes
    return vegamag,counts





def compare_template(obs_mag_array,Filters,template_wave,template_flux):

    ###Input data types###
        #obs_mag_array : numpy array
        #Filters : string list
        #template_wave : numpy array
        #template_flux: numpy array

    #Import the filter details#
    specmag_array = []
    flux_ratios = []
    pivot_array = []
    counts_list = []
    for l in range(len(Filters)):
        filter_wave,filter_tp = np.loadtxt(Filters[l], dtype = float, usecols=(0,1), unpack=True)
        nonzero =  np.where(filter_tp > 0.0)

        filter_wave = filter_wave[nonzero]
        filter_tp = filter_tp[nonzero]

        spec_mag,counts = vegaspecphot(template_wave,template_flux,Filters[l])
        specmag_array.append(spec_mag)
        counts_list.append(counts)
                                    
        value = 10**(-.4*(obs_mag_array[l]-specmag_array[l]))
                     
        flux_ratios.append(value)
        pivot_array.append(pivot_wavelength(Filters[l]))

    specmag_array = np.asarray(specmag_array)
    flux_ratios = np.asarray(flux_ratios)
    pivot_array = np.asarray(pivot_array)
    
    ratiocurve = np.interp(template_wave,pivot_array,flux_ratios)
    mangledspec = ratiocurve*template_flux

    sed = countsin_sedout(counts_list)
    sed_wave = np.linspace(1600,8000,num = 641,endpoint =True)
    sed_interp = np.interp(template_wave,sed_wave,sed)

    #make a plot of the ratiocurve
##    
##    fig,ax1 = plt.subplots()
##
##    ##Configure plot axes for first plot
##    ax1.plot(template_wave,mangledspec,'k--')
##    ax1.set_xlabel('Wavelength')
##    ax1.set_ylabel('Mangled Spectrum Flux', color ='k')
##    ax1.tick_params('y',colors='k')
##
##    #Since I want to plot 2 different scaled plots on one plot, configure axes for second plot
##    
##    ax2 = ax1.twinx() #Both plots share same x axis
##    ax2.plot(template_wave,template_flux,'b--') #The plot with different color than first plot
##    ax2.set_ylabel('Template flux', color = 'b') #set a y label and give it a color
##    ax2.tick_params('y', colors='b') #Give y ticks of second plot a corresponding color for clarity
##    ax2.yaxis.tick_right()
##
##    ax3 = ax1.twinx()
##    ax3.plot(template_wave,sed_interp,'r-',label='SED')
##    ax3.set_yticks([])
##        
##    ax3.legend(loc = 'upper right')
##    plt.show()
    
    #Make flux ratio vs wavelength plot#
##    plt.plot(template_wave,ratiocurve,'b--')
##    plt.xlabel('wavelength')
##    plt.ylabel(r'Flux Ratio ($\frac{F_{obs}}{F_{spec}}$)')
##    plt.title('Flux Ratio v.s. Wavelength')
##    plt.show()
    return mangledspec,template_wave

mangled,wave = compare_template(obs_mags,filters,vega_wave,vega_flux)

###run vegaspecphot on the mangled spectrum###
w2mag = vegaspecphot(wave,mangled,filters[0])
m2mag = vegaspecphot(wave,mangled,filters[1])
w1mag = vegaspecphot(wave,mangled,filters[2])
umag = vegaspecphot(wave,mangled,filters[3])
bmag = vegaspecphot(wave,mangled,filters[4])
vmag = vegaspecphot(wave,mangled,filters[5])

mags = [w2mag[0],m2mag[0],w1mag[0],umag[0],bmag[0],vmag[0]]

flux_ratio2 = 10**(-.4*np.subtract(obs_mags,mags))
    
        