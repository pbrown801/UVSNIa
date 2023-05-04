
#function for weighting the overlap
def overlapweight(list1, list2):
   start = 0
   #global variables
   global cbinf, numbin, cbind, cbinu, cbinf2, numbin2, cbind2, cbinu2, overlap, floverlap, uncoverlap, dqoverlap, wv3, efflen, cbinw, cbinw2 
   #getting the largest number of indicies in the list
   if(len(list1)>len(list2)):
    efflen = len(list1)
   if(len(list2)>len(list1)):
    efflen = len(list2)
   if(len(list2)==len(list1)):
    efflen = len(list1) 
   #aranging the wavelength 
   if(min(list1)<min(list2)):
    wv3 = (np.arange(int(min(list1)), int(max(list2)), 10))
   if(min(list2) < min(list1)): 
    wv3 = (np.arange(int(min(list2)), int(max(list1)), 10))
   for i in range(len(wv3)):
      #w3 placeholder for current element of wv3, probably uncessary
      w3 = wv3[i]
      #for x in list1, get indicides of values between x and y, sum and average, then add to the thing, do the same with list2, then weight
        #binning data
      for x in range(start, efflen):
          print(str(list1[x]) + " " + str(w3))
          if(list1[x]<w3):
           cbinw += list1[x]
           cbinf += fl1[x]
           cbinu += unc1[x] #need to set as floats
           cbind += dq1[x] #need to get index of the wv to add same idx of flux and other data
           numbin += 1
           print("1")
           print(str(list2[x]) + " " + str(w3))  
          if(list2[x]<w3):
           cbinw2 += list2[x]
           cbinf2 += fl2[x]
           cbinu2 += unc1[x]
           cbind2 += dq1[x] #need to get index of the wv to add same idx of flux and other data
           numbin2 += 1  
           print("2")
      #start = int((numbin+numbin2)/2): trying to limit the time use of the for loop by changing the start index to the one it ended on last time, doesnt work for ? reason
      try:  
        #adding the (binned) overlapping data in order of list1, list 2
        #from list 1
        overlap.append(float(cbinw/numbin))
        floverlap.append(float(cbinf/numbin)) 
        uncoverlap.append(float(cbinu/numbin))
        dqoverlap.append(float(cbind/numbin))
        #from list 2
        overlap.append(float(cbinw2/numbin2))
        floverlap.append(float(cbinf2/numbin2)) 
        uncoverlap.append(float(cbinu2/numbin2))
        dqoverlap.append(float(cbind2/numbin2))
      except:
        #print("excepted")
        pass
   #begin weighting 
   print("weighting")
   print(overlap)
  #sorting
   #turning into arrays
   overlap = np.array(overlap)
   floverlap = np.array(floverlap)
   uncoverlap = np.array(uncoverlap)
   dqoverlap = np.array(dqoverlap)
   ind = overlap.argsort()
   #sorting by argsorted indicies
   overlap = overlap[ind]
   floverlap = floverlap[ind]
   uncoverlap = uncoverlap[ind]
   dqoverlap = dqoverlap[ind]
   #turning them back into lists? hopefully
   overlap = overlap.tolist()
   floverlap = floverlap.tolist()
   uncoverlap = uncoverlap.tolist()
   dqoverlap = dqoverlap.tolist()
   #begin weighting function
   if((max(list1) > min(list2)) & (min(list1) < min(list2))): #end of list1 and beginning of list2 overlap
    #because stuff is added to overlap in order, i could take it from overlap in order and put it into list1 and list2 arrays (binned), and then use those for the weighting?
    #find where overlap begins and ends               
    for x in range(int(min(overlap)),int(max(overlap))): #from the min of the overlap to the max of the overlap, find the list indicies?
        idxo1 = list1.index(x)
        idxo2 = list2.index(x)
       
        if(list1(idxo1) < median(overlap)):
         weighting2 = (list2(idxo2)/median(overlap))
         weighting1 = 1-weighting2 #(list1(idxo1)/median(overlap))
         weightedflux.append(weighting1*fl1[idxo1]+weighting2*fl2[idxo2])
       
        if(list1(idxo1) > median(overlap)):
         weighting1 = ((2*median(overlap)-list1(idxo1))/median(overlap))
         weighting2 = 1-weighting1#((2*median(overlap)-list2(idxo2))/median(overlap))
         weightedflux.append(weighting1*fl1[idxo1]+weighting2*fl2[idxo2])
   if((max(list2) > min (list2)) & (min(list2) < min(list1))): #end of list2 and beginning of list1 overlap    
    for x in range(int(min(overlap)),int(max(overlap))):
        idxo1 = list1.index(x)
        idxo2 = list2.index(x)
       
        if(list1(idxo1) < median(overlap)):
         weighting1 = (list1(idxo1)/median(overlap))
         weighting2 = 1-weighting1#(list2(idxo2)/median(overlap))
         weightedflux.append(weighting1*fl1[idxo1]+weighting2*fl2[idxo2])
       
        if(list1(idxo1) > median(overlap)):
         weighting2 = ((2*median(overlap)-list2(idxo2))/median(overlap))
         weighting1 = 1-weighting2#((2*median(overlap)-list1(idxo1))/median(overlap))    
         weightedflux.append(weighting1*fl1[idxo1]+weighting2*fl2[idxo2])
   return

def combinespectra(combfiles,opticalfile,filename):
	print('start')
	from re import X
	from matplotlib.pyplot import xlabel
	import matplotlib.pyplot as plt
	import numpy as np
	#import matplotlib.ticker as ticker
	import statistics
	from statistics import median
	import math
	from numpy.lib.function_base import append
	import astropy.units as u
	from specutils import Spectrum1D
	from specutils.manipulation import FluxConservingResampler
	from astropy.nddata import StdDevUncertainty
	import pandas as pd
	import os

	#size of the bins
	binsize = 10
	#min and max wv considered (in latter part)
	wvmin = 1600
	wvmax= 5650

	optwave = []
	optflux = []
	wavelength1 = []
	flux1 = []
	fluxmean = []
	wavelengthcompile = []
	fluxcompile = []
	uncertaintycompile = []
	dataqualitycompile = []
	dataquality1 = []
	wavelength = []
	wavelengthb = []
	flux = []
	dataquality = []

	#uncertainties for weights
	uncertainty1 = []
	uncertainty = []
	uncertaintyb = []
	#can have a string of spectra so once complete, it can run through it with all the files
	# create a third array which is from the minimum of wavelength1 (the smaller one?) to the max of wavelength2(bigger range?), with spaces of 10 units    
	#wavelength3 = np.arange(int(min(wavelength1)), int(max(wavelength2)), 10)

	#create empty flux array c:\Users\Evan\Desktop\Portfolio\Code\DrBProj\combinespectra.py", line 156
	flux3 = []
	fileprogress = 0
		
	#empty wavelength arrays
	wv1, wv2, wv3 = [], [], [] 
	#empty flux arrays
	fl1, fl2, weightedflux = [], [], []
	#empty uncertainity arrays
	unc1, unc2, unc3 = [], [], []
	#empty data quality arrays
	dq1, dq2, dq3 = [], [], []
	#declare overlap array so it can be used outside overlap function
	overlap, floverlap, uncoverlap, dqoverlap = [], [], [], [] 
	w3 = 0.0
	cbinw, cbinw2 = 0, 0
	cbinf, cbinf2 = 0, 0
	cbinu, cbinu2 = 0, 0
	cbind, cbind2 = 0, 0
	numbin, numbin2 = 0, 0
	efflen = 0
	#default weights
	weighting1 = .5
	weighting2 = .5

	#if there is more than one file with overlaps
	print("begin reading")
	if(len(combfiles) > 1):
		print("length > 1")
		#getting data from the files
		for x in range(len(combfiles)):
			print("file 1 working")
			with open(combfiles[x],'r') as f1:
				for line in f1:
					content = line.split()
					for num in range(0, (len(content)-3)):
						wv1.append(float(content[0]))
						fl1.append(float(content[1]))
						if(float(content[2]) != None):  #there may not be uncertainity or data quality in these files
							unc1.append(float(content[2])) #add a weighting of one to uncertainity, if spectra looks reasonable its probably okay.
						else:
							unc1.append(1)
						dq1.append(float(content[3]))
			print("file 2 working")
			with open(combfiles[x+1],'r') as f2:       #opening the second file in advance in order to compare the two datasets
				for line in f2:
					content = line.split()
					for num in range(0, (len(content)-3)):
						wv2.append(float(content[0]))
						fl2.append(float(content[1]))
						if(float(content[2]) != None):
							unc2.append(float(content[2]))
						else:
							unc2.append(1)
						dq2.append(float(content[3]))
				if(fileprogress < 1):
						#for the first 2 files, compare data from the files
					print("begin function ")
					overlapweight(wv1,wv2) #function to generate weighted flux
				if(fileprogress > 1): #for the next files, compare one data set from files, and one we already have genereated based off the data from previous files
					overlapweight(wv3, wv2)
	print("end")
  
	#flux3 = np.array(flux3) 
	#a case for if there is only one file
	if(len(combfiles) == 1): 
		for x in range(len(combfiles)):
			with open(x,'r') as f1:
				for line in f1:
					content = line.split()
					for num in range(0, (len(content)-3)):
						wv1.append(float(content[0]))
						fl1.append(float(content[1]))
						unc1.append(float(content[2]))
						dq1.append(float(content[3])) 
					wv3 = np.arange(int(min(wv1)), int(max(wv1)), 10)
					for i in range(len(wv3)):
						#w3 placeholder
						w3 = wv3[i]
						idx1, idx2 = None, None
						#checking if there is flux at the index of w3, if there is at both, average them, 
						#if there is one but not another, it just takes that flux
						try:
							idx1 = wv1.tolist().index(w3)
						except:
							pass
						if idx2 is None:
							flux3.append(fl1[idx1])
							unc3.append(unc1[idx2])
							dq3.append(dq1[idx2])
							continue

	#getting data from the optical file
	with open(opticalfile,'r') as f:
		for line in f:
			print(line)
			content = line.split()
			optwave.append(float(content[0]))
			optflux.append(float(content[1]))

	#getting data from files, putting into arrays
	for file in files:
		with open(file,'r') as f:
			for line in f:
				content = line.split()
				for num in range(0, (len(content)-3)):
					wavelengthcompile.append(float(content[0]))
					fluxcompile.append(float(content[1]))
					uncertaintycompile.append(float(content[2]))
					dataqualitycompile.append(float(content[3]))
       
	#making data into floats
	#taking only data within a specified wv range
	length = wavelengthcompile
	for x in range(0,len(length)):
    		if((wavelengthcompile[x] > wvmin) and (wavelengthcompile[x] < wvmax) and (dataqualitycompile[x] != 16912)):
        		wavelength.append(wavelengthcompile[x])
        		flux.append(fluxcompile[x])
        		uncertainty.append(uncertaintycompile[x])
        		dataquality.append(dataqualitycompile[x])

	#sorting them by wavelength
	#adding the other data
	wavelength.append(wv3)
	flux.append(weightedflux)
	uncertainty.append(unc3)
	dataquality.append(dq3)
	#turning them into arrays
	wavelength = np.array(wavelength)
	flux = np.array(flux)
	uncertainty = np.array(uncertainty)
	dataquality = np.array(dataquality)
	#arranging them
	#using sorting, basing them on the index of the wavelengths
	inds = wavelength.argsort()
	flux = flux[inds]
	wavelength = wavelength[inds]
	uncertainty = uncertainty[inds]
	dataquality = dataquality[inds]
					

	#getting min max of wavelengths in the binsize
	min=math.ceil(min(wavelength)/binsize)*binsize
	max=math.floor(max(wavelength/binsize))*binsize
	print(min, max)
	#creating array of evenly spaced bins for wavelength
	newwave=[]
	for w in range(min, max, binsize):
    		newwave.append(w)

	#using Spectrum1D to get the flux right, converting to angstroms
	input_spectra = Spectrum1D( flux=flux*u.erg/u.s/u.cm/u.cm/u.angstrom, spectral_axis=np.array(wavelength)*u.angstrom, uncertainty=StdDevUncertainty(uncertainty) )

	# , uncertainty=uncertainty 
	#rebinning the flux using a flux conserving resapmler
	flux_rebinned = FluxConservingResampler()
	newflux = flux_rebinned(input_spectra, newwave*u.angstrom) 
	#plotting the bins of wavelength and flux
	plt.plot(newwave, newflux.flux, 'b-', label = "HST", alpha = 1, linewidth = .9)
	#plotting optical view stuff
	plt.plot(optwave, optflux, 'c-', label = "LCO", alpha = 1, linewidth = .7)
	plt.yscale('log')
	plt.legend()
	#axes = plt.axes()
	#axes.set_xlim([1600,9000])
	#axes.set_ylim([10**(-20),10**(-16)])
	plt.savefig((str(filename) + ".png"))
	plt.show()

	#potential other way to save text file
	#data = np.column_stack([newwave, newflux])
	#datafile_path = "combined_spectra.txt"
	#np.savetxt(datafile_path , data)

	#writing data to file
	f = open(filename + ".dat", "x")
	for x in range(0, len(newwave)):
  		f.write(str(newwave[x]) + " " + str(float(newflux[x].flux.value)))
  		f.write("\n")  
	f.close()

"""
 # try:
      #  idx2 = int(list2.index(w3))
      #  print(idx1)
      #except:
      #   pass
      if(idx1 is None and (idx2 is not None)):
        print(idx2)
        weightedflux.append(fl2[idx2])
        unc3.append(unc2[idx2])
        dq3.append(dq2[idx2])
        continue
      if(idx2 is None and(idx1 is not None)):
        weightedflux.append(fl1[idx1])
        unc3.append(unc1[idx2])
        dq3.append(dq1[idx2])
        continue
      #appending the fluxes in order of flux from list 1, flux from list 2, in order to weight them later
      if(idx1 is not None and idx2 is not None):
       flux3.append(fl1[idx1]) #have something thats equal to 1 where it doesnt overlap with other spectrum, and decreases linear to 0 at its furthurst wavelength bin, so their weighting things will add up to one, if they both have uncertainties
       flux3.append(fl2[idx2])
       unc3.append(.5*(unc1[idx1]+unc2[idx2]))                  #for weighting, find the overlap, make weight add up to one, .5 and .5 at center of overlap
       dq3.append(.5*(dq1[idx1]+dq2[idx2]))
       #adding the overlapping wavelengths in the bins

"""
