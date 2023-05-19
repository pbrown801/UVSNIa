;;;;; given a spectrum, output the spectrophotometry and count rates in UVOT system
;;;;; as well as the number of counts that would be in each filters wavelength range
;;;;; (to deal with overlapping filters)

pro uvot_speccounts, spectrum, specmags, speccounts, w2allcounts, m2allcounts, w1allcounts

h  = 1D*6.6260755E-27
c  = 1D*2.99792458E18
hc = 1D*1.986447E-8
e  = 1D*2.71828
k  = 1D*1.380658E-16

;;Read in the filter Effective Area curves
; can be obtained over the web as well
; on orbit calibrations of the UVOT optical filters as presented in Poole et al. 2008
; The curves are effective areas with units of cm^2
;
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/V_P08.txt",copyfile='V_P08.txt') 
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/B_P08.txt",copyfile='B_P08.txt') 
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/U_P08.txt",copyfile='U_P08.txt') 
;
; Breeveld et al. 2011 updated the UV filters
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/UVW1_B11.txt",copyfile='UVW1_B11.txt') 
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/UVM2_B11.txt",copyfile='UVM2_B11.txt') 
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/UVW2_B11.txt",copyfile='UVW2_B11.txt') 
;
; these are "red-leak corrected" filters similar to that used in Brown et al. 2010
; but using the final version of the Breeveld et al. 2011 filter curves
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/UVW1_B11rc.txt",copyfile='UVW1_B11rc.txt') 
; filters=webget("http://people.physics.tamu.edu/pbrown/UVOT/UVW2_B11rc.txt",copyfile='UVW2_B11rc.txt') 
;
; put the filters in a folder with a shortcut defined as $UVOTFILTERS
;

; old naming convention
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
readcol,"$SNSCRIPTS/B_UVOT.txt",lambda,B_EA,/silent
readcol,"$SNSCRIPTS/U_UVOT.txt",lambda,U_EA,/silent
readcol,"$SNSCRIPTS/UVW1_2010.txt",lambda,W1_EA,/silent
readcol,"$SNSCRIPTS/UVM2_2010.txt",lambda,M2_EA,/silent
readcol,"$SNSCRIPTS/UVW2_2010.txt",lambda,W2_EA,/silent
readcol,"$SNSCRIPTS/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"$SNSCRIPTS/UVW2-rc.txt",lambda,W2rc_EA,/silent

readcol,"$UVOTFILTERS/V_P08.txt",lambda,V_EA,/silent
readcol,"$UVOTFILTERS/B_P08.txt",lambda,B_EA,/silent
readcol,"$UVOTFILTERS/U_P08.txt",lambda,U_EA,/silent
readcol,"$UVOTFILTERS/UVW1_B11.txt",lambda,W1_EA,/silent
readcol,"$UVOTFILTERS/UVM2_B11.txt",lambda,M2_EA,/silent
readcol,"$UVOTFILTERS/UVW2_B11.txt",lambda,W2_EA,/silent
readcol,"$UVOTFILTERS/UVW1_B11rc.txt",lambda,W1rc_EA,/silent
readcol,"$UVOTFILTERS/UVW2_B11rc.txt",lambda,W2rc_EA,/silent

;;; the wavelength regions below represent the intersections
;;; of the filter curves when normalized by the effective area
;;; under the curve, effectively deweighting broader filters
w2lambda=where(lambda ge 1600 and lambda lt 2030)
m2lambda=where(lambda ge 2030 and lambda lt 2460)
w1lambda=where(lambda ge 2460 and lambda lt 3050)
uulambda=where(lambda ge 3050 and lambda lt 3870)
bblambda=where(lambda ge 3870 and lambda lt 4960)
vvlambda=where(lambda ge 4960 and lambda le 8000)

;;Read in the spectrum and interpolate to the Filter Curves

; check to see if it is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse
;zero=where(finite(sp_flux) eq 0)
;sp_flux[zero]=0.0

;  note: this assumes a flux spectrum versus wavelength in Angstroms
flux=interpol(sp_flux,sp_wave,lambda)

;calculate the count rates and magnitudes from the spectrum
v_counts=V_EA*flux*(10*lambda/hc)
b_counts=B_EA*flux*(10*lambda/hc)
u_counts=U_EA*flux*(10*lambda/hc)
w1_counts=W1_EA*flux*(10*lambda/hc)
m2_counts=M2_EA*flux*(10*lambda/hc)
w2_counts=W2_EA*flux*(10*lambda/hc)
w1rc_counts=W1rc_EA*flux*(10*lambda/hc)
w2rc_counts=W2rc_EA*flux*(10*lambda/hc)

v_spec=-2.5*alog10(total(v_counts))+17.89
b_spec=-2.5*alog10(total(b_counts))+19.11
u_spec=-2.5*alog10(total(u_counts))+18.34
w1_spec=-2.5*alog10(total(w1_counts))+17.44
m2_spec=-2.5*alog10(total(m2_counts))+16.85
w2_spec=-2.5*alog10(total(w2_counts))+17.38
w1rc_spec=-2.5*alog10(total(w1rc_counts))+17.34
w2rc_spec=-2.5*alog10(total(w2rc_counts))+17.25

; compute the fraction of the uvw1 and uvw2 counts coming through the red leak corrected versions
uvfrac=[total(w1rc_counts)/total(w1_counts),total(w2rc_counts)/total(w2_counts)]

speccounts=[total(w2_counts), total(m2_counts), total(w1_counts), total(u_counts), total(b_counts), total(v_counts), total(w1rc_counts), total(w2rc_counts)]
specmags=[w2_spec, m2_spec, w1_spec, u_spec, b_spec, v_spec, w1rc_spec, w2rc_spec]


;  determine the amount of counts coming through the different regions of the filter curves 
;  in order to correct for overlap of the filters
w2allcounts=[total(w2_counts[w2lambda]),total(w2_counts[m2lambda]),total(w2_counts[w1lambda]),total(w2_counts[uulambda]),total(w2_counts[bblambda]),total(w2_counts[vvlambda])]

m2allcounts=[total(m2_counts[w2lambda]),total(m2_counts[m2lambda]),total(m2_counts[w1lambda]),total(m2_counts[uulambda]),total(m2_counts[bblambda]),total(m2_counts[vvlambda])]

w1allcounts=[total(w1_counts[w2lambda]),total(w1_counts[m2lambda]),total(w1_counts[w1lambda]),total(w1_counts[uulambda]),total(w1_counts[bblambda]),total(w1_counts[vvlambda])]

end



pro SN2011by_mangle

magfile='$SWIFT_FINALDATA/PhotArchive/SN2005cf_uvotB11.dat'


uvotmags=[ 16.9140 , 18.3670 , 15.0640 , 13.3900 , 13.6390 , 13.6250]

spectrum='SN2011by_0509_full.dat'

zeropoints=[17.38,16.85,17.44,18.34,19.11,17.89]
uvotcounts=10.0^((zeropoints-uvotmags)/2.5)

; read this in just to set the lambda scale for the interpolation of the spectrum
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent

;;Read in the Spectra and interpolate to the
;;Filter Curves

; check to see if it is a string (ie a filename) or an array
s=size(spectrum)
if (s[1] eq 7) then begin
	readcol,spectrum,sp_wave,sp_flux,/silent
endif else begin
	sp_wave=spectrum[0,*]
	sp_flux=spectrum[1,*]
endelse

originalflux=interpol(sp_flux,sp_wave,lambda)
newspec=originalflux
newflux=originalflux

uvot_speccounts, spectrum, specmags, speccounts, w2allcounts, m2allcounts, w1allcounts

	filterratios=2.512^(specmags[0:5]-uvotmags)

plot, lambda, originalflux, /ylog
newflux[where(lambda lt 3000)]=newflux[where(lambda lt 3000)]/filterratios[1]
newflux[where(lambda ge 3000)]=newflux[where(lambda ge 3000)]/filterratios[4]
oplot, lambda, newflux

uvot_speccounts, [[lambda],[newflux]], specmags, speccounts, w2allcounts, m2allcounts, w1allcounts

uvot_speccounts, [[lambda],[originalflux]], specmags, speccounts, w2allcounts, m2allcounts, w1allcounts
uvot_speccounts, [[sp_wave],[sp_flux]], specmags, speccounts, w2allcounts, m2allcounts, w1allcounts
uvot_speccounts, spectrum, specmags, speccounts, w2allcounts, m2allcounts, w1allcounts


stop
end