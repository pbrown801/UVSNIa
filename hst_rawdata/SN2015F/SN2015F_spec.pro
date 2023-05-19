pro SN2015F_spec



spec1name='SN2015F'
spec1file='SN2015F_full_20150326.dat'
distance1=31.511   ;  Riess et al. 2016
distance1err=0.053
spec1_ebv=3.1
redshift1=0.004890
;;;;;;


;;;;;;;
spec2name='SN2011fe'
spec2file='../SN2011fe/SN2011fe_sullivan/ptf11kly_20110910.obs.dat'
distance2=29.135 ; Riess et al. 2016
distance2err=0.035 ; 
color2='blue'
redshift2=0.000804

;;;;;;;
spec3name='SN2011by'
spec3file='../SN2011by/SN2011by_0509_full.dat'
distance3=31.587  ; Riess et al. 2016
distance3err=0.070 ; 
color3='dark green'
redshift3= 0.002843 
ytop=0.5

;;;;;;;;;;;;

pjb_uvotspec_all, spec1file, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres, fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux, zeropoints=zeropoints, effwavelength=effwavelength, sp_wave=sp_wave, sp_flux=sp_flux,smooth_wave=smooth_wave, smooth_flux=smooth_flux

spec1_wave=sp_wave
spec1_flux=sp_flux
spec1_uvotres=spectrum_uvotres
spec1_smooth=smooth_flux

spec1_mags=mag_array

pjb_uvotspec_all, spec2file, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres, fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux, zeropoints=zeropoints, effwavelength=effwavelength, sp_wave=sp_wave, sp_flux=sp_flux,smooth_wave=smooth_wave, smooth_flux=smooth_flux

spec2_wave=sp_wave
spec2_flux=sp_flux
spec2_uvotres=spectrum_uvotres
spec2_smooth=smooth_flux

spec2_mags=mag_array

;;;;;;;;;;;;;;
pjb_uvotspec_all, spec3file, mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres, fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux, zeropoints=zeropoints, effwavelength=effwavelength, sp_wave=sp_wave, sp_flux=sp_flux,smooth_wave=smooth_wave, smooth_flux=smooth_flux

spec3_wave=sp_wave
spec3_flux=sp_flux

spec3_uvotres=spectrum_uvotres
spec3_smooth=smooth_flux

spec3_mags=mag_array

spec3_mags[5]=12.95
;;;;;;;;;;;;;;
max=floor(alog10(max([spec1_smooth,spec2_smooth,spec3_smooth])))

ebv2=(spec1_mags[4]-spec1_mags[5]) - (spec2_mags[4]-spec2_mags[5])
;ebv3=(spec1_mags[4]-spec1_mags[5]) - (spec3_mags[4]-spec3_mags[5])
ebv3=(spec1_mags[4]-spec1_mags[5]) - (12.92-12.95)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fm_unred, spec1_wave, spec1_flux, 0.556/3.1, funred, R_V = 3.1


;plot, spec1_wave/1.00489, funred*4.0*3.14159265*(10.0^(distance1/5.0))^2.0, xrange=[2000,6000]
;oplot, spec2_wave, spec2_flux*4.0*3.14159265*(10.0^(distance2/5.0))^2.0

;plot, spec1_wave/1.00489, funred*4.0*3.14159265*(10.0^(distance1/5.0))^2.0, xrange=[2000,6000]
;oplot, spec3_wave, spec3_flux*4.0*3.14159265*(10.0^(distance3/5.0))^2.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xdata=[5,6,7]
ydata=[4,5,6]

nxplots=1
nyplots=1
nplots=nxplots*nyplots

; from http://www.iluvatar.org/~dwijn/idlfigures

; to get R_lambda to work
;!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a single journal column width
xsize = 8.8

; putting these in real size (centimeters) 
; so they are the same regardless of the number of plots
wall = 0.04*8.8 ; these values were originally made for the 8.8 width
margin=0.16*8.8 ; a little big if the y-axis label is single digits 
                ; but accomodates more digits

; width and height of single plots
a = (xsize - margin - wall)/nxplots
; for square plots
b = a 
; for golden ratio plots
b = a * 2d / (1 + sqrt(5))

ysize= b * nyplots + margin + wall

ysize = (margin + nyplots*(b + wall ) )
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

xrange0=[1500,6000]
nxticks=9

xtitle0='Wavelength [Angstroms]'

fontsize=12


readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent

figurename=spec1name+'_'+spec2name+'.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=0
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='F!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,ytop], ystyle=1, xrange=[1500,5500], xstyle=1, $
xticks=4, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

fm_unred, spec1_wave, spec1_flux, 0.18, funred, R_V = 3.1

cgoplot,  spec1_wave/(1.0+redshift1), funred*2.512^(distance1), color='red', line=0
cgoplot,  spec2_wave/(1.0+redshift2), spec2_flux*2.512^(distance2)*0.7, color=color2, line=1
cgoplot,  spec3_wave/(1.0+redshift3), spec3_flux*2.512^(distance3)*1.0, color=color3, line=2

al_legend, [spec2name+' x0.7', spec1name, spec3name],   linestyle=[1,0,2],  box=0, $
pos=[0.34,0.95], color=[cgcolor(color2),cgcolor('red'),cgcolor(color3)],/norm, charsize=1


device, /close
SET_PLOT, 'X'


;;;;;
figurename=spec1name+'_'+spec2name+'_walker.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=0
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='F!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,0.4], ystyle=1, xrange=[1500,5500], xstyle=1, $
xticks=4, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

fm_unred, spec1_wave, spec1_flux, 0.18, funred, R_V = 3.1

cgoplot,  spec1_wave/(1.0+redshift1), funred*2.512^(distance1), color='red', line=0
cgoplot,  spec2_wave/(1.0+redshift2), spec2_flux*2.512^(distance2)*0.7, color=color2, line=2

scale=11.0*10^(12.0)
walkerpath='/Users/pbrown/Desktop/SN/sniamodels/walker/'
walkerspec=walkerpath+'spct-9.4-0.20.dat
readcol,walkerspec,wave,flux,/silent
cgoplot, wave, flux*scale, color='green', linestyle=1




walkerpath='/Users/pbrown/Desktop/SN/sniamodels/walker/'
walkerspec=walkerpath+'spct-9.4-2.00.dat
readcol,walkerspec,wave,flux,/silent
cgoplot, wave, flux*scale, color='orange', linestyle=3


al_legend, ['W -0.2','W -2.0'],   linestyle=[1,3],  box=0, $
pos=[3000,0.1], color=['green','orange'], charsize=1




al_legend, [spec2name+' x0.7', spec1name],   linestyle=[2,0],  box=0, $
pos=[1500,0.4], color=[cgcolor(color2),cgcolor('red')], charsize=1


device, /close
SET_PLOT, 'X'

stop

;;;;;

figurename=spec1name+'_'+spec2name+'_R.eps'


SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=0
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='R!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,20], ystyle=1, xrange=[1500,5500], xticks=4, xstyle=1, $
xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

EBV=ebv2

rv=3.1
	lambda_av_mw31=sne_mw_reddening(lambda,Ebv,rv=rv)
rv=1.7	
	lambda_av_mw17=sne_mw_reddening(lambda,Ebv,rv=rv)
	lambda_av_mw14=sne_mw_reddening(lambda,Ebv,rv=1.4)
	lambda_av_mw11=sne_mw_reddening(lambda,Ebv,rv=1.1)
	lambda_av_mw26=sne_mw_reddening(lambda,Ebv,rv=2.6)
	lambda_av_smc=sne_smc_reddening(lambda,Ebv)
	lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv)
	lambda_av_salt=sne_salt_reddening(lambda,Ebv)
	lambda_av_comb=sne_goobarlmc_reddening(lambda,0.5*Ebv)+sne_mw_reddening(lambda,0.5*Ebv,rv=2.6)


oplot, lambda, 	lambda_av_mw31/EBV, linestyle=0
oplot, lambda, 	lambda_av_mw17/EBV, linestyle=1
oplot, lambda, 	lambda_av_smc/EBV, linestyle=2
oplot, lambda, 	lambda_av_cslmc/EBV, linestyle=3
;oplot, lambda, 	lambda_av_salt/EBV, linestyle=4
;oplot, lambda, 	lambda_av_comb/EBV, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14/EBV, linestyle=0, thick=4

compare=where(lambda gt max([spec1_wave[0],spec2_wave[0]]) and lambda lt min([max(spec1_wave),max(spec2_wave)]) ) 

comparesmooth=where(smooth_wave gt max([spec1_wave[0],spec2_wave[0]]) and lambda lt min([max(spec1_wave),max(spec2_wave)]) ) 

cgoplot, smooth_wave[comparesmooth], ( (distance2-distance1) +2.5*alog10(spec2_smooth[comparesmooth]/spec1_smooth[comparesmooth]))/Ebv, thick=5, color='red'

;cgoplot, lambda[compare], ( (distance2-distance1) +2.5*alog10(spec2_uvotres[compare]/spec1_uvotres[compare]))/Ebv, thick=5, color='red'

;al_legend, ['MW R_V=3.1', 'CCM R_V=1.7', 'SMC', 'G08LMC', 'SALT'], background_color='white',pos=[0.42,0.94], /norm, charsize=1, box=0
   
al_legend, ['MW R_V=3.1', 'CCM R_V=1.7', 'SMC', 'G08LMC'],   linestyle=[0,1,2,3],  background_color='white', pos=[0.4,0.93], /norm, charsize=1, box=0


device, /close
SET_PLOT, 'X'


$open SN2015F_SN2011fe_R.eps
;;;;;;;;;;;;;;
figurename=spec1name+'_'+spec3name+'_R.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=0
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='R!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,20], ystyle=1, xrange=[1500,5500], xstyle=1, $
xticks=4, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

EBV=ebv3

rv=3.1
	lambda_av_mw31=sne_mw_reddening(lambda,Ebv,rv=rv)
rv=1.7	
	lambda_av_mw17=sne_mw_reddening(lambda,Ebv,rv=rv)
	lambda_av_mw14=sne_mw_reddening(lambda,Ebv,rv=1.4)
	lambda_av_mw11=sne_mw_reddening(lambda,Ebv,rv=1.1)
	lambda_av_mw26=sne_mw_reddening(lambda,Ebv,rv=2.6)
	lambda_av_smc=sne_smc_reddening(lambda,Ebv)
	lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv)
	lambda_av_salt=sne_salt_reddening(lambda,Ebv)
	lambda_av_comb=sne_goobarlmc_reddening(lambda,0.5*Ebv)+sne_mw_reddening(lambda,0.5*Ebv,rv=2.6)


oplot, lambda, 	lambda_av_mw31/EBV, linestyle=0
oplot, lambda, 	lambda_av_mw17/EBV, linestyle=1
oplot, lambda, 	lambda_av_smc/EBV, linestyle=2
oplot, lambda, 	lambda_av_cslmc/EBV, linestyle=3
;oplot, lambda, 	lambda_av_salt/EBV, linestyle=4
;oplot, lambda, 	lambda_av_comb/EBV, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14/EBV, linestyle=0, thick=4

compare=where(lambda gt max([spec1_wave[0],spec3_wave[0]]) and lambda lt min([max(spec1_wave),max(spec3_wave)]) ) 

comparesmooth=where(smooth_wave gt max([spec1_wave[0],spec3_wave[0]]) and lambda lt min([max(spec1_wave),max(spec3_wave)]) ) 

cgoplot, smooth_wave[comparesmooth], ( (distance3-distance1) +2.5*alog10(spec3_smooth[comparesmooth]/spec1_smooth[comparesmooth]))/Ebv, thick=5, color='red'

;cgoplot, lambda[compare], ( (distance3-distance1) +2.5*alog10(spec3_uvotres[compare]/spec1_uvotres[compare]))/Ebv, thick=5, color='red'

al_legend, ['MW R_V=3.1', 'CCM R_V=1.7', 'SMC', 'G08LMC'],   linestyle=[0,1,2,3],  background_color='white', pos=[0.4,0.93], /norm, charsize=1, box=0



device, /close
SET_PLOT, 'X'


$open SN2015F_SN2011by_R.eps

;;;;;;;;;;;
print, 'final stop'
stop
end


