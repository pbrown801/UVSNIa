pro comparespec_15F_11feby


spec1name='SN2015F'
spec1file='SN2015F_full_20150326.dat'
distance1=31.511   ;  Riess et al. 2016
distance1err=0.053
spec1_ebv=3.1
;;;;;;


;;;;;;;
spec2name='SN2011fe'
spec2file='../SN2011fe/SN2011fe_sullivan/ptf11kly_20110910.obs.dat'
distance2=29.135 ; Riess et al. 2016
distance2err=0.035 ; 
color2='blue'

;;;;;;;
spec3name='SN2011by'
spec3file='../SN2011by/SN2011by_0509_full.dat'
distance3=31.587  ; Riess et al. 2016
distance3err=0.070 ; 
color3='green'

ytop=1

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


plot, spec1_wave/1.00489, funred*4.0*3.14159265*(10.0^(distance1/5.0))^2.0, xrange=[2000,6000]
oplot, spec2_wave, spec2_flux*4.0*3.14159265*(10.0^(distance2/5.0))^2.0

plot, spec1_wave/1.00489, funred*4.0*3.14159265*(10.0^(distance1/5.0))^2.0, xrange=[2000,6000]
oplot, spec3_wave, spec3_flux*4.0*3.14159265*(10.0^(distance3/5.0))^2.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xdata=[5,6,7]
ydata=[4,5,6]

nxplots=1
nyplots=4
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
wall = 0.03*8.8 ; these values were originally made for the 8.8 width
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

figurename=spec1name+'_'+spec2name+'_specex.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=3
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='F!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,ytop], ystyle=1, xrange=xrange0, xstyle=1, $
;yrange=[0,ceil(max([spec1_smooth,spec2_smooth]+1.0)/10.0^max)], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

cgoplot,  smooth_wave[0:43], spec1_smooth[0:43]/10.0^(max[0]), color='red', line=0
cgoplot,  smooth_wave[0:43], spec2_smooth[0:43]/10.0/10.0^(max[0]), color=color2, line=2

al_legend, [spec2name+' / 10', spec1name],   linestyle=[2,0],  box=0, $
pos=[0.14,0.95], color=[cgcolor(color2),cgcolor('red')],/norm, charsize=1


;;;;;


nx=0
ny=2
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='R!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,20], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

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
oplot, lambda, 	lambda_av_salt/EBV, linestyle=4
;oplot, lambda, 	lambda_av_comb/EBV, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14/EBV, linestyle=0, thick=4

cgoplot, smooth_wave[0:43], ( (distance2-distance1) +2.5*alog10(spec2_smooth[0:43]/spec1_smooth[0:43]))/Ebv, thick=5, color='red'

al_legend, ['MW R_V=3.1', 'CCM R_V=1.7', 'SMC', 'G08LMC', 'SALT'],   linestyle=[0,1,2,3,4],  background_color='white', $
pos=[0.42,0.72], /norm, charsize=1, box=0

;;;;;;;;;;;

nx=0
ny=1
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='A!l!7k!3!N - A!DV!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,5], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

rv=3.1
	lambda_av_mw31=sne_mw_reddening(lambda,Ebv,rv=rv)
rv=1.7	
lambda_av_mw17=sne_mw_reddening(lambda,Ebv,rv=rv)
	lambda_av_smc=sne_smc_reddening(lambda,Ebv)
	lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv)
	lambda_av_salt=sne_salt_reddening(lambda,Ebv)

;oplot, lambda, 	lambda_av_mw31-3.1*EBV, linestyle=0
;oplot, lambda, 	lambda_av_mw17-1.7*EBV, linestyle=1
;oplot, lambda, 	lambda_av_smc-2.9*EBV, linestyle=2
;oplot, lambda, 	lambda_av_cslmc-1.7*EBV, linestyle=3
;oplot, lambda, 	lambda_av_salt-1.16*EBV, linestyle=4

av31=lambda_av_mw31[where(lambda eq 5500)]
av11=lambda_av_mw11[where(lambda eq 5500)]
av14=lambda_av_mw14[where(lambda eq 5500)]
av17=lambda_av_mw17[where(lambda eq 5500)]
av26=lambda_av_mw26[where(lambda eq 5500)]
avsmc=lambda_av_smc[where(lambda eq 5500)]
avcslmc=lambda_av_cslmc[where(lambda eq 5500)]
avsalt=lambda_av_salt[where(lambda eq 5500)]
avcomb=lambda_av_comb[where(lambda eq 5500)]

oplot, lambda, 	lambda_av_mw31-av31[0], linestyle=0
oplot, lambda, 	lambda_av_mw17-av17[0], linestyle=1
oplot, lambda, 	lambda_av_smc-avsmc[0], linestyle=2
oplot, lambda, 	lambda_av_cslmc-avcslmc[0], linestyle=3
oplot, lambda, 	lambda_av_salt-avsalt[0], linestyle=4
;oplot, lambda, 	lambda_av_comb-avcomb[0], linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14-av14[0], linestyle=0, thick=4

cgoplot, smooth_wave[0:43], (2.5*alog10(spec2_smooth[0:43]/spec1_smooth[0:43]))-(spec1_mags[5]-spec2_mags[5]), thick=5, color='red'

;;;;;;;;;;;

nx=0
ny=0
cgplot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=xtitle0, ytitle='A!l!7k!3!N', charsize=1.0,  $
 yrange=[0,5], ystyle=1, xrange=xrange0, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, $
xtickname=[' ','2000',' ','3000',' ','4000',' ','5000',' ',' ']

cgoplot, smooth_wave[0:43], (distance2-distance1) +2.5*alog10(spec2_smooth[0:43]/spec1_smooth[0:43]), thick=5, color='red'

oplot, lambda, 	lambda_av_mw31, linestyle=0
oplot, lambda, 	lambda_av_mw17, linestyle=1
oplot, lambda, 	lambda_av_smc, linestyle=2
oplot, lambda, 	lambda_av_cslmc, linestyle=3
oplot, lambda, 	lambda_av_salt, linestyle=4
;oplot, lambda, 	lambda_av_mw11, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_comb, linestyle=0, thick=4

device, /close
SET_PLOT, 'X'

;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

figurename=spec1name+'_'+spec3name+'_specex.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

;;;   x 10^'+strtrim(max,1)
nx=0
ny=3
;;;;;;;;;;;;;;; plotting the spectra

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='F!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,ytop], ystyle=1, xrange=xrange0, xstyle=1, $
;yrange=[0,ceil(max([spec1_smooth,spec2_smooth]+1.0)/10.0^max)], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

cgoplot,  smooth_wave[0:43], spec1_smooth[0:43]/10.0^(max[0]), color='red', line=0
cgoplot,  smooth_wave[0:43], spec3_smooth[0:43]/5.0/10.0^(max[0]), color=color2, line=2

al_legend, [spec3name+' / 5', spec1name],   linestyle=[2,0],  box=0, $
pos=[0.14,0.95], color=[cgcolor(color2),cgcolor('red')],/norm, charsize=1


;;;;;


nx=0
ny=2
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='R!l!7k!3!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,20], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

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
oplot, lambda, 	lambda_av_salt/EBV, linestyle=4
;oplot, lambda, 	lambda_av_comb/EBV, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14/EBV, linestyle=0, thick=4

cgoplot, smooth_wave[0:43], ( (distance3-distance1) +2.5*alog10(spec3_smooth[0:43]/spec1_smooth[0:43]))/Ebv, thick=5, color='red'

al_legend, ['MW R_V=3.1', 'CCM R_V=1.7', 'SMC', 'G08LMC', 'SALT'],   linestyle=[0,1,2,3,4],  background_color='white', $
pos=[0.42,0.72], /norm, charsize=1, box=0

;;;;;;;;;;;

nx=0
ny=1
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle='A!l!7k!3!N - A!DV!N', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0,5], ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

rv=3.1
	lambda_av_mw31=sne_mw_reddening(lambda,Ebv,rv=rv)
rv=1.7	
lambda_av_mw17=sne_mw_reddening(lambda,Ebv,rv=rv)
	lambda_av_smc=sne_smc_reddening(lambda,Ebv)
	lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv)
	lambda_av_salt=sne_salt_reddening(lambda,Ebv)

;oplot, lambda, 	lambda_av_mw31-3.1*EBV, linestyle=0
;oplot, lambda, 	lambda_av_mw17-1.7*EBV, linestyle=1
;oplot, lambda, 	lambda_av_smc-2.9*EBV, linestyle=2
;oplot, lambda, 	lambda_av_cslmc-1.7*EBV, linestyle=3
;oplot, lambda, 	lambda_av_salt-1.16*EBV, linestyle=4

av31=lambda_av_mw31[where(lambda eq 5500)]
av11=lambda_av_mw11[where(lambda eq 5500)]
av14=lambda_av_mw14[where(lambda eq 5500)]
av17=lambda_av_mw17[where(lambda eq 5500)]
av26=lambda_av_mw26[where(lambda eq 5500)]
avsmc=lambda_av_smc[where(lambda eq 5500)]
avcslmc=lambda_av_cslmc[where(lambda eq 5500)]
avsalt=lambda_av_salt[where(lambda eq 5500)]
avcomb=lambda_av_comb[where(lambda eq 5500)]

oplot, lambda, 	lambda_av_mw31-av31[0], linestyle=0
oplot, lambda, 	lambda_av_mw17-av17[0], linestyle=1
oplot, lambda, 	lambda_av_smc-avsmc[0], linestyle=2
oplot, lambda, 	lambda_av_cslmc-avcslmc[0], linestyle=3
oplot, lambda, 	lambda_av_salt-avsalt[0], linestyle=4
;oplot, lambda, 	lambda_av_comb-avcomb[0], linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14-av14[0], linestyle=0, thick=4

cgoplot, smooth_wave[0:43], (2.5*alog10(spec3_smooth[0:43]/spec1_smooth[0:43]))-(spec1_mags[5]-spec3_mags[5]), thick=5, color='red'

;;;;;;;;;;;

nx=0
ny=0
cgplot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=xtitle0, ytitle='A!l!7k!3!N', charsize=1.0,  $
 yrange=[0,5], ystyle=1, xrange=xrange0, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, $
xtickname=[' ','2000',' ','3000',' ','4000',' ','5000',' ',' ']

cgoplot, smooth_wave[0:43], (distance3-distance1) +2.5*alog10(spec3_smooth[0:43]/spec1_smooth[0:43]), thick=5, color='red'

oplot, lambda, 	lambda_av_mw31, linestyle=0
oplot, lambda, 	lambda_av_mw17, linestyle=1
oplot, lambda, 	lambda_av_smc, linestyle=2
oplot, lambda, 	lambda_av_cslmc, linestyle=3
oplot, lambda, 	lambda_av_salt, linestyle=4
;oplot, lambda, 	lambda_av_mw11, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_mw14, linestyle=0, thick=4
;oplot, lambda, 	lambda_av_comb, linestyle=0, thick=4

device, /close
SET_PLOT, 'X'


$open SN2015F_SN2011fe_specex.eps

$open SN2015F_SN2011by_specex.eps



stop
end
