;;;;;

pro specphotcounts, spectrum, specmags, speccounts, redfrac=redfrac

h = 1D*6.626E-27
c = 1D*2.998E18
hc = 1D*1.986E-8
e = 1D*2.71828
k = 1D*1.38E-16


;;Read in the filter Effective Area curves
readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent
readcol,"$SNSCRIPTS/B_UVOT.txt",lambda,B_EA,/silent
readcol,"$SNSCRIPTS/U_UVOT.txt",lambda,U_EA,/silent
readcol,"$SNSCRIPTS/UVW1_2010.txt",lambda,W1_EA,/silent
readcol,"$SNSCRIPTS/UVM2_2010.txt",lambda,M2_EA,/silent
readcol,"$SNSCRIPTS/UVW2_2010.txt",lambda,W2_EA,/silent
readcol,"$SNSCRIPTS/UVW1-rc.txt",lambda,W1rc_EA,/silent
readcol,"$SNSCRIPTS/UVW2-rc.txt",lambda,W2rc_EA,/silent

;filters=fltarr(9,n_elements(lambda)
;filters[0,*]=lambda


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
zero=where(finite(sp_flux) eq 0)
sp_flux[zero]=0.0

flux=interpol(sp_flux,sp_wave,lambda)


v_counts=V_EA*flux*(10*lambda/hc)
b_counts=B_EA*flux*(10*lambda/hc)
u_counts=U_EA*flux*(10*lambda/hc)
w1_counts=W1_EA*flux*(10*lambda/hc)
m2_counts=M2_EA*flux*(10*lambda/hc)
w2_counts=W2_EA*flux*(10*lambda/hc)
w1rc_counts=W1rc_EA*flux*(10*lambda/hc)
w2rc_counts=W2rc_EA*flux*(10*lambda/hc)
;window, 0
;plot, lambda, m2_counts
;oplot, lambda, w2rc_counts

v_spec=-2.5*alog10(total(v_counts))+17.89
b_spec=-2.5*alog10(total(b_counts))+19.11
u_spec=-2.5*alog10(total(u_counts))+18.34
w1_spec=-2.5*alog10(total(w1_counts))+17.44
m2_spec=-2.5*alog10(total(m2_counts))+16.85
w2_spec=-2.5*alog10(total(w2_counts))+17.38
w1rc_spec=-2.5*alog10(total(w1rc_counts))+17.34
w2rc_spec=-2.5*alog10(total(w2rc_counts))+17.25

redfrac=[total(w1rc_counts)/total(w1_counts),total(w2rc_counts)/total(w2_counts)]
;specmags=[v_spec, b_spec, u_spec, w1_spec, m2_spec, w2_spec, w1rc_spec, w2rc_spec]

speccounts=[total(w2_counts), total(m2_counts), total(w1_counts), total(u_counts), total(b_counts), total(v_counts), total(w1rc_counts), total(w2rc_counts)]
specmags=[w2_spec, m2_spec, w1_spec, u_spec, b_spec, v_spec, w1rc_spec, w2rc_spec]
;print, "specmags", specmags


end


;;;;;;;;;;
pro addsedendpoints, inputmags, sedwavelengths, sedfluxes, outputsedwavelengths, outputsedfluxes

bluewavelength=1000
redwavelength=9000

;zeropoints=[17.35, 16.82, 17.49, 18.34, 19.11, 17.89]
;inputmagcounts=10.0^((zeropoints-inputmags)/2.5)

readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent

;print, 'sed fluxes ', sedfluxes
lowred=10.0^(-30)
highred=10.0*sedfluxes[n_elements(sedfluxes)-1]
; print, 'orig high,low red ', lowred, highred
nsteps=10.0
diffv=fltarr(nsteps)

for i=0,3 do begin

	for f=0.0,nsteps-1 do begin
		rededge=lowred+f/nsteps*(highred-lowred)
		newfluxspectrum=interpol([sedfluxes,rededge],[sedwavelengths,redwavelength], lambda)
		spectrum=[transpose(lambda),transpose(newfluxspectrum)]
		specphotcounts, spectrum, specmags, speccounts
		diffv[f]=(specmags[5]-inputmags[5])
	endfor

	bestred=where(abs(diffv) eq min(abs(diffv)))
	rededge=lowred+float(bestred[0])/nsteps*(highred-lowred)

	lowrednew=lowred+(bestred[0]-1)/nsteps*(highred-lowred) 
	highrednew=lowred+(bestred[0]+1)/nsteps*(highred-lowred)
 
;	print, "lowred, highred", lowrednew, highrednew
	lowred=max([10.0^(-30),lowrednew])
	highred=highrednew

;print, "red ", bestred, rededge, min(abs(diffv))

endfor


diffw2=fltarr(nsteps)
diffm2=fltarr(nsteps)
diffw2m2=fltarr(nsteps)

lowblue=10.0^(-30)
highblue=10.0*sedfluxes[0]
 
for i=0,3 do begin
	for f=0,nsteps-1 do begin
		blueedge=lowblue+f/nsteps*(highblue-lowblue)
		newfluxspectrum=interpol([blueedge,sedfluxes,rededge],[bluewavelength,sedwavelengths,redwavelength], lambda)
		spectrum=[transpose(lambda),transpose(newfluxspectrum)]
		specphotcounts, spectrum, specmags, speccounts
		diffw2[f]=(specmags[0]-inputmags[0])
		diffm2[f]=(specmags[1]-inputmags[1])
		diffw2m2[f]=abs(specmags[0]-inputmags[0]) + abs(specmags[1]-inputmags[1])
		diffv[f]=(specmags[5]-inputmags[5])
	endfor

	bestblue=where(abs(diffw2m2) eq min(abs(diffw2m2)))
	blueedge=lowblue+float(bestblue[0])/nsteps*(highblue-lowblue)

	lowbluenew=lowblue+(bestblue[0]-1)/nsteps*(highblue-lowblue) 
	highbluenew=lowblue+(bestblue[0]+1)/nsteps*(highblue-lowblue)
 
	lowblue=max([10.0^(-30),lowbluenew])
	highblue=highbluenew
endfor


outputsedwavelengths=[bluewavelength,sedwavelengths,redwavelength]
outputsedfluxes=[blueedge,sedfluxes,rededge]

newfluxspectrum=interpol([blueedge,sedfluxes,rededge],[bluewavelength,sedwavelengths,redwavelength],lambda)

spectrum=[transpose(lambda),transpose(newfluxspectrum)]

end

pro hst_uvot_plot, SNname, dt, mjdepochs, uvotepochmags, uvotepochmagerrs, hstepochmags

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  DOUBLE PLOT two/one third ;;;;;;;;;;;;;;;;;;;;
;

nplots=2
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.03
margin=0.12
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
xrange=[55680.0,55700.0]
xrange=[55675.0-55500.0,55710.0-55500.0]
yrange2=[max(dt.mag_array,/nan)+1.0,min(dt.mag_array,/nan)-1.0]
yrange2=[20.0,12.0]
yrange1=[-0.1,0.3]
xtitle1='Modified Julian Date-55500.0'
ytitle1='UVOT-HST'
ytitle2='Vega Mags'
nxticks=9
nyticks=9
nyticks2=5
figurename=SNname+'_res.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=10, bits_per_pixel=8, /color

plot, dt.vv[0,*]-55500.0, dt.vv[1,*], /nodata, /noerase, position=[x1,y1+(y2-y1)/2,x2,y1+(y2-y1)/2+y2-y1], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange2, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

oploterror, dt.w2[0,*]-55500.0, dt.w2[1,*], dt.w2[2,*], psym=cgsymcat(5), color=cgcolor('black')
oploterror, dt.m2[0,*]-55500.0, dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(1), color=cgcolor('red')
oploterror, dt.w1[0,*]-55500.0, dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(2), color=cgcolor('maroon')
oploterror, dt.uu[0,*]-55500.0, dt.uu[1,*], dt.uu[2,*], psym=cgsymcat(4), color=cgcolor('purple')
oploterror, dt.bb[0,*]-55500.0, dt.bb[1,*], dt.bb[2,*], psym=cgsymcat(6), color=cgcolor('blue')
oploterror, dt.vv[0,*]-55500.0, dt.vv[1,*], dt.vv[2,*], psym=cgsymcat(9), color=cgcolor('dark green')


oplot, mjdepochs-55500.0, hstepochmags[0,*], symsize=2, psym=cgsymcat(5), color=cgcolor('black')
oplot, mjdepochs-55500.0, hstepochmags[1,*], symsize=2, psym=cgsymcat(1), color=cgcolor('red')
oplot, mjdepochs-55500.0, hstepochmags[2,*], symsize=2, psym=cgsymcat(2), color=cgcolor('maroon')
oplot, mjdepochs-55500.0, hstepochmags[3,*], symsize=2, psym=cgsymcat(4), color=cgcolor('purple')
oplot, mjdepochs-55500.0, hstepochmags[4,*], symsize=2, psym=cgsymcat(6), color=cgcolor('blue')
oplot, mjdepochs-55500.0, hstepochmags[5,*], symsize=2, psym=cgsymcat(9), color=cgcolor('dark green')




plot, mjdepochs-55500.0, uvotepochmags[0,*]-hstepochmags[0,*], /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle1, ytitle=ytitle1, charsize=1.0,  $
 yrange=[-0.2,0.3], ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks/2+1.0, ytickv=ytickvalues

oploterror, mjdepochs-55500.0+0.2, uvotepochmags[0,*]-hstepochmags[0,*], uvotepochmagerrs[0,*], psym=cgsymcat(5), color=cgcolor('black')
oploterror, mjdepochs-55500.0+0.4, uvotepochmags[1,*]-hstepochmags[1,*], uvotepochmagerrs[1,*], psym=cgsymcat(1), color=cgcolor('red')
oploterror, mjdepochs-55500.0-1.2, uvotepochmags[2,*]-hstepochmags[2,*], uvotepochmagerrs[2,*], psym=cgsymcat(2), color=cgcolor('maroon')
oploterror, mjdepochs-55500.0-0.8, uvotepochmags[3,*]-hstepochmags[3,*], uvotepochmagerrs[3,*], psym=cgsymcat(4), color=cgcolor('purple')
oploterror, mjdepochs-55500.0-0.4, uvotepochmags[4,*]-hstepochmags[4,*], uvotepochmagerrs[4,*], psym=cgsymcat(6), color=cgcolor('blue')
oploterror, mjdepochs-55500.0, uvotepochmags[5,*]-hstepochmags[5,*], uvotepochmagerrs[5,*], psym=cgsymcat(9), color=cgcolor('dark green')

al_legend, ['w2','m2','w1','u','b','v'], color=[cgcolor('black'), cgcolor('red'),cgcolor('maroon') ,cgcolor('purple') , cgcolor('blue'), cgcolor('dark green') ], psym=[5,1,2,4,6,9], pos=[202.0, 0.3], box=0

device, /close
SET_PLOT, 'X'
$open SN2011by_res.eps 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end


pro hst_uvot

;;;; read in uvot data
SNname='SN2011by'
snphot_array, '$SNFOLDER/SwiftSNarchive/PhotArchive/'+SNname+'_uvotB13.dat', '$SNFOLDER/UVOTDATA/'+SNname+'_uvotB13.fits', dt=dt
;snphot_array, '$SNFOLDER/SwiftSNarchive/'+SNname+'/'+SNname+'_uvotB13.dat', '$SNFOLDER/UVOTDATA/'+SNname+'_uvotB13.fits', dt=dt

dates=['2011-04-30 11:33:41','2011-05-05 03:00:06','2011-05-09 06:15:47','2011-05-13 17:39:31','2011-05-18 04:32:31']
spectra=["SN2011by/SN2011by_0430_low.dat","SN2011by/SN2011by_0505_opt.dat","SN2011by/SN2011by_0509_full.dat","SN2011by/SN2011by_0513_opt.dat","SN2011by/SN2011by_0518_opt.dat"]
mjdepochs=dblarr(n_elements(dates))
uvotepochmags=fltarr(6,n_elements(dates))
uvotepochmagerrs=fltarr(6,n_elements(dates))
hstepochmags=fltarr(6,n_elements(dates))
for n=0,n_elements(dates)-1 do mjdepochs[n]=date_conv(dates[n],'MODIFIED')

for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmags[f,n]=interpol(dt.mag_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])
for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmagerrs[f,n]=interpol(dt.magerr_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])

wave_short=[1810.00 , 2080.00 , 2320.00 , 3240.00 , 4070.00 , 5210.00    ]
wave_long=[ 2270.00 , 2420.00 , 2850.00 , 3700.00 , 4630.00 , 5650.00    ]

for n=0,n_elements(dates)-1 do begin
	spectrum=spectra[n]
	readcol,spectrum,sedwavelengths,sedfluxes,/silent
	addsedendpoints, uvotepochmags[*,n], sedwavelengths, sedfluxes, outputsedwavelengths, outputsedfluxes
	uvot_specphot, [transpose(outputsedwavelengths),transpose(outputsedfluxes)], specmags, speccounts  

for f=0,5 do if sedwavelengths[0] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-1] lt wave_long[f] then speccounts[f]= !Values.F_NAN
for f=0,5 do if sedwavelengths[0] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-1] lt wave_long[f] then specmags[f]= !Values.F_NAN
	hstepochmags[*,n]=specmags[0:5]
endfor

hst_uvot_plot, SNname, dt, mjdepochs, uvotepochmags, uvotepochmagerrs, hstepochmags


print, 'final stop'
stop
end

