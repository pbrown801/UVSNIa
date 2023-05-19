pro SN2011fe_uvot_hst_comparison

;;;; read in uvot data
SNname='SN2011fe'

Bpeaktime=55814.51

pjb_phot_array_B141, '$SNFOLDER/SwiftSNarchive/SN2011fe/'+SNname+'_uvotB14.1.dat', dt=dt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SN2011fe sullivan files
readcol,'SN2011fe_sullivanfiles.txt', sullivanspectrafiles, format='A', /silent
sullivandates=[ '2011-08-28 12:00:00', '2011-08-31 06:57:46','2011-09-03 10:00:01','2011-09-07 09:58:46','2011-09-10 10:02:09','2011-09-13 16:04:10','2011-09-19 14:09:10','2011-10-01 06:46:19','2011-10-07 08:24:35','2011-10-21 05:02:26' ]
sullivanmjd=fltarr(n_elements(sullivandates))
sullivanmags=fltarr(6,n_elements(sullivandates))

for n=0,n_elements(sullivandates)-1 do sullivanmjd[n] = date_conv(sullivandates[n],'MODIFIED')

for n=0,n_elements(sullivandates)-1 do begin
	spectrum=sullivanspectrafiles[n]
	readcol,spectrum,sedwavelengths,sedfluxes,/silent

	uvot_specphot, [transpose(sedwavelengths),transpose(sedfluxes)], specmags, speccounts=speccounts

	sullivanmags[*,n]=specmags[0:5]
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

bolofiles=['SN2011fe_bolo/M001522.dat','SN2011fe_bolo/M001424.dat','SN2011fe_bolo/M001324.dat','SN2011fe_bolo/M001223.dat', 'SN2011fe_bolo/M001124.dat','SN2011fe_bolo/M001025.dat','SN2011fe_bolo/M000924.dat','SN2011fe_bolo/M000825.dat', 'SN2011fe_bolo/M000723.dat','SN2011fe_bolo/M000625.dat','SN2011fe_bolo/M000526.dat','SN2011fe_bolo/M000125.dat', 'SN2011fe_bolo/M000027.dat','SN2011fe_bolo/P000072.dat','SN2011fe_bolo/P000172.dat','SN2011fe_bolo/P000271.dat', 'SN2011fe_bolo/P000372.dat', 'SN2011fe_bolo/P000671.dat','SN2011fe_bolo/P000871.dat','SN2011fe_bolo/P001170.dat', 'SN2011fe_bolo/P001370.dat','SN2011fe_bolo/P001669.dat','SN2011fe_bolo/P001869.dat','SN2011fe_bolo/P002169.dat','SN2011fe_bolo/P002368.dat']
boloepoch=[-15.22,-14.24,-13.24,-12.23,-11.24,-10.25,-9.24,-8.25,-7.23,-6.25,-5.26,-1.25,-0.27,0.72,1.72,2.71,3.72,6.71,8.71,11.70,13.70,16.69,18.69,21.69,23.68]
bolomjd=boloepoch+Bpeaktime

nbolospectra=n_elements(bolofiles)

bolomags=fltarr(6,nbolospectra)

for n=0,nbolospectra-1 do begin

	spectrum=bolofiles[n]
	readcol,spectrum,sedwavelengths,sedfluxes,/silent

	Ebvmw=0.0088
	;; reddening the spectrum to get it as observed
	Amw=sne_mw_reddening(sedwavelengths,-Ebvmw)
	
	specred=sedfluxes*10^(Amw/2.5)

	uvot_specphot, [transpose(sedwavelengths),transpose(specred)], specmags, speccounts=speccounts

	bolomags[*,n]=specmags[0:5]
;	print, mag_array[*,n]
endfor





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


dates=[ '2011-08-31 06:57:46','2011-09-03 10:00:01','2011-09-07 09:58:46','2011-09-10 10:02:09','2011-09-13 16:04:10','2011-09-19 14:09:10','2011-09-19 15:27:39','2011-10-01 06:46:19','2011-10-07 08:24:35','2011-10-21 05:02:26' ]
mypath='$SNFOLDER/UVspectra/hst/SN2011fe/myspec/'
spectra= mypath+['SN2011fe_110831_full.dat', 'SN2011fe_110903_full.dat', 'SN2011fe_110907_full.dat', 'SN2011fe_110910_full.dat', 'SN2011fe_110913_full.dat', 'SN2011fe_110919_ccd.dat', 'SN2011fe_110919_mama.dat', 'SN2011fe_111001_full.dat', 'SN2011fe_111007_full.dat', 'SN2011fe_111021_full.dat' ]
myspecdet=[ 'mama', 'ccd', 'ccd', 'ccd', 'ccd', 'ccd', 'mama', 'mama' , 'mama', 'mama' ]

ccd=where(myspecdet eq 'ccd')
mama=where(myspecdet eq 'mama')


mjdepochs=dblarr(n_elements(dates))
uvotepochmags=fltarr(6,n_elements(dates))
uvotepochmagerrs=fltarr(6,n_elements(dates))
hstepochmags=fltarr(6,n_elements(dates))
for n=0,n_elements(dates)-1 do mjdepochs[n]=date_conv(dates[n],'MODIFIED')

for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmags[f,n]=interpol(dt.mag_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])

uvotepochmags[5,1]=!Values.F_NAN
uvotepochmags[5,2]=!Values.F_NAN
uvotepochmags[5,3]=!Values.F_NAN
uvotepochmags[5,4]=!Values.F_NAN
uvotepochmags[5,5]=!Values.F_NAN
uvotepochmags[4,0]=!Values.F_NAN
uvotepochmags[4,1]=!Values.F_NAN
uvotepochmags[4,2]=!Values.F_NAN
uvotepochmags[4,3]=!Values.F_NAN
uvotepochmags[4,4]=!Values.F_NAN
uvotepochmags[4,5]=!Values.F_NAN
uvotepochmags[4,6]=!Values.F_NAN
uvotepochmags[4,7]=!Values.F_NAN
uvotepochmags[4,8]=!Values.F_NAN
uvotepochmags[3,0]=!Values.F_NAN
uvotepochmags[3,1]=!Values.F_NAN
uvotepochmags[3,2]=!Values.F_NAN
uvotepochmags[3,3]=!Values.F_NAN
uvotepochmags[3,4]=!Values.F_NAN
uvotepochmags[3,5]=!Values.F_NAN
uvotepochmags[3,6]=!Values.F_NAN

for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmagerrs[f,n]=interpol(dt.magerr_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])

wave_short=[1810.00 , 2080.00 , 2320.00 , 3240.00 , 4070.00 , 5210.00    ]
wave_long=[ 2270.00 , 2420.00 , 2850.00 , 3700.00 , 4630.00 , 5650.00    ]

for n=0,n_elements(dates)-1 do begin
	spectrum=spectra[n]
	readcol,spectrum,sedwavelengths,sedfluxes,/silent
;	addsedendpoints, uvotepochmags[*,n], sedwavelengths, sedfluxes, outputsedwavelengths, outputsedfluxes
;	uvot_specphot, [transpose(outputsedwavelengths),transpose(outputsedfluxes)], specmags, speccounts  

if sedwavelengths[0] gt 1600 then sedfluxes=[sedfluxes[0],sedfluxes]
if sedwavelengths[n_elements(sedwavelengths)-1] lt 8000 then sedfluxes=[sedfluxes,sedfluxes[n_elements(sedfluxes)-1]]
if sedwavelengths[0] gt 1600 then sedwavelengths=[1500,sedwavelengths]
if sedwavelengths[n_elements(sedwavelengths)-1] lt 8000 then sedwavelengths=[sedwavelengths,9000]
	uvot_specphot, [transpose(sedwavelengths),transpose(sedfluxes)], specmags, speccounts=speccounts



for f=0,5 do if sedwavelengths[1] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-2] lt wave_long[f] then speccounts[f]= !Values.F_NAN
for f=0,5 do if sedwavelengths[1] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-2] lt wave_long[f] then specmags[f]= !Values.F_NAN
	hstepochmags[*,n]=specmags[0:5]
endfor


maghigh=floor(max(dt.mag_array,/nan))+1.0
maglow=floor(min(dt.mag_array,/nan))-1.0
mjdhigh=floor(max(dt.time_array,/nan)/10.0)*10.0
mjdlow=floor(min(dt.time_array,/nan)/10.0)*10.0

orchid  = FSC_COLOR('Orchid',!D.TABLE_SIZE-2)
blue  = FSC_COLOR('Blue',!D.TABLE_SIZE-3)
cyan  = FSC_COLOR('Cyan',!D.TABLE_SIZE-4)
green  = FSC_COLOR('Green',!D.TABLE_SIZE-5)
orange  = FSC_COLOR('Orange',!D.TABLE_SIZE-6)
red  = FSC_COLOR('Red',!D.TABLE_SIZE-7)
black  = FSC_COLOR('Black',!D.TABLE_SIZE-8)
white  = FSC_COLOR('White',!D.TABLE_SIZE-8)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  DOUBLE PLOT two/one third ;;;;;;;;;;;;;;;;;;;;
;

nplots=1.5
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
margin=0.16
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
xrange=[mjdlow-Bpeaktime[0],mjdhigh-Bpeaktime[0]]
xrange=[-15,45]
yrange2=[max(dt.mag_array,/nan)+1.0,min(dt.mag_array,/nan)-1.0]
yrange2=[17.0,9.0]
yrange1=[-0.2,0.3]
xtitle1='Days from Maximum Light'

ytitle1='m2 Residuals'
ytitle2='Vega Mags'
nxticks=4
nyticks=4
nyticks2=4
figurename='SN2011fe_res.eps'

fontsize=14
SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, dt.vv[0,*]-Bpeaktime[0], dt.vv[1,*], /nodata, /noerase, position=[x1,y1+(y2-y1)/2,x2,y1+(y2-y1)/2+y2-y1], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange2, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

oploterror, dt.w2[0,*]-Bpeaktime[0], dt.w2[1,*], dt.w2[2,*], psym=cgsymcat(5), color=cgcolor('black'), errcolor=cgcolor('black')
oploterror, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(1), color=cgcolor('red'), errcolor=cgcolor('red')
oploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(2), color=cgcolor('maroon'), errcolor=cgcolor('maroon')
oploterror, dt.uu[0,*]-Bpeaktime[0], dt.uu[1,*], dt.uu[2,*], psym=cgsymcat(4), color=cgcolor('purple'), errcolor=cgcolor('purple')
oploterror, dt.bb[0,*]-Bpeaktime[0], dt.bb[1,*], dt.bb[2,*], psym=cgsymcat(6), color=cgcolor('blue'), errcolor=cgcolor('blue')
oploterror, dt.vv[0,*]-Bpeaktime[0], dt.vv[1,*], dt.vv[2,*], psym=cgsymcat(9), color=cgcolor('green'), errcolor=cgcolor('green')

;oplot, mjdepochs-Bpeaktime[0], hstepochmags[0,*], symsize=2, psym=cgsymcat(5), color=cgcolor('black')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*], symsize=2, psym=cgsymcat(1), color=cgcolor('red')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*], symsize=2, psym=cgsymcat(2), color=cgcolor('maroon')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[3,*], symsize=2, psym=cgsymcat(4), color=cgcolor('purple')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[4,*], symsize=2, psym=cgsymcat(6), color=cgcolor('blue')
;oplot, mjdepochs-Bpeaktime[0], hstepochmags[5,*], symsize=2, psym=cgsymcat(9), color=cgcolor('green')



templatefiles=['$SNSCRIPTS/SNIa_w2_template_fe.dat', $
'$SNSCRIPTS/SNIa_m2_template_fe.dat', $
'$SNSCRIPTS/SNIa_w1_template_fe.dat', $
'$SNSCRIPTS/SNIa_w1_template_fe.dat', $
'$SNFOLDER/mlcs/vectors_B.dat', '$SNFOLDER/mlcs/vectors_V.dat' ]

;; account for differences in template definitions
templatemagshift=[0.0,0.0,0.0,0.0,19.570837, 19.511642]
templateepochshift=[0.0,0.0,0.0,0.0,0.0, -1.0]


f=1
	readcol, templatefiles[f], templateepoch, templatemags, /silent
	templateepoch=templateepoch+templateepochshift[f]
	templatemags=templatemags+templatemagshift[f]

interpedtemplatemags=fltarr(n_elements(dt.m2time))
for n=0,n_elements(interpedtemplatemags)-1 do interpedtemplatemags[n]=interpol(templatemags+13.06-0.1,templateepoch+55814.3-0.3,dt.m2time[n])

;plot, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*], xrange=[-15,0], yrange=[18,12]
;oplot, dt.m2[0,*]-Bpeaktime[0],interpedtemplatemags

;;;;;;;;;;
ploterror, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*]-interpedtemplatemags, dt.m2[2,*], xrange=[-15,45], yrange=[-0.2,0.2], /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle1, ytitle=ytitle1, charsize=1.0,  $
ystyle=1,xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



hstinterpedtemplatemags=fltarr(n_elements(dates))
for n=0,n_elements(hstinterpedtemplatemags)-1 do hstinterpedtemplatemags[n]=interpol(templatemags+13.06-0.1,templateepoch+55814.3-0.3,mjdepochs[n])


oplot, mjdepochs[ccd]-Bpeaktime[0], hstepochmags[1,ccd]-hstinterpedtemplatemags[ccd], symsize=1, psym=cgsymcat(9), color=cgcolor('red')
oplot, mjdepochs[mama]-Bpeaktime[0], hstepochmags[1,mama]-hstinterpedtemplatemags[mama], symsize=1, psym=cgsymcat(16), color=cgcolor('red')

sullivaninterpedtemplatemags=fltarr(n_elements(sullivandates))
for n=0,n_elements(sullivaninterpedtemplatemags)-1 do sullivaninterpedtemplatemags[n]=interpol(templatemags+13.06-0.1,templateepoch+55814.3-0.3,sullivanmjd[n])


oplot, sullivanmjd-Bpeaktime[0], sullivanmags[1,*]-sullivaninterpedtemplatemags, symsize=1, psym=5, color=cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;
bolointerpedtemplatemags=fltarr(n_elements(bolomjd))
for n=0,n_elements(bolointerpedtemplatemags)-1 do bolointerpedtemplatemags[n]=interpol(templatemags+13.06-0.1,templateepoch+55814.3-0.3,bolomjd[n])


oplot, bolomjd-Bpeaktime[0], bolomags[1,*]-bolointerpedtemplatemags, symsize=1, psym=5, color=cgcolor('green')

;;;;;;;;;;;;;;;;;;;;;

; .run SN2011fe_uvot_hst_comparison
; SN2011fe_uvot_hst_comparison


;al_legend, ['w2','m2','w1','u','b','v'], color=[cgcolor('black'), ;cgcolor('red'),cgcolor('maroon') ,cgcolor('purple') , cgcolor('blue'), cgcolor('green') ], psym=[5,1, 2,4,6, 9 ], position=[845.0,0.3], box=0, charsize=0.7

device, /close
SET_PLOT, 'X'
$open SN2011fe_res.eps 



xrange=[mjdlow-Bpeaktime[0],mjdhigh-Bpeaktime[0]]
xrange=[-15,45]
yrange2=[max(dt.mag_array,/nan)+1.0,min(dt.mag_array,/nan)-1.0]
yrange2=[17.0,9.0]
yrange1=[-0.2,0.3]
xtitle1='Days from Maximum Light'

ytitle1='w1  Residuals'
ytitle2='Vega Mags'
nxticks=4
nyticks=4
nyticks2=4
figurename='SN2011fe_w1res.eps'

fontsize=14
SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, dt.vv[0,*]-Bpeaktime[0], dt.vv[1,*], /nodata, /noerase, position=[x1,y1+(y2-y1)/2,x2,y1+(y2-y1)/2+y2-y1], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange2, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

oploterror, dt.w2[0,*]-Bpeaktime[0], dt.w2[1,*], dt.w2[2,*], psym=cgsymcat(5), color=cgcolor('black'), errcolor=cgcolor('black')
oploterror, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(1), color=cgcolor('red'), errcolor=cgcolor('red')
oploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(2), color=cgcolor('maroon'), errcolor=cgcolor('maroon')
oploterror, dt.uu[0,*]-Bpeaktime[0], dt.uu[1,*], dt.uu[2,*], psym=cgsymcat(4), color=cgcolor('purple'), errcolor=cgcolor('purple')
oploterror, dt.bb[0,*]-Bpeaktime[0], dt.bb[1,*], dt.bb[2,*], psym=cgsymcat(6), color=cgcolor('blue'), errcolor=cgcolor('blue')
oploterror, dt.vv[0,*]-Bpeaktime[0], dt.vv[1,*], dt.vv[2,*], psym=cgsymcat(9), color=cgcolor('green'), errcolor=cgcolor('green')

;oplot, mjdepochs-Bpeaktime[0], hstepochmags[0,*], symsize=2, psym=cgsymcat(5), color=cgcolor('black')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*], symsize=2, psym=cgsymcat(1), color=cgcolor('red')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*], symsize=2, psym=cgsymcat(2), color=cgcolor('maroon')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[3,*], symsize=2, psym=cgsymcat(4), color=cgcolor('purple')
oplot, mjdepochs-Bpeaktime[0], hstepochmags[4,*], symsize=2, psym=cgsymcat(6), color=cgcolor('blue')
;oplot, mjdepochs-Bpeaktime[0], hstepochmags[5,*], symsize=2, psym=cgsymcat(9), color=cgcolor('green')



templatefiles=['$SNSCRIPTS/SNIa_w2_template_fe.dat', $
'$SNSCRIPTS/SNIa_m2_template_fe.dat', $
'$SNSCRIPTS/SNIa_w1_template_fe.dat', $
'$SNSCRIPTS/SNIa_w1_template_fe.dat', $
'$SNFOLDER/mlcs/vectors_B.dat', '$SNFOLDER/mlcs/vectors_V.dat' ]

;; account for differences in template definitions
templatemagshift=[0.0,0.0,0.0,0.0,19.570837, 19.511642]
templateepochshift=[0.0,0.0,0.0,0.0,0.0, -1.0]


f=2
	readcol, templatefiles[f], templateepoch, templatemags, /silent
	w1templateepoch=templateepoch+templateepochshift[f]
	w1templatemags=templatemags+templatemagshift[f]

w1interpedtemplatemags=fltarr(n_elements(dt.w1time))
for n=0,n_elements(w1interpedtemplatemags)-1 do w1interpedtemplatemags[n]=interpol(w1templatemags+11.02-0.08,w1templateepoch+55812.2-0.3,dt.w1time[n])

;ploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], xrange=[-5,5], yrange=[12,10.8]
;oplot, dt.w1[0,*]-Bpeaktime[0],w1interpedtemplatemags

;;;;;;;;;;
ploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*]-w1interpedtemplatemags, dt.w1[2,*], xrange=[-15,45], yrange=[-0.2,0.2], /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle1, ytitle=ytitle1, charsize=1.0,  $
ystyle=1,xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



w1hstinterpedtemplatemags=fltarr(n_elements(dates))
for n=0,n_elements(w1hstinterpedtemplatemags)-1 do w1hstinterpedtemplatemags[n]=interpol(w1templatemags+11.02-0.1,w1templateepoch+55812.2-0.0,mjdepochs[n])


oplot, mjdepochs[ccd]-Bpeaktime[0], hstepochmags[2,ccd]-w1hstinterpedtemplatemags[ccd], symsize=1, psym=cgsymcat(9), color=cgcolor('red')
oplot, mjdepochs[mama]-Bpeaktime[0], hstepochmags[2,mama]-w1hstinterpedtemplatemags[mama], symsize=1, psym=cgsymcat(16), color=cgcolor('red')

w1sullivaninterpedtemplatemags=fltarr(n_elements(sullivandates))
for n=0,n_elements(w1sullivaninterpedtemplatemags)-1 do w1sullivaninterpedtemplatemags[n]=interpol(w1templatemags+11.02-0.08,w1templateepoch+55812.2-0.3,sullivanmjd[n])


oplot, sullivanmjd-Bpeaktime[0], sullivanmags[2,*]-w1sullivaninterpedtemplatemags, symsize=1, psym=5, color=cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;
w1bolointerpedtemplatemags=fltarr(n_elements(bolomjd))
for n=0,n_elements(w1bolointerpedtemplatemags)-1 do w1bolointerpedtemplatemags[n]=interpol(w1templatemags+11.02-0.08,w1templateepoch+55812.2-0.3,bolomjd[n])


oplot, bolomjd-Bpeaktime[0], bolomags[2,*]-w1bolointerpedtemplatemags, symsize=1, psym=5, color=cgcolor('green')

;;;;;;;;;;;;;;;;;;;;;

; .run SN2011fe_uvot_hst_comparison
; SN2011fe_uvot_hst_comparison


;al_legend, ['w2','m2','w1','u','b','v'], color=[cgcolor('black'), ;cgcolor('red'),cgcolor('maroon') ,cgcolor('purple') , cgcolor('blue'), cgcolor('green') ], psym=[5,1, 2,4,6, 9 ], position=[845.0,0.3], box=0, charsize=0.7

device, /close
SET_PLOT, 'X'
$open SN2011fe_w1res.eps 



print, 'final stop'

stop

;;;;;;;;;;;;;;;;;;;;;


;oploterror, mjdepochs[where(myspecdet eq 'mama')]-Bpeaktime[0]-0.6, uvotepochmags[0,where(myspecdet eq 'mama')]-hstepochmags[0,where(myspecdet eq 'mama')], uvotepochmagerrs[0,where(myspecdet eq 'mama')],psym=cgsymcat(1), color=cgcolor('black'), errcolor=cgcolor('black'), symsize=1.5
oploterror, mjdepochs[where(myspecdet eq 'mama')]-Bpeaktime[0]+0.2, uvotepochmags[1,where(myspecdet eq 'mama')]-hstepochmags[1,where(myspecdet eq 'mama')], uvotepochmagerrs[1,where(myspecdet eq 'mama')],psym=cgsymcat(17,thick=3), color=cgcolor('red'), errcolor=cgcolor('red'), thick=2, symsize=1.5
oploterror, mjdepochs[where(myspecdet eq 'mama')]-Bpeaktime[0]-0.2, uvotepochmags[2,where(myspecdet eq 'mama')]-hstepochmags[2,where(myspecdet eq 'mama')], uvotepochmagerrs[2,where(myspecdet eq 'mama')],psym=cgsymcat(2,thick=3), color=cgcolor('maroon'), errcolor=cgcolor('maroon'), thick=2, symsize=1.5


;oploterror, mjdepochs[where(myspecdet eq 'ccd')]-Bpeaktime[0]-0.6, uvotepochmags[0,where(myspecdet eq 'ccd')]-hstepochmags[0,where(myspecdet eq 'ccd')], uvotepochmagerrs[0,where(myspecdet eq 'ccd')],psym=cgsymcat(1), color=cgcolor('black'), errcolor=cgcolor('black')
oploterror, mjdepochs[where(myspecdet eq 'ccd')]-Bpeaktime[0]+0.2, uvotepochmags[1,where(myspecdet eq 'ccd')]-hstepochmags[1,where(myspecdet eq 'ccd')], uvotepochmagerrs[1,where(myspecdet eq 'ccd')],psym=cgsymcat(17), color=cgcolor('red'), errcolor=cgcolor('red')
oploterror, mjdepochs[where(myspecdet eq 'ccd')]-Bpeaktime[0]-0.2, uvotepochmags[2,where(myspecdet eq 'ccd')]-hstepochmags[2,where(myspecdet eq 'ccd')], uvotepochmagerrs[2,where(myspecdet eq 'ccd')],psym=cgsymcat(2), color=cgcolor('maroon'), errcolor=cgcolor('maroon')

oplot, [mjdepochs[0]-Bpeaktime[0]-5.0, max(mjdepochs)-Bpeaktime[0]-5.0], [0.0, 0.0]

oploterror, mjdepochs-Bpeaktime[0]+0.4, uvotepochmags[3,*]-hstepochmags[3,*], uvotepochmagerrs[3,*],psym=cgsymcat(4), color=cgcolor('purple'), errcolor=cgcolor('purple')
oploterror, mjdepochs-Bpeaktime[0], uvotepochmags[4,*]-hstepochmags[4,*], uvotepochmagerrs[4,*],psym=cgsymcat(6), color=cgcolor('blue'), errcolor=cgcolor('blue')
;oploterror, mjdepochs-Bpeaktime[0]-0.4, uvotepochmags[5,*]-hstepochmags[5,*], uvotepochmagerrs[5,*],psym=cgsymcat(9), color=cgcolor('green'), errcolor=cgcolor('green')







plot, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*]-interpedtemplatemags, /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle1, ytitle=ytitle1, charsize=1.0,  $
 yrange=yrange1, ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



ploterror, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*], dt.m2[2,*], xrange=[800,860], yrange=[18,12.5]
oplot, dt.m2[0,*]-Bpeaktime[0], interpedtemplatemags                        

oplot, sullivanmjd-Bpeaktime[0], sullivanmags[1,*], psym=5, symsize=4

stop
end
