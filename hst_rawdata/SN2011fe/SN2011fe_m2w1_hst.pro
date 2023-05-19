pro SN2011fe_m2w1_hst

;;;; read in uvot data
SNname='SN2011fe'

Bpeaktime=55814.51

pjb_phot_array_B141, '$SNFOLDER/SwiftSNarchive/SN2011fe/'+SNname+'_uvotB14.1.dat', dt=dt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dates=[ '2011-08-28 04:00:00','2011-08-31 06:57:46','2011-09-03 10:00:01','2011-09-07 09:58:46','2011-09-10 10:02:09','2011-09-13 16:04:10','2011-09-19 14:09:10','2011-09-19 15:27:39','2011-10-01 06:46:19','2011-10-07 08:24:35','2011-10-21 05:02:26' ]
spectra=[ '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110828_firstuvepoch.dat',  '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110831_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110903_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110907_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110910_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110913_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110919_ccd.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_110919_mama.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_111001_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_111007_full.dat', '$SNFOLDER/hst/SN2011fe/myspec/SN2011fe_111021_full.dat' ]
myspecdet=[ 'mama',  'mama', 'ccd', 'ccd', 'ccd', 'ccd', 'ccd', 'mama', 'mama' , 'mama', 'mama' ]

ccd=where(myspecdet eq 'ccd')
mama=where(myspecdet eq 'mama')


mjdepochs=dblarr(n_elements(dates))
uvotepochmags=fltarr(6,n_elements(dates))
uvotepochmagerrs=fltarr(6,n_elements(dates))
hstepochmags=fltarr(6,n_elements(dates))
for n=0,n_elements(dates)-1 do mjdepochs[n]=date_conv(dates[n],'MODIFIED')

for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmags[f,n]=interpol(dt.mag_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])


for n=0,n_elements(dates)-1 do for f=0,5 do uvotepochmagerrs[f,n]=interpol(dt.magerr_array[f,where(finite(dt.mag_array[f,*]) eq 1)],dt.time_array[where(finite(dt.mag_array[f,*]) eq 1)],mjdepochs[n])


for n=0,n_elements(dates)-1 do begin
	uvot_specphot,spectra[n], specmags, speccounts=speccounts
	hstepochmags[*,n]=specmags[0:5]
endfor


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
xrange=[-15,45]
yrange1=[17.5,10.5]
yrange2=[-0.2,0.2]
yrange3=[-0.2,0.2]
xtitle1='Days from Maximum Light'

ytitle3='m2 '
ytitle2='w1 residuals'
ytitle1='Vega Mags'
nxticks=4
nyticks=4
nyticks1=7
nyticks2=4
nyticks3=4
figurename='SN2011fe_m2w1res.eps'

fontsize=14
SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

ploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], /nodata, /noerase, position=[x1,y1+(y2-y1),x2,y1+(y2-y1)+y2-y1], $
xtitle=' ',   ytitle=ytitle1, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange1, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks1, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1), symsize=0.2

oploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(16), color=cgcolor('pink'), errcolor=cgcolor('pink'), symsize=0.4
oploterror, dt.m2[0,*]-Bpeaktime[0], dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(15), color=cgcolor('powder blue'), errcolor=cgcolor('powder blue'), symsize=0.4

;oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*], symsize=1, psym=cgsymcat(6), color=cgcolor('blue')
;oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*], symsize=1, psym=cgsymcat(9), color=cgcolor('red')


oplot, mjdepochs[ccd]-Bpeaktime[0], hstepochmags[1,ccd], symsize=1, psym=cgsymcat(6), color=cgcolor('blue')
oplot, mjdepochs[mama]-Bpeaktime[0], hstepochmags[1,mama], symsize=1, psym=cgsymcat(15), color=cgcolor('blue')



oplot, mjdepochs[ccd]-Bpeaktime[0], hstepochmags[2,ccd], symsize=1, psym=cgsymcat(9), color=cgcolor('red')
oplot, mjdepochs[mama]-Bpeaktime[0], hstepochmags[2,mama], symsize=1, psym=cgsymcat(16), color=cgcolor('red')


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
ploterror, dt.w1[0,*]-Bpeaktime[0], w1interpedtemplatemags-dt.w1[1,*], dt.w1[2,*], xrange=[-15,45], yrange=[-0.2,0.2], /nodata, /noerase, position=[x1,y1+(y2-y1)/2,x2,y1+(y2-y1)], $
xtitle=' ', ytitle=ytitle2, charsize=1.0,  $
ystyle=1,xstyle=1, color=cgcolor('pink'), errcolor=cgcolor('pink'), $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


oploterror, dt.w1[0,*]-Bpeaktime[0], w1interpedtemplatemags-dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(16), color=cgcolor('pink'), errcolor=cgcolor('pink'), symsize=0.4

w1hstinterpedtemplatemags=fltarr(n_elements(dates))
for n=0,n_elements(w1hstinterpedtemplatemags)-1 do w1hstinterpedtemplatemags[n]=interpol(w1templatemags+11.02-0.08,w1templateepoch+55812.2-0.3,mjdepochs[n])


oplot, mjdepochs[ccd]-Bpeaktime[0], w1hstinterpedtemplatemags[ccd]-hstepochmags[2,ccd], symsize=1, psym=cgsymcat(9), color=cgcolor('red')
oplot, mjdepochs[mama]-Bpeaktime[0], w1hstinterpedtemplatemags[mama]-hstepochmags[2,mama], symsize=1, psym=cgsymcat(16), color=cgcolor('red')

; ploterror, dt.w1[0,*]-Bpeaktime[0], dt.w1[1,*], dt.w1[2,*], xrange=[-15,-10]
;oplot, mjdepochs-Bpeaktime[0], w1hstinterpedtemplatemags


f=1
	readcol, templatefiles[f], m2templateepoch, m2templatemags, /silent
	m2templateepoch=m2templateepoch+templateepochshift[f]
	m2templatemags=m2templatemags+templatemagshift[f]

m2interpedtemplatemags=fltarr(n_elements(dt.m2time))
for n=0,n_elements(m2interpedtemplatemags)-1 do m2interpedtemplatemags[n]=interpol(m2templatemags+13.06-0.1,m2templateepoch+55814.3-0.3,dt.m2time[n])


;;;;;;;;;;
ploterror, dt.m2[0,*]-Bpeaktime[0], m2interpedtemplatemags-dt.m2[1,*], dt.m2[2,*], xrange=[-15,45], yrange=[-0.2,0.2], /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle, ytitle=ytitle3, charsize=1.0,  $
ystyle=1,xstyle=1, color=cgcolor('powder blue'), errcolor=cgcolor('powder blue'),$
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

oploterror, dt.m2[0,*]-Bpeaktime[0], m2interpedtemplatemags-dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(15), color=cgcolor('powder blue'), errcolor=cgcolor('powder blue'), symsize=0.4


hstm2interpedtemplatemags=fltarr(n_elements(dates))
for n=0,n_elements(hstm2interpedtemplatemags)-1 do hstm2interpedtemplatemags[n]=interpol(m2templatemags+13.06-0.1,m2templateepoch+55814.3-0.3,mjdepochs[n])


oplot, mjdepochs[ccd]-Bpeaktime[0], hstm2interpedtemplatemags[ccd]-hstepochmags[1,ccd], symsize=1, psym=cgsymcat(6), color=cgcolor('blue')
oplot, mjdepochs[mama]-Bpeaktime[0], hstm2interpedtemplatemags[mama]-hstepochmags[1,mama], symsize=1, psym=cgsymcat(15), color=cgcolor('blue')



al_legend, ['UVOT uvw1','STIS-CCD w1','STIS-MAMA w1','UVOT uvm2','STIS-CCD m2','STIS-MAMA m2'], color=[cgcolor('pink'), cgcolor('red'),cgcolor('red') ,cgcolor('powder blue') , cgcolor('blue'), cgcolor('blue') ], psym=[cgcolor(16),cgcolor(9),cgcolor(16),cgcolor(15),cgcolor(6),cgcolor(15)], symsize=[0.4,1,1,0.4,1,1], position=[24.0,1.4], box=0, charsize=0.7

device, /close
SET_PLOT, 'X'
$open SN2011fe_m2w1res.eps 

stop
end
