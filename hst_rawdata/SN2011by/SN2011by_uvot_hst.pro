
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro SN2011by_uvot_hst


;;;; read in uvot data
SNname='SN2011by'
snphot_array, '$SNFOLDER/SwiftSNarchive/PhotArchive/'+SNname+'_uvotB13.dat', '$SNFOLDER/UVOTDATA/'+SNname+'_uvotB13.fits', dt=dt

dates=['2011-04-30 11:33:41','2011-05-05 03:00:06','2011-05-09 06:15:47','2011-05-13 17:39:31','2011-05-18 04:32:31']
spectra=["SN2011by_0430_full.dat","SN2011by_0505_opt.dat","SN2011by_0509_full.dat","SN2011by_0513_opt.dat","SN2011by_0518_opt.dat"]
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
;	addsedendpoints, uvotepochmags[*,n], sedwavelengths, sedfluxes, outputsedwavelengths, outputsedfluxes
;	uvot_specphot, [transpose(outputsedwavelengths),transpose(outputsedfluxes)], specmags, speccounts  

if sedwavelengths[0] gt 1600 then sedfluxes=[sedfluxes[0],sedfluxes]
if sedwavelengths[n_elements(sedwavelengths)-1] lt 8000 then sedfluxes=[sedfluxes,sedfluxes[n_elements(sedfluxes)-1]]
if sedwavelengths[0] gt 1600 then sedwavelengths=[1500,sedwavelengths]
if sedwavelengths[n_elements(sedwavelengths)-1] lt 8000 then sedwavelengths=[sedwavelengths,9000]
uvot_specphot, [transpose(sedwavelengths),transpose(sedfluxes)], specmags, speccounts  


for f=0,5 do if sedwavelengths[1] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-2] lt wave_long[f] then speccounts[f]= !Values.F_NAN
for f=0,5 do if sedwavelengths[1] gt wave_short[f] or sedwavelengths[n_elements(sedwavelengths)-2] lt wave_long[f] then specmags[f]= !Values.F_NAN
	hstepochmags[*,n]=specmags[0:5]


endfor


maghigh=floor(max(dt.mag_array,/nan))+1.0
maglow=floor(min(dt.mag_array,/nan))
mjdhigh=floor(max(dt.time_array,/nan)/10.0)*10.0
mjdlow=floor(min(dt.time_array,/nan)/10.0)*10.0


; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
wall = 0.03
margin=0.12
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall)*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
y3 = y2 + wall*8.8/ysize
y4 = y3 + b*8.8/ysize
yc = y4 + wall*8.8/ysize
fontsize=14

xdata=[1,2,3,4]
ydata=[1,2,3,4]
Bpeaktime=dt.time_array[ where(dt.mag_array[4,*] eq min(dt.mag_array[4,*], /nan) ) ]


SET_PLOT, 'PS'

device, filename=SNname+'_uvot_hst.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, xdata,ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle=' ',   ytitle='UVOT Vega Mags', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[mjdlow-Bpeaktime[0],mjdhigh-Bpeaktime[0]],yrange=[maghigh,maglow], ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



;oploterror, dt.w2[0,where(finite(dt.w2[1,*]) eq 1)]-Bpeaktime[0],dt.w2[1,where(finite(dt.w2[1,where(finite(dt.w2[1,*]) eq 1)]) eq 1)],dt.w2[2,*], psym=3
oploterror, dt.m2[0,where(finite(dt.m2[1,*]) eq 1)]-Bpeaktime[0],dt.m2[1,where(finite(dt.m2[1,where(finite(dt.m2[1,*]) eq 1)]) eq 1)],dt.m2[2,*], psym=-5 
oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-Bpeaktime[0],dt.w1[1,where(finite(dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]) eq 1)],dt.w1[2,*], psym=-6
oploterror, dt.uu[0,where(finite(dt.uu[1,*]) eq 1)]-Bpeaktime[0],dt.uu[1,where(finite(dt.uu[1,where(finite(dt.uu[1,*]) eq 1)]) eq 1)],dt.uu[2,*], psym=-4
oploterror, dt.bb[0,where(finite(dt.bb[1,*]) eq 1)]-Bpeaktime[0],dt.bb[1,where(finite(dt.bb[1,where(finite(dt.bb[1,*]) eq 1)]) eq 1)],dt.bb[2,*], psym=-5
;oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-Bpeaktime[0],dt.vv[1,where(finite(dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]) eq 1)],dt.vv[2,*], psym=3


oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*], psym=-5, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*], psym=-6, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[3,*], psym=-4, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[4,*], psym=-5, symsize=2


al_legend, ['b','u','w1','m2'],   psym=[5,4,6,5],  background_color='white',$
pos=[0.8,0.9], /norm, charsize=1.0, box=0


device, /close
SET_PLOT, 'X'
$open SN2011by_uvot_hst.eps


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Bpeaktime=dt.time_array[ where(dt.mag_array[4,*] eq min(dt.mag_array[4,*], /nan) ) ]
	w1b=where( finite(dt.mag_array[2,*]) eq 1 and finite(dt.mag_array[4,*]) eq 1 )
	m2w1=where( finite(dt.mag_array[1,*]) eq 1 and finite(dt.mag_array[2,*]) eq 1 )
colorhigh=floor(max([reform(dt.mag_array[1,m2w1]-dt.mag_array[2,m2w1]),reform(dt.mag_array[2,w1b]-dt.mag_array[4,w1b])],/nan)*20.0)/20.0+0.2
colorlow=floor(min([reform(dt.mag_array[1,m2w1]-dt.mag_array[2,m2w1]),reform(dt.mag_array[2,w1b]-dt.mag_array[4,w1b])],/nan)*20.0)/20.0-0.2
SET_PLOT, 'PS'

device, filename=SNname+'_uvot_hst_colors.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=10, bits_per_pixel=8, /color

plot, xdata,ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Days from Maximum Light',   ytitle='UVOT colors', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[mjdlow-Bpeaktime,mjdhigh-Bpeaktime],yrange=[colorlow,colorhigh], ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



oploterror, dt.time_array[w1b]-Bpeaktime[0], dt.mag_array[2,w1b]-dt.mag_array[4,w1b], sqrt(dt.magerr_array[2,w1b]^2.0+dt.magerr_array[4,w1b]^2.0), linestyle=2, psym=-5



oploterror, dt.time_array[m2w1]-Bpeaktime[0], dt.mag_array[1,m2w1]-dt.mag_array[2,m2w1], sqrt(dt.magerr_array[1,m2w1]^2.0+dt.magerr_array[2,m2w1]^2.0), linestyle=2, psym=-6


oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*]-hstepochmags[4,*], psym=-5, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*]-hstepochmags[2,*], psym=-6, symsize=2

al_legend, ['w1-b','m2-w1'],   psym=[5,6],  background_color='white',$
pos=[0.8,0.9], /norm, charsize=1.0, box=0


device, /close
SET_PLOT, 'X'
$open SN2011by_uvot_hst_colors.eps

;;;;;;;;;;;;;; light curves and colors
nplots=2
ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a
nxticks=10
nyticks=10

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize

nplots=2
figurename='SN2011by_lc_clrs.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
xtitle='Days from Maximum Light',   ytitle='UVOT Vega Mags', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[mjdlow-Bpeaktime[0],mjdhigh-Bpeaktime[0]],yrange=[maghigh,maglow], ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



;oploterror, dt.w2[0,where(finite(dt.w2[1,*]) eq 1)]-Bpeaktime[0],dt.w2[1,where(finite(dt.w2[1,where(finite(dt.w2[1,*]) eq 1)]) eq 1)],dt.w2[2,*], psym=3
oploterror, dt.m2[0,where(finite(dt.m2[1,*]) eq 1)]-Bpeaktime[0],dt.m2[1,where(finite(dt.m2[1,where(finite(dt.m2[1,*]) eq 1)]) eq 1)],dt.m2[2,*], psym=-5 
oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-Bpeaktime[0],dt.w1[1,where(finite(dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]) eq 1)],dt.w1[2,*], psym=-6
oploterror, dt.uu[0,where(finite(dt.uu[1,*]) eq 1)]-Bpeaktime[0],dt.uu[1,where(finite(dt.uu[1,where(finite(dt.uu[1,*]) eq 1)]) eq 1)],dt.uu[2,*], psym=-4
oploterror, dt.bb[0,where(finite(dt.bb[1,*]) eq 1)]-Bpeaktime[0],dt.bb[1,where(finite(dt.bb[1,where(finite(dt.bb[1,*]) eq 1)]) eq 1)],dt.bb[2,*], psym=-5
;oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-Bpeaktime[0],dt.vv[1,where(finite(dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]) eq 1)],dt.vv[2,*], psym=3


oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*], psym=-5, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*], psym=-6, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[3,*], psym=-4, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[4,*], psym=-5, symsize=2


al_legend, ['b','u','w1','m2'],   psym=[5,4,6,5],  background_color='white',$
pos=[0.8,0.9], /norm, charsize=1.0, box=0

plot, xdata,ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Days from Maximum Light',   ytitle='UVOT colors', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[mjdlow-Bpeaktime,mjdhigh-Bpeaktime],yrange=[colorlow,colorhigh], ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues


	w1b=where( finite(dt.mag_array[2,*]) eq 1 and finite(dt.mag_array[4,*]) eq 1 )

oploterror, dt.time_array[w1b]-Bpeaktime[0], dt.mag_array[2,w1b]-dt.mag_array[4,w1b], sqrt(dt.magerr_array[2,w1b]^2.0+dt.magerr_array[4,w1b]^2.0), linestyle=2, psym=-5


	m2w1=where( finite(dt.mag_array[1,*]) eq 1 and finite(dt.mag_array[2,*]) eq 1 )

oploterror, dt.time_array[m2w1]-Bpeaktime[0], dt.mag_array[1,m2w1]-dt.mag_array[2,m2w1], sqrt(dt.magerr_array[1,m2w1]^2.0+dt.magerr_array[2,m2w1]^2.0), linestyle=2, psym=-6


oplot, mjdepochs-Bpeaktime[0], hstepochmags[2,*]-hstepochmags[4,*], psym=-5, symsize=2
oplot, mjdepochs-Bpeaktime[0], hstepochmags[1,*]-hstepochmags[2,*], psym=-6, symsize=2

al_legend, ['w1-b','m2-w1'],   psym=[5,6],  background_color='white',$
pos=[0.8,0.9], /norm, charsize=1.0, box=0


device, /close
SET_PLOT, 'X'
$open SN2011by_lc_clrs.eps 




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
xrange=[55680.0,55700.0]
xrange=[55675.0-55000.0,55710.0-55000.0]
xrange=[55675.0-55000.0,55735.0-55000.0]
yrange2=[max(dt.mag_array,/nan)+1.0,min(dt.mag_array,/nan)-1.0]
yrange2=[21.0,11.0]
yrange1=[-0.2,0.3]
xtitle1='Modified Julian Date-55000'
ytitle1='UVOT-HST Residuals'
ytitle2='Vega Mags'
nxticks=6
nyticks=5
nyticks2=5
figurename='SN2011by_res.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, dt.vv[0,*]-55000.0, dt.vv[1,*], /nodata, /noerase, position=[x1,y1+(y2-y1)/2,x2,y1+(y2-y1)/2+y2-y1], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange2, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

oploterror, dt.w2[0,*]-55000.0, dt.w2[1,*], dt.w2[2,*], psym=cgsymcat(5), color=cgcolor('black'), errcolor=cgcolor('black')
oploterror, dt.m2[0,*]-55000.0, dt.m2[1,*], dt.m2[2,*], psym=cgsymcat(1), color=cgcolor('red'), errcolor=cgcolor('red')
oploterror, dt.w1[0,*]-55000.0, dt.w1[1,*], dt.w1[2,*], psym=cgsymcat(2), color=cgcolor('maroon'), errcolor=cgcolor('maroon')
oploterror, dt.uu[0,*]-55000.0, dt.uu[1,*], dt.uu[2,*], psym=cgsymcat(4), color=cgcolor('purple'), errcolor=cgcolor('purple')
oploterror, dt.bb[0,*]-55000.0, dt.bb[1,*], dt.bb[2,*], psym=cgsymcat(6), color=cgcolor('blue'), errcolor=cgcolor('blue')
oploterror, dt.vv[0,*]-55000.0, dt.vv[1,*], dt.vv[2,*], psym=cgsymcat(9), color=cgcolor('dark green'), errcolor=cgcolor('dark green')

oplot, mjdepochs-55000.0, hstepochmags[0,*], symsize=2, psym=cgsymcat(5), color=cgcolor('black')
oplot, mjdepochs-55000.0, hstepochmags[1,*], symsize=2, psym=cgsymcat(1), color=cgcolor('red')
oplot, mjdepochs-55000.0, hstepochmags[2,*], symsize=2, psym=cgsymcat(2), color=cgcolor('maroon')
oplot, mjdepochs-55000.0, hstepochmags[3,*], symsize=2, psym=cgsymcat(4), color=cgcolor('purple')
oplot, mjdepochs-55000.0, hstepochmags[4,*], symsize=2, psym=cgsymcat(6), color=cgcolor('blue')
oplot, mjdepochs-55000.0, hstepochmags[5,*], symsize=2, psym=cgsymcat(9), color=cgcolor('dark green')



plot, mjdepochs-55000.0, uvotepochmags[0,*]-hstepochmags[0,*], /nodata, /noerase, position=[x1,y1,x2,y1+(y2-y1)/2], $
xtitle=xtitle1, ytitle=ytitle1, charsize=1.0,  $
 yrange=yrange1, ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues


oploterror, mjdepochs-55000.0-0.6, uvotepochmags[0,*]-hstepochmags[0,*], uvotepochmagerrs[0,*],psym=cgsymcat(17), color=cgcolor('black'), errcolor=cgcolor('black'), thick=2, symsize=2
oploterror, mjdepochs-55000.0+0.2, uvotepochmags[1,*]-hstepochmags[1,*], uvotepochmagerrs[1,*],psym=cgsymcat(1), color=cgcolor('red'), errcolor=cgcolor('red'), thick=2, symsize=2
oploterror, mjdepochs-55000.0-0.2, uvotepochmags[2,*]-hstepochmags[2,*], uvotepochmagerrs[2,*],psym=cgsymcat(2), color=cgcolor('maroon'), errcolor=cgcolor('maroon'), thick=2, symsize=2
oploterror, mjdepochs-55000.0+0.4, uvotepochmags[3,*]-hstepochmags[3,*], uvotepochmagerrs[3,*],psym=cgsymcat(4), color=cgcolor('purple'), errcolor=cgcolor('purple')
oploterror, mjdepochs-55000.0, uvotepochmags[4,*]-hstepochmags[4,*], uvotepochmagerrs[4,*],psym=cgsymcat(6), color=cgcolor('blue'), errcolor=cgcolor('blue')
oploterror, mjdepochs-55000.0-0.4, uvotepochmags[5,*]-hstepochmags[5,*], uvotepochmagerrs[5,*],psym=cgsymcat(9), color=cgcolor('dark green'), errcolor=cgcolor('dark green')


al_legend, ['w2','m2','w1','u','b','v'], color=[cgcolor('black'), cgcolor('red'),cgcolor('maroon') ,cgcolor('purple') , cgcolor('blue'), cgcolor('dark green') ], psym=[5,1, 2,4,6, 9 ], position=[715.0,0.3], box=0, charsize=0.7

device, /close
SET_PLOT, 'X'
$open SN2011by_res.eps 





print, 'final stop'
stop
end

