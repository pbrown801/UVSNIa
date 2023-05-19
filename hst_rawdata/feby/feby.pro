pro feby
;; apply extinction to see if they look similar

byspec='SN2011by_0509_full.dat'
fespec='SN2011fe_110910_full.dat'

readcol,"$SNSCRIPTS/V_UVOT.txt",lambda,V_EA,/silent

	readcol,byspec,by_wave,by_flux,/silent
	readcol,fespec,fe_wave,fe_flux,/silent

byfluxlambda=interpol(by_flux,by_wave,lambda)
fefluxlambda=interpol(fe_flux,fe_wave,lambda)

plot, lambda, byfluxlambda
oplot, lambda, fefluxlambda

;function sne_mw_reddening,lambda,Ebv,rv=rv,z=z,av=av

bymwa=sne_mw_reddening(lambda,0.012)
femwa=sne_mw_reddening(lambda,0.008)
byfluxderedmw=byfluxlambda*10^(-bymwa/2.5)
fefluxderedmw=fefluxlambda*10^(-femwa/2.5)

plot, lambda, byfluxderedmw
oplot, lambda, fefluxderedmw


uvot_specphot, [transpose(lambda),transpose(fefluxderedmw)], femags, fecounts
uvot_specphot, [transpose(lambda),transpose(byfluxderedmw)], bymags, bycounts


bymwhosta=sne_mw_reddening(lambda,-0.15)
byfluxderedmwhost=byfluxderedmw*10^(-bymwhosta/2.5)

bymwhostg=sne_goobarlmc_reddening(lambda,-0.15)
byfluxderedmwhostg=byfluxderedmw*10^(-bymwhosta/2.5)


fefinalflux=fefluxlambda
byfinalflux=byfluxlambda
fefinalflux=fefluxderedmw
byfinalflux=byfluxderedmw
byfinalflux=byfluxderedmwhost
byfinalflux=byfluxderedmwhostg


factor=fefinalflux[where (lambda ge 4500 and lambda le 5000)]/byfinalflux[where (lambda ge 4500 and lambda le 5000)]

byfinalfluxnorm=byfinalflux*mean(factor)

plot, lambda, fefinalflux, /ylog, yrange=[10.0^(-15.0),10.0^(-12.0)], /ystyle
oplot, lambda, byfinalfluxnorm
;oplot, lambda, byfluxderedmwhost*100


plot, lambda, fefinalflux/byfinalfluxnorm

restore, '$SNFOLDER/sniamodels/sauerall/sauerallmodels.sav'
;neg plot, densitywave, densityspec_array[*,2]/densityspec_array[*,1], xrange=[1600,4200] 


;;; something similar
plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]

oplot, densitywave, densityspec_array[*,1]/densityspec_array[*,3]

plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,2]/densityspec_array[*,5] 


plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,3]/densityspec_array[*,6] 






plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,12]/densityspec_array[*,8]

;;; this gets the 3700 A dip right but the mid-UV is in the wrong direction
plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,8]/densityspec_array[*,12]


plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,2]/densityspec_array[*,10]


; okay 
plot, lambda, fefinalflux/byfinalfluxnorm, xrange=[1600,4200], yrange=[0.0,5]
oplot, densitywave, densityspec_array[*,2]/densityspec_array[*,12]

;;;;;;;;;;;;;




restore, filename='$SNFOLDER/sniamodels/lentz/lentzmodels.sav'
lentzdates=datenums-17.0
lentz_mag_array=mag_array
lentz_spec_array=spec_array
lentz_wave=wave
lentzsymbol=17
lentzblank=5

lentz7000=where(lentz_wave-7000.0 eq min(abs(lentz_wave - 7000.0)))

restore, filename='$SNFOLDER/sniamodels/walker/walkermodels.sav'
walker_n=n
walker_mag_array=mag_array
walker_spec_array=spec_array
walker_wave=wave
walkersymbol=15
walkerblank=6
walker7000=where(walker_wave-7000.0 eq min(abs(walker_wave - 7000.0)))

restore, filename='$SNFOLDER/sniamodels/sauer/sauermodels.sav'
sauer_mag_array=mag_array
sauer_spec_array=spec_array
sauer_wave=wave
sauersymbol=16
sauerblank=9
sauer7000=where(sauer_wave-7000.0 eq min(abs(sauer_wave - 7000.0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  2 x 3 PLOT  ;;;;;;;;;;;;;;;;;;;;

xplots=3
yplots=2
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
xsize = 18.0
wall = 0.03
; this was the margin normalized for 8.8 cm and font 14
margin=0.16
; bumped up for log scale
margin=0.18
a = (xsize - (margin*8.8 + xplots*wall*8.8))/3.0
b = a * 2d / (1 + sqrt(5))

ysize = (margin*8.8 + yplots*(b + wall*8.8 ) )
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a
fontsize=14

x1 = margin
x2 = x1 + a
x3 = x2 + a 
x4 = x3 + a 
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
xrange=[1500.0,7500.0]
xrangeleft=[1500.0,7500.0]
xrangemiddle=[1500.0,7500.0]
xrangeright=[1500.0,7500.0]

nxticks=6
nyticks=5

nytickstop=3
yrangetop=[0.01,10.0]
yrangebottom=[0.0,5.0]

figurename='spectralplots.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

plot, walker_wave, walker_spec_array[*,0,3], /nodata, /noerase, $
position=[(margin)*8.8/xsize,(margin)*8.8/ysize+1.0*(b+wall*8.8/ysize)/ysize,(margin)*8.8/xsize+a/xsize,(margin)*8.8/ysize+2.0*(b+wall*8.8/ysize)/ysize], $
xtitle=' ',   ytitle='Normalized Flux', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=yrangetop,  /ylog,ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nytickstop, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)
walkernorm=walker_spec_array[walker7000,0,3]
for n=0,5 do oplot, walker_wave, walker_spec_array[*,n,3]/walkernorm[0]


xyouts, 3500, 0.1, 'Walker'

good=where(lentz_wave lt 9000)

plot, walker_wave, walker_spec_array[*,0,3]/walker_spec_array[*,6,3], /nodata, /noerase, $
position=[(margin)*8.8/xsize, (margin)*8.8/ysize+0.0*(b+wall)/ysize,(margin)*8.8/xsize+a/xsize, (margin)*8.8/ysize+1.0*(b)/ysize], $
xtitle=' ', ytitle='Flux Ratio', charsize=1.0,  $
yrange=yrangebottom, ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

for n=0,3 do oplot, walker_wave, walker_spec_array[*,n,3]/walker_spec_array[*,6,3]
cgoplot, lambda, fefinalflux/byfinalfluxnorm, color='red'

;;; middle



plot, lentz_wave[good],lentz_spec_array[good,2,0], /nodata, /noerase, $
position=[(margin)*8.8/xsize+1.0*(a+wall)/xsize,(margin)*8.8/ysize+1.0*(b+wall)/ysize,(margin)*8.8/xsize+2.0*(a+wall)/xsize,(margin)*8.8/ysize+2.0*(b+wall)/ysize], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrangetop, /ylog, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nytickstop, ytickv=ytickvalues , $
 xtickname=replicate(' ',nxticks+1), ytickname=replicate(' ',nytickstop+1)
lentznorm=lentz_spec_array[lentz7000,2,0]

for n=0,4 do oplot, lentz_wave[good],lentz_spec_array[good,2,n]/lentznorm[0]

xyouts, 3500, 0.1, 'Lentz'


plot, lentz_wave[good],lentz_spec_array[good,2,0]/lentz_spec_array[good,2,4], /nodata, /noerase, $
position=[(margin)*8.8/xsize+1.0*(a+wall)/xsize,(margin)*8.8/ysize+0.0*(b+wall)/ysize,(margin)*8.8/xsize+2.0*(a+wall)/xsize, (margin)*8.8/ysize+1.0*(b+wall)/ysize], $
xtitle='Wavelength [Angstroms] ', ytitle=ytitle2, charsize=1.0,  $
 yrange=yrangebottom, ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, ytickname=replicate(' ',nyticks+1)

for n=0,3 do oplot, lentz_wave[good],lentz_spec_array[good,2,n]/lentz_spec_array[good,2,4]
cgoplot, lambda, fefinalflux/byfinalfluxnorm, color='red'

;;; right

plot, sauer_wave, sauer_spec_array[*,0], /nodata, /noerase, $
position=[(margin)*8.8/xsize+2.0*(a+wall)/xsize,(margin)*8.8/ysize+1.0*(b+wall)/ysize,(margin)*8.8/xsize+3.0*(a+wall)/xsize,(margin)*8.8/ysize+2.0*(b+wall)/ysize], $
xtitle=xtitle2,   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrangetop, /ylog, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nytickstop, ytickv=ytickvalues , $
xtickname=replicate(' ',nxticks+1), ytickname=replicate(' ',nytickstop+1)

sauernorm=sauer_spec_array[sauer7000,0]
for n=0,8 do oplot, sauer_wave, sauer_spec_array[*,n]/sauernorm[0]

xyouts, 3500, 0.1, 'Sauer'


plot, sauer_wave, sauer_spec_array[*,0]/sauer_spec_array[*,0], /nodata, /noerase, position=[margin*8.8/xsize+2.0*(a+wall)/xsize, (margin)*8.8/ysize+0.0*(b+wall)/ysize,(margin)*8.8/xsize+3.0*(a+wall)/xsize,(margin)*8.8/ysize+1.0*(b+wall)/ysize], $
xtitle=xtitle, ytitle=ytitle2, charsize=1.0,  $
 yrange=yrangebottom, ystyle=1, xrange=xrange, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, ytickname=replicate(' ',nyticks+1)

for n=0,8 do oplot, sauer_wave, sauer_spec_array[*,n]/sauer_spec_array[*,0]


cgoplot, lambda, fefinalflux/byfinalfluxnorm, color='red'

device, /close
SET_PLOT, 'X'
$open spectralplots.eps 





stop
end