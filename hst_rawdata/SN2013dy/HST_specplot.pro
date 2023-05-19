pro HST_specplot



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
xsize = 17.6

; putting these in real size (centimeters) 
; so they are the same regardless of the number of plots
wall = 0.1*xsize ; these values were originally made for the 8.8 width
margin=0.1*xsize ; a little big if the y-axis label is single digits 
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
yticklen = ticklen/a*100.0

xdata=[0,0,0,0]
ydata=[0,0,0,0]


fontsize=12



figurename='SN2013dy_specplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

nx=0
ny=0
;;;;;;;;;;;


readcol, 'sn2013dy-visit1-hst.flm', wavelength, flux



cgplot, xdata, ydata, /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='log(Observed Flux)+Constant', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[1400,8000], yrange=[-10,22]

axis, yaxis=1, ytitle='Days from Maximum Light', yrange=[-10,22]

cgoplot, wavelength, alog10(flux)-alog10(mean(flux[where(abs(wavelength-8000.0) lt 50)]))-6.6, color='red'

readcol, 'sn2013dy-visit2-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))-2.5, color='dark green'

readcol, 'sn2013dy-visit3-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))-0.8, color='blue'


readcol, 'sn2013dy-visit4-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+1.2, color='blue'

readcol, 'sn2013dy-visit5-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+4.5, color='red'

readcol, 'sn2013dy-visit6-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+8.4, color='dark green'

readcol, 'sn2013dy-visit7-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+12.0, color='orange'

readcol, 'sn2013dy-visit8-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+14.0, color='violet'

readcol, 'sn2013dy-visit9-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+17.9, color='red'

readcol, 'sn2013dy-visit10-hst.flm', wavelength, flux

cgoplot, wavelength, alog10(flux)-alog10(mean(flux(where(abs(wavelength-8000.0) lt 50))))+20.8, color='blue'

;cgtext, -15, 1.47, 'SN1992A', color='red', charsize=0.8

;;;;;;;;;;;
device, /close
SET_PLOT, 'X'
$open SN2013dy_specplot.eps


print, 'final stop'
stop
end