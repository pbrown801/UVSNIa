pro HST_2specplot



nxplots=2
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
yticklen = ticklen/a*100.0

xdata=[0,0,0,0]
ydata=[0,0,0,0]


fontsize=12



figurename='SN2013dy_2specplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'
;;  'R!l!7k!3!N'

nx=0
ny=0
;;;;;;;;;;;


readcol, 'sn2013dy-visit1-hst.flm', wavelength, flux



cgplot, wavelength, flux, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='Observed Flux', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[1400,8000], /ylog, yrange=[10.0^(-18), 10^(-13.0)]

cgoplot, wavelength, flux, color='red'

readcol, 'sn2013dy-visit2-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='dark green'

readcol, 'sn2013dy-visit3-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='blue'



nx=1

readcol, 'sn2013dy-visit3-hst.flm', wavelength, flux


cgplot, wavelength, flux, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='Observed Flux', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[1400,8000], yrange=[10.0^(-18), 10^(-13.0)], /ylog

readcol, 'sn2013dy-visit3-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='blue'

readcol, 'sn2013dy-visit4-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='red'

readcol, 'sn2013dy-visit5-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='dark green'

readcol, 'sn2013dy-visit6-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='orange'

readcol, 'sn2013dy-visit7-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='violet'

readcol, 'sn2013dy-visit8-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='red'

readcol, 'sn2013dy-visit9-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='dark green'

readcol, 'sn2013dy-visit10-hst.flm', wavelength, flux

cgoplot, wavelength, flux, color='blue'

;cgtext, -15, 1.47, 'SN1992A', color='red', charsize=0.8

;;;;;;;;;;;
device, /close
SET_PLOT, 'X'
$open SN2013dy_2specplot.eps


print, 'final stop'
stop
end