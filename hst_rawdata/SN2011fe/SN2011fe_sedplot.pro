pro SN2011fe_sedplot

spec='myspec/SN2011fe_110831_full.dat'


	readcol,spec,specwavelengths,specfluxes,/silent

	pjb_uvotspec_all, spec, mag_array=mag_array, counts_array=counts_array

	filtereffwavelength=[1930,2200,2600,3450,4350,5460]

	stellarfluxfactors=[7.024, 9.179, 4.027, 1.666, 1.391, 2.728]*10.0^(-16.0)

	sedfluxes=counts_array[0:5]*stellarfluxfactors

plot, specwavelengths, specfluxes, /ylog, yrange=[10.0^(-16.0),10.0^(-12.0)], xrange=[1500,10000], /xlog
oplot, filtereffwavelength, sedfluxes, psym=4, symsize=2



nplots=1
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.04
margin=0.12
margin=0.17
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

figurename='SN2011fe_sed.eps'

fontsize=14
SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

f=0


;;;;;;;;;;
cgplot, specwavelengths/1000, specfluxes,position=[x1,y1,x2,y2], $
xtitle='Wavelength [10^3 Angstroms]', ytitle='Flux Density', charsize=1.0,  $
/ylog, yrange=[10.0^(-16.0),10.0^(-12.0)], xrange=[1.5,11], /xlog, ystyle=1,xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=5, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues

oplot, filtereffwavelength/1000.0, sedfluxes, psym=4, symsize=2


device, /close
SET_PLOT, 'X'
$open SN2011fe_sed.eps 





stop
end