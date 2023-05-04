pro plot_spectra, SNname,  spectrum_array, color_array, epoch_array, yrange, uv_scale, optical_scale


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


;;;;;;;;;;;;;;;;;;;;;;;;;

figurename=SNname+'_fullspecplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'

nx=0
ny=0
;;;;;;;;;;;


cgplot, xdata, ydata, /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='Constantxlog(Observed Flux)+Constant', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[1400,8000], yrange=yrange

axis, yaxis=1, ytitle='Days from Maximum Light', yrange=[-10,22]

for y=0,n_elements(spectrum_array)-1 do begin
	readcol,spectrum_array[y] , wavelength, flux
	cgoplot, wavelength, uv_scale*alog10(flux)-uv_scale*alog10(mean(flux[where(abs(wavelength-8000.0) lt 50)]))+epoch_array[y], color=color_array[y]

endfor

;;;;;;;;;;;
device, /close
SET_PLOT, 'X'
spawn, 'open '+SNname+'_fullspecplot.eps'

figurename=SNname+'_uvspecplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'

nx=0
ny=0
;;;;;;;;;;;


cgplot, xdata, ydata, /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='2xlog(Observed Flux)+Constant', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[1400,4000], yrange=yrange

axis, yaxis=1, ytitle='Days from Maximum Light', yrange=[-10,22]

for y=0,n_elements(spectrum_array)-1 do begin
	readcol,spectrum_array[y] , wavelength, flux
	cgoplot, wavelength, uv_scale*alog10(flux)-uv_scale*alog10(mean(flux[where(abs(wavelength-8000.0) lt 50)]))+epoch_array[y], color='grey'

	newwave=fltarr( (max(wavelength)-min(wavelength))/10.0)
	newflux=fltarr( (max(wavelength)-min(wavelength))/10.0)
	for x=0,n_elements(newwave)-1 do newwave[x]=ceil(min(wavelength)/5.0)*5.0+x*10.0
	for x=0,n_elements(newwave)-1 do newflux[x]=mean(flux[where(wavelength lt newwave[x]+5.0 and wavelength gt newwave[x]-5.0)])

cgoplot, newwave, uv_scale*alog10(newflux)-uv_scale*alog10(mean(newflux[where(abs(newwave-8000.0) lt 50)]))+epoch_array[y], color=color_array[y]



endfor

;;;;;;;;;;;
device, /close
SET_PLOT, 'X'
spawn, 'open '+SNname+'_uvspecplot.eps'
figurename=SNname+'_opticalspecplot.eps'

SET_PLOT, 'PS'

device, filename=figurename, xsize=xsize, ysize=ysize, /encapsulated, $ 
font_size=fontsize, bits_per_pixel=8, /color, /tt_font, set_font='Times'

nx=0
ny=0
;;;;;;;;;;;


cgplot, xdata, ydata, /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle='Wavelength [Angstroms]',   ytitle='4xlog(Observed Flux)+Constant', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen/10, yticklen=yticklen/10, $
ystyle=1, xrange=[4000,8000], yrange=yrange

axis, yaxis=1, ytitle='Days from Maximum Light', yrange=[-10,22]

for y=0,n_elements(spectrum_array)-1 do begin
	readcol,spectrum_array[y] , wavelength, flux
	cgoplot, wavelength, optical_scale*alog10(flux)-optical_scale*alog10(mean(flux[where(abs(wavelength-8000.0) lt 50)]))+epoch_array[y], color=color_array[y]

endfor

;;;;;;;;;;;
device, /close
SET_PLOT, 'X'
spawn, 'open '+SNname+'_opticalspecplot.eps'

end

pro HST_specplot

SNname='SN2013dy'
spectrum_array=['sn2013dy-visit1-hst.flm','sn2013dy-visit2-hst.flm','sn2013dy-visit3-hst.flm','sn2013dy-visit4-hst.flm','sn2013dy-visit5-hst.flm','sn2013dy-visit6-hst.flm','sn2013dy-visit7-hst.flm','sn2013dy-visit8-hst.flm','sn2013dy-visit9-hst.flm','sn2013dy-visit10-hst.flm']
color_array=['red', 'dark green', 'blue', 'orange', 'violet', 'red', 'dark green', 'blue', 'orange', 'violet', 'red', 'blue']
epoch_array=[-6.6,-2.5,-0.8,1.2,4.5,8.4, 12, 14, 17.9, 20.8]
yrange=[floor(min(epoch_array))-5.0,floor(max(epoch_array))+4.0]
uv_scale=2.0
optical_scale=4.0

plot_spectra, SNname,  spectrum_array, color_array, epoch_array, yrange, uv_scale, optical_scale

;;;;;;;;;;;;;;

SNname='SN2011iv'
spectrum_array=['sn2011iv-visit1-hst.flm','sn2011iv-visit2-hst.flm','sn2011iv-visit3-hst.flm','sn2011iv-visit4-hst.flm','sn2011iv-visit5-hst.flm','sn2011iv-visit6-hst.flm','sn2011iv-visit7-hst.flm']
color_array=['red', 'dark green', 'blue', 'orange', 'violet', 'red', 'dark green', 'blue', 'orange', 'violet', 'red', 'blue']
epoch_array=[0.4,4.4,9.4,13.4,17.4,21.4,29.4]
yrange=[floor(min(epoch_array))-5.0,floor(max(epoch_array))+4.0]
uv_scale=2.0
optical_scale=4.0

plot_spectra, SNname,  spectrum_array, color_array, epoch_array, yrange, uv_scale, optical_scale

;;;;;;;;;;;;;;

SNname='SN2016ccj'
spectrum_array=['SN2016ccj_hst_20160514.dat','SN2016ccj_hst_20160521.dat','SN2016ccj_hst_20160529.dat']
color_array=['red', 'dark green', 'blue', 'orange', 'violet', 'red', 'dark green', 'blue', 'orange', 'violet', 'red', 'blue']
epoch_array=[0.0,7,15]
yrange=[floor(min(epoch_array))-5.0,floor(max(epoch_array))+4.0]
uv_scale=4.0
optical_scale=8.0

plot_spectra, SNname,  spectrum_array, color_array, epoch_array, yrange, uv_scale, optical_scale

;;;;;;;;;;;;;;

SNname='SN2017erp'
spectrum_array=['SN2017erp_hst_20170629.dat','SN2017erp_hst_20170702.dat','SN2017erp_hst_20170707.dat','SN2017erp_hst_20170712.dat']
color_array=['red', 'dark green', 'blue', 'orange', 'violet', 'red', 'dark green', 'blue', 'orange', 'violet', 'red', 'blue']
epoch_array=[-1.0,2,7,14]
yrange=[floor(min(epoch_array))-5.0,floor(max(epoch_array))+4.0]
uv_scale=1.50
optical_scale=4.0

plot_spectra, SNname,  spectrum_array, color_array, epoch_array, yrange, uv_scale, optical_scale

;;;;;;;;;;;;;;






print, 'final stop'
stop
end