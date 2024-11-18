PRO output, z, muvspectra, nuvspectra 

;z=0.0094

;muvspectra=['SN2021fxy_muv_20210329_1.dat', 'SN2021fxy_muv_20210329_2.dat', 'SN2021fxy_muv_20210329_3.dat','SN2021fxy_muv_20210330_1.dat','SN2021fxy_muv_20210330_2.dat' ]

;nuvspectra=['SN2021fxy_nuv_20210329_1.dat']

;output='SN2021fxy_hst_20210329_rest_'

n_muv=n_elements(muvspectra)

muvflux=fltarr(n_muv)

wavestart=1650
wavestop=5600

;for n=0, n_muv-1 do begin
;readcol, muvspectra[n], wave_uv1_1,  flux_uv1_1,  fluxerr_uv1_1,  dq_uv1_1,  format=('F,F,F,I')
;endfor

wave_uv1_1=0
flux_uv1_1=0
fluxerr_uv1_1=0
dq_uv1_1=0
dqbits_uv1_1=0

wave_uv1_2=0
flux_uv1_2=0
fluxerr_uv1_2=0
dq_uv1_2=0
dqbits_uv1_2=0


wave_uv1_3=0
flux_uv1_3=0
fluxerr_uv1_3=0
dq_uv1_3=0
dqbits_uv1_3=0


wave_uv1_4=0
flux_uv1_4=0
fluxerr_uv1_4=0
dq_uv1_4=0
dqbits_uv1_4=0

wave_uv1_5=0
flux_uv1_5=0
fluxerr_uv1_5=0
dq_uv1_5=0
dqbits_uv1_5=0



;;;;; use error column/ data quality

                   readcol, nuvspectra[0], wave_opt1_1, flux_opt1_1, fluxerr_opt1_1, dq_opt1_1, format=('F,F,F,I')

                   readcol, muvspectra[0], wave_uv1_1,  flux_uv1_1,  fluxerr_uv1_1,  dq_uv1_1,  format=('F,F,F,I')

if n_muv gt 1 then readcol, muvspectra[1], wave_uv1_2,  flux_uv1_2,  fluxerr_uv1_2,  dq_uv1_2,  format=('F,F,F,I')

if n_muv gt 2 then readcol, muvspectra[2], wave_uv1_3,  flux_uv1_3,  fluxerr_uv1_3,  dq_uv1_3,  format=('F,F,F,I')

if n_muv gt 3 then readcol, muvspectra[3], wave_uv1_4,  flux_uv1_4,  fluxerr_uv1_4,  dq_uv1_4,  format=('F,F,F,I')

if n_muv gt 4 then readcol, muvspectra[4], wave_uv1_5,  flux_uv1_5,  fluxerr_uv1_5,  dq_uv1_5,  format=('F,F,F,I')




                   dqbits_opt1_1=strarr(16,n_elements(dq_opt1_1))
                   dqbits_uv1_1 =strarr(16,n_elements(dq_uv1_1))
if n_muv gt 1 then dqbits_uv1_2 =strarr(16,n_elements(dq_uv1_2))
if n_muv gt 2 then dqbits_uv1_3 =strarr(16,n_elements(dq_uv1_3))
if n_muv gt 3 then dqbits_uv1_4 =strarr(16,n_elements(dq_uv1_4))
if n_muv gt 4 then dqbits_uv1_5 =strarr(16,n_elements(dq_uv1_5))

                   for n=0,n_elements(dq_opt1_1)-1 do dqbits_opt1_1[*,n]=binary(dq_opt1_1[n])
                   for n=0,n_elements(dq_uv1_1)-1  do dqbits_uv1_1[*,n]=binary(dq_uv1_1[n])
if n_muv gt 1 then for n=0,n_elements(dq_uv1_2)-1  do dqbits_uv1_2[*,n]=binary(dq_uv1_2[n])
if n_muv gt 2 then for n=0,n_elements(dq_uv1_3)-1  do dqbits_uv1_3[*,n]=binary(dq_uv1_3[n])
if n_muv gt 3 then for n=0,n_elements(dq_uv1_4)-1  do dqbits_uv1_4[*,n]=binary(dq_uv1_4[n])
if n_muv gt 4 then for n=0,n_elements(dq_uv1_5)-1  do dqbits_uv1_5[*,n]=binary(dq_uv1_5[n])



avgflux_uv1_1=flux_uv1_1
for f=3,n_elements(flux_uv1_1)-3 do avgflux_uv1_1[f]=mean(flux_uv1_1[f-2:f+2])

medflux_uv1_1=flux_uv1_1
for f=3,n_elements(flux_uv1_1)-3 do medflux_uv1_1[f]=median(flux_uv1_1[f-2:f+2])

stdflux_uv1_1=fluxerr_uv1_1
for f=3,n_elements(flux_uv1_1)-3 do stdflux_uv1_1[f]=stddev(flux_uv1_1[f-2:f+2])



                   uv1good=[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')]
if n_muv gt 1 then uv2good=[where(dqbits_uv1_2[4,*] ne '1' and dqbits_uv1_2[6,*] ne '1' and dqbits_uv1_2[13,*] ne '1')]
if n_muv gt 2 then uv3good=[where(dqbits_uv1_3[4,*] ne '1' and dqbits_uv1_3[6,*] ne '1' and dqbits_uv1_3[13,*] ne '1')]
if n_muv gt 3 then uv4good=[where(dqbits_uv1_4[4,*] ne '1' and dqbits_uv1_4[6,*] ne '1' and dqbits_uv1_4[13,*] ne '1')]
if n_muv gt 4 then uv5good=[where(dqbits_uv1_5[4,*] ne '1' and dqbits_uv1_5[6,*] ne '1' and dqbits_uv1_5[13,*] ne '1')]

opt1good=where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )


;;;;;;;;  this gives a good spectrum
cgplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )],   flux_opt1_1[[where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )]], psym=14, color='blue', xrange=[2800,6000]
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1' or dqbits_opt1_1[6,*] eq '1' )],   flux_opt1_1[where(dqbits_opt1_1[11,*] eq '1' or dqbits_opt1_1[6,*] eq '1' )], psym=15, color='red'

oplot, [wavestart,wavestart], [0,10^(-12.0)], psym=-5,color='black'
oplot, [wavestop,wavestop], [0,10^(-12.0)], psym=-5, color='black'
print, "change wavestart and wavestop if you want"

stop

cgplot,   wave_uv1_1[uv1good],    flux_uv1_1[uv1good],   psym=16, color='green'
oplot, [wavestart,wavestart], [0,10^(-12.0)], psym=-5,color='black'
oplot, [wavestop,wavestop], [0,10^(-12.0)], psym=-5, color='black'
print, "change wavestart and wavestop if you want"

stop

if n_muv gt 1 then print, 'second muv spectrum '
if n_muv gt 1 then cgplot,   wave_uv1_2[uv2good],    flux_uv1_2[uv2good],   psym=16, color='green'
if n_muv gt 1 then stop

if n_muv gt 2 then print, 'third muv spectrum '
if n_muv gt 2 then cgplot,   wave_uv1_3[uv3good],    flux_uv1_3[uv3good],   psym=16, color='green'
if n_muv gt 2 then stop




;cgplot,   wave_uv1_1[uv1good and abs(flux_uv1_1-avgflux_uv1_1)<2*stdflux_uv1_1],    flux_uv1_1[uv1good and abs(flux_uv1_1-avgflux_uv1_1)<2*stdflux_uv1_1],   psym=16, color='green'


;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1' and abs(flux_uv1_1-avgflux_uv1_1) lt 2*stdflux_uv1_1)],    flux_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1' and abs(flux_uv1_1-avgflux_uv1_1) lt 2*stdflux_uv1_1)],   color='red'





 mg2= [2796.0,2803] 
;; weak mg2 from MW, nothing obvious between z=0.02-0.03



;1600 - 5680 combining both spectra


; 1 0000 0000 0000 0001 Error in the Reed-Solomon decoding (an algorithm for error correction in digital communications).
; 2 0000 0000 0000 0010 Lost data replaced by fill values.
; 4 0000 0000 0000 0100 Bad detector pixel (e.g., bad column or row, mixed science and bias for overscan, or beyond aperture).
;;; 8 0000 0000 0000 1000 Data masked by occulting bar. 
; 16 0000 0000 0001 0000 Pixel having dark rate > 5 σ times the median dark level.
;;;;; 32 0000 0000 0010 0000 Large blemish, depth > 40% of the normalized p-flat (repeller wire).
; 64 0000 0000 0100 0000 Vignetted pixel
; 128 0000 0000 1000 0000 Pixel in the overscan region.
; 256 0000 0001 0000 0000 Saturated pixel, count rate at 90% of max possible—local non-linearity turns over and is multi-valued; 
;				pixels within 10% of turnover and all pixels within 4 pixels of that pixel are flagged.
; 512 0000 0010 0000 0000 Bad pixel in reference file.
; 1024 0000 0100 0000 0000 Small blemish, depth between 40% and 70% of the normalized flat. Applies to MAMA and CCD p-flats.
; 2048 0000 1000 0000 0000 >30% of background pixels rejected by sigma-clip, or flagged, during 1-D spectral extraction.
;;;; 4096 0001 0000 0000 0000 Extracted flux affected by bad input data.
; 8192 0010 0000 0000 0000 ata rejected in input pixel during image combination for cosmic ray rejection.
; 16384 0100 0000 0000 0000 Extracted flux not CTI corrected because gross counts are ≤ 0.




if n_muv eq 1 then flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good] ]

if n_muv gt 1 then flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good] ]

if n_muv gt 2 then flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good], flux_uv1_3[uv3good] ]

if n_muv gt 3 then flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good], flux_uv1_3[uv3good], flux_uv1_4[uv4good] ]

if n_muv gt 4 then flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good], flux_uv1_3[uv3good], flux_uv1_4[uv4good], flux_uv1_5[uv5good]]


if n_muv eq 1 then fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good] ]

if n_muv gt 1 then fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good] ]

if n_muv gt 2 then fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good] ]

if n_muv gt 3 then fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good], fluxerr_uv1_4[uv4good] ]

if n_muv gt 4 then fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good], fluxerr_uv1_4[uv4good], fluxerr_uv1_5[uv5good]]




;fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good], fluxerr_uv1_4[uv4good], fluxerr_uv1_5[uv5good]]


if n_muv eq 1 then wave=[wave_opt1_1[opt1good], wave_uv1_1[uv1good] ]

if n_muv gt 1 then wave=[wave_opt1_1[opt1good], wave_uv1_1[uv1good], wave_uv1_2[uv2good] ]

if n_muv gt 2 then wave=[wave_opt1_1[opt1good], wave_uv1_1[uv1good], wave_uv1_2[uv2good], wave_uv1_3[uv3good]]

if n_muv gt 3 then wave=[wave_opt1_1[opt1good], wave_uv1_1[uv1good], wave_uv1_2[uv2good], wave_uv1_3[uv3good], wave_uv1_4[uv4good]]

if n_muv gt 4 then wave=[wave_opt1_1[opt1good], wave_uv1_1[uv1good], wave_uv1_2[uv2good], wave_uv1_3[uv3good], wave_uv1_4[uv4good], wave_uv1_5[uv4good]]


;wave=[wave_opt1_1[opt1good], wave_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')], wave_uv1_2[uv2good], wave_uv1_3[uv3good], wave_uv1_4[uv4good], wave_uv1_5[uv4good]]

;flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good], flux_uv1_3[uv3good], flux_uv1_4[uv4good], flux_uv1_5[uv5good]]
;fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good], fluxerr_uv1_4[uv4good], fluxerr_uv1_5[uv5good]]


cgplot, wave, flux, color='blue', /ylog, yrange=[min(flux[where(flux gt 0.0)]),max(flux[where(flux gt 0.0)])], psym=-3
oplot, [wavestart,wavestart], [10^(-18.0),10^(-12.0)], psym=-5,color='black'
oplot, [wavestop,wavestop], [10^(-18.0),10^(-12.0)], psym=-5, color='black'

print, "last chance to change wavestart and wavestop if you want"
stop
;;;
output=output+strcompress(uint(wavestart), /REMOVE_ALL)+'_'+strcompress(uint(wavestop), /REMOVE_ALL)+'.dat'
nbins=(wavestop-wavestart)/5+1

newwave=fltarr(nbins)
newflux=make_array(nbins,/double, value=!Values.F_NAN)
newfluxerr=make_array(nbins,/double, value=!Values.F_NAN)


for i=0,n_elements(newwave)-1 do newwave[i]=wavestart+i*5.0

weights=2D
for i=0,n_elements(newwave)-1 do begin

	match=where(abs(wave/(1+z) - newwave[i]) le 5 and finite(flux) eq 1,count)
	if count gt 0 then	weights=1.0/(fluxerr[match]^2.0)
	if count gt 0 then	newflux[i]=(1.0+z)*total( flux[match]*weights ) / total(weights)
	if count gt 0 then	newfluxerr[i]=(1.0+z)*sqrt(1.0/total( weights) )
endfor

empty=where(finite(newflux) eq 0)
good =where(finite(newflux) eq 1)



cgplot, newwave, newflux, color='blue', /ylog, yrange=[min(newflux[where(newflux gt 0.0)]),max(newflux[where(newflux gt 0.0)])], psym=-3
oplot, newwave, newflux, psym=-0, color='black'
;for n=0,n_elements(empty)-1 do newflux[empty[n]]=interpol(newflux[good], newwave[good], newwave[empty[n]])
;for n=0,n_elements(empty)-1 do newfluxerr[empty[n]]=2.0*interpol(newfluxerr[good], newwave[good], newwave[empty[n]])

;cgoplot, newwave[empty], newflux[empty], psym=14, color='red'


stop

print, 'begin first epoch'
for i=0,n_elements(newwave)-1 do print, newwave[i], newflux[i], newfluxerr[i]
print, 'end first epoch'




writecol, output, newwave, newflux, newfluxerr




print, 'all done'
stop

end