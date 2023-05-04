pro SN2021fxy_hst


z=0.0094

;;;;; use error column/ data quality

readcol, 'SN2021fxy_nuv_20210329_1.dat', wave_opt1_1, flux_opt1_1, fluxerr_opt1_1, dq_opt1_1, format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210329_1.dat', wave_uv1_1,  flux_uv1_1,  fluxerr_uv1_1,  dq_uv1_1,  format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210329_2.dat', wave_uv1_2,  flux_uv1_2,  fluxerr_uv1_2,  dq_uv1_2,  format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210329_3.dat', wave_uv1_3,  flux_uv1_3,  fluxerr_uv1_3,  dq_uv1_3,  format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210330_1.dat', wave_uv1_4,  flux_uv1_4,  fluxerr_uv1_4,  dq_uv1_4,  format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210330_2.dat', wave_uv1_5,  flux_uv1_5,  fluxerr_uv1_5,  dq_uv1_5,  format=('F,F,F,I')


dqbits_uv1_1=strarr(16,n_elements(dq_uv1_1))
dqbits_uv1_2=strarr(16,n_elements(dq_uv1_2))
dqbits_uv1_3=strarr(16,n_elements(dq_uv1_3))
dqbits_uv1_4=strarr(16,n_elements(dq_uv1_4))
dqbits_uv1_5=strarr(16,n_elements(dq_uv1_5))
dqbits_opt1_1=strarr(16,n_elements(dq_opt1_1))

for n=0,n_elements(dq_opt1_1)-1 do dqbits_opt1_1[*,n]=binary(dq_opt1_1[n])
for n=0,n_elements(dq_uv1_1)-1 do dqbits_uv1_1[*,n]=binary(dq_uv1_1[n])
for n=0,n_elements(dq_uv1_2)-1 do dqbits_uv1_2[*,n]=binary(dq_uv1_2[n])
for n=0,n_elements(dq_uv1_3)-1 do dqbits_uv1_3[*,n]=binary(dq_uv1_3[n])
for n=0,n_elements(dq_uv1_4)-1 do dqbits_uv1_4[*,n]=binary(dq_uv1_4[n])
for n=0,n_elements(dq_uv1_5)-1 do dqbits_uv1_5[*,n]=binary(dq_uv1_5[n])

cgplot,   wave_opt1_1[where(dq_opt1_1 ge 0)],   flux_opt1_1[where(dq_opt1_1 ge 0)],   yrange=[0,2.5*10.0^(-14)], xrange=[2500,5300]
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=2
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[5,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[5,*] eq '1')]],   psym=16, color='red', symsize=1.5
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[6,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[6,*] eq '1')]],   psym=14, color='green'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[7,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[7,*] eq '1')]],   psym=15, color='orange'

cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] ne '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] ne '1')]], color='blue'

cgoplot,   wave_opt1_1[where(dqbits_opt1_1[0,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[0,*] eq '1')]],   psym=14, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[1,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[1,*] eq '1')]],   psym=15, color='red'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[2,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[2,*] eq '1')]],   psym=16, color='green'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[3,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[3,*] eq '1')]],   psym=14, color='orange'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[4,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[4,*] eq '1')]],   psym=15, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[8,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[8,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[9,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[9,*] eq '1')]],   psym=14, color='red'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[10,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[10,*] eq '1')]],   psym=14, color='green'
;cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=2
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[12,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[12,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[13,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[13,*] eq '1')]],   psym=14, color='black'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[14,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[14,*] eq '1')]],   psym=15, color='black'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[15,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[15,*] eq '1')]],   psym=16, color='black'



;;;;;;;;  this gives a good spectrum
cgplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )],   flux_opt1_1[[where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )]], psym=14, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1' or dqbits_opt1_1[6,*] eq '1' )],   flux_opt1_1[where(dqbits_opt1_1[11,*] eq '1' or dqbits_opt1_1[6,*] eq '1' )], psym=15, color='red'

cgoplot,   wave_opt1_1,   flux_opt1_1







cgplot,   wave_uv1_1,   flux_uv1_1, color='blue', yrange=[-1.0*10.0^(-15.0),1.2*10.0^(-14.0)]



print, wave_uv1_1[where(dqbits_uv1_1[5] eq '1')]

;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[0,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[0,*]  eq '1')]],   psym=14, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[1,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[1,*]  eq '1')]],   psym=15, color='red'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[2,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[2,*]  eq '1')]],   psym=16, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[3,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[3,*]  eq '1')]],   psym=14, color='orange'
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[4,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[4,*]  eq '1')]],   psym=15, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[5,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[5,*]  eq '1')]],   psym=16, color=red, symsize=2
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[6,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[6,*]  eq '1')]],   psym=14, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[7,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[7,*]  eq '1')]],   psym=15, color='orange'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[8,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[8,*]  eq '1')]],   psym=16, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[9,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[9,*]  eq '1')]],   psym=14, color='red'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[10,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[10,*] eq '1')]],   psym=14, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[11,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=3
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[12,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[12,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[13,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[13,*] eq '1')]],   psym=14, color='black'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[14,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[14,*] eq '1')]],   psym=15, color='black'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[15,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[15,*] eq '1')]],   psym=16, color='black'

print, size(where(dqbits_uv1_1[0,*] eq '1'))
print, size(where(dqbits_uv1_1[1,*] eq '1'))
print, size(where(dqbits_uv1_1[2,*] eq '1'))
print, size(where(dqbits_uv1_1[3,*] eq '1'))
print, size(where(dqbits_uv1_1[4,*] eq '1'))
print, size(where(dqbits_uv1_1[5,*] eq '1'))
print, size(where(dqbits_uv1_1[6,*] eq '1'))
print, size(where(dqbits_uv1_1[7,*] eq '1'))
print, size(where(dqbits_uv1_1[8,*] eq '1'))
print, size(where(dqbits_uv1_1[9,*] eq '1'))
print, size(where(dqbits_uv1_1[10,*] eq '1'))
print, size(where(dqbits_uv1_1[11,*] eq '1'))
print, size(where(dqbits_uv1_1[12,*] eq '1'))
print, size(where(dqbits_uv1_1[13,*] eq '1'))
print, size(where(dqbits_uv1_1[14,*] eq '1'))
print, size(where(dqbits_uv1_1[15,*] eq '1'))


print, where(dqbits_uv1_1[4,*] eq '1')
print, where(dqbits_uv1_1[6,*] eq '1')
print, where(dqbits_uv1_1[13,*] eq '1')

cgoplot,   wave_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')],    flux_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')],   psym=16, color='green'


print, size(where(dqbits_uv1_2[0,*] eq '1'))
print, size(where(dqbits_uv1_2[1,*] eq '1'))
print, size(where(dqbits_uv1_2[2,*] eq '1'))
print, size(where(dqbits_uv1_2[3,*] eq '1'))
print, size(where(dqbits_uv1_2[4,*] eq '1'))
print, size(where(dqbits_uv1_2[5,*] eq '1'))
print, size(where(dqbits_uv1_2[6,*] eq '1'))
print, size(where(dqbits_uv1_2[7,*] eq '1'))
print, size(where(dqbits_uv1_2[8,*] eq '1'))
print, size(where(dqbits_uv1_2[9,*] eq '1'))
print, size(where(dqbits_uv1_2[10,*] eq '1'))
print, size(where(dqbits_uv1_2[11,*] eq '1'))
print, size(where(dqbits_uv1_2[12,*] eq '1'))
print, size(where(dqbits_uv1_2[13,*] eq '1'))
print, size(where(dqbits_uv1_2[14,*] eq '1'))
print, size(where(dqbits_uv1_2[15,*] eq '1'))





print, size(where(dqbits_uv1_5[0,*] eq '1'))
print, size(where(dqbits_uv1_5[1,*] eq '1'))
print, size(where(dqbits_uv1_5[2,*] eq '1'))
print, size(where(dqbits_uv1_5[3,*] eq '1'))
print, size(where(dqbits_uv1_5[4,*] eq '1'))
print, size(where(dqbits_uv1_5[5,*] eq '1'))
print, size(where(dqbits_uv1_5[6,*] eq '1'))
print, size(where(dqbits_uv1_5[7,*] eq '1'))
print, size(where(dqbits_uv1_5[8,*] eq '1'))
print, size(where(dqbits_uv1_5[9,*] eq '1'))
print, size(where(dqbits_uv1_5[10,*] eq '1'))
print, size(where(dqbits_uv1_5[11,*] eq '1'))
print, size(where(dqbits_uv1_5[12,*] eq '1'))
print, size(where(dqbits_uv1_5[13,*] eq '1'))
print, size(where(dqbits_uv1_5[14,*] eq '1'))
print, size(where(dqbits_uv1_5[15,*] eq '1'))




cgplot,   wave_uv1_1,   flux_uv1_1/fluxerr_uv1_1, xrange=[1600,2000]




 mg2= [2796.0,2803] 
;; weak mg2 from MW, nothing obvious between z=0.02-0.03



;1600 - 5680 combining both spectra


uv1good=[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')]
uv2good=[where(dqbits_uv1_2[4,*] ne '1' and dqbits_uv1_2[6,*] ne '1' and dqbits_uv1_2[13,*] ne '1')]
uv3good=[where(dqbits_uv1_3[4,*] ne '1' and dqbits_uv1_3[6,*] ne '1' and dqbits_uv1_3[13,*] ne '1')]
uv4good=[where(dqbits_uv1_4[4,*] ne '1' and dqbits_uv1_4[6,*] ne '1' and dqbits_uv1_4[13,*] ne '1')]
uv5good=[where(dqbits_uv1_5[4,*] ne '1' and dqbits_uv1_5[6,*] ne '1' and dqbits_uv1_5[13,*] ne '1')]

opt1good=where(dqbits_opt1_1[11,*] ne '1' and dqbits_opt1_1[6,*] ne '1' )


wave=[wave_opt1_1[opt1good], wave_uv1_1[where(dqbits_uv1_1[4,*] ne '1' and dqbits_uv1_1[6,*] ne '1' and dqbits_uv1_1[13,*] ne '1')], wave_uv1_2[where(dqbits_uv1_2[4,*] ne '1' and dqbits_uv1_2[6,*] ne '1' and dqbits_uv1_2[13,*] ne '1')], wave_uv1_3[where(dqbits_uv1_3[4,*] ne '1' and dqbits_uv1_3[6,*] ne '1' and dqbits_uv1_3[13,*] ne '1')], wave_uv1_4[where(dqbits_uv1_4[4,*] ne '1' and dqbits_uv1_4[6,*] ne '1' and dqbits_uv1_4[13,*] ne '1')], wave_uv1_5[where(dqbits_uv1_5[4,*] ne '1' and dqbits_uv1_5[6,*] ne '1' and dqbits_uv1_5[13,*] ne '1')]]

flux=[flux_opt1_1[opt1good], flux_uv1_1[uv1good], flux_uv1_2[uv2good], flux_uv1_3[uv3good], flux_uv1_4[uv4good], flux_uv1_5[uv5good]]
fluxerr=[fluxerr_opt1_1[opt1good], fluxerr_uv1_1[uv1good], fluxerr_uv1_2[uv2good], fluxerr_uv1_3[uv3good], fluxerr_uv1_4[uv4good], fluxerr_uv1_5[uv5good]]




;;;
newwave=fltarr(800+1)
newflux=make_array(800+1,/double, value=!Values.F_NAN)
newfluxerr=make_array(800+1,/double, value=!Values.F_NAN)


for i=0,n_elements(newwave)-1 do newwave[i]=1600+i*5.0

weights=2D
for i=0,n_elements(newwave)-1 do begin

	match=where(abs(wave/(1+z) - newwave[i]) le 5 and finite(flux) eq 1,count)
	if count gt 0 then	weights=1.0/(fluxerr[match]^2.0)
	if count gt 0 then	newflux[i]=(1.0+z)*total( flux[match]*weights ) / total(weights)
	if count gt 0 then	newfluxerr[i]=(1.0+z)*sqrt(1.0/total( weights) )
endfor

empty=where(finite(newflux) eq 0)
good=where(finite(newflux) eq 1)



cgplot, newwave, newflux, color='blue'

for n=0,n_elements(empty)-1 do newflux[empty[n]]=interpol(newflux[good], newwave[good], newwave[empty[n]])
for n=0,n_elements(empty)-1 do newfluxerr[empty[n]]=2.0*interpol(newfluxerr[good], newwave[good], newwave[empty[n]])

cgoplot, newwave[empty], newflux[empty], psym=14, color='red'




print, 'begin first epoch'
for i=0,n_elements(newwave)-1 do print, newwave[i], newflux[i], newfluxerr[i]
print, 'end first epoch'


stop







;;;;; epoch 2

readcol, 'SN2021fxy_nuv_20210401_1.dat', wave_opt2_1, flux_opt2_1, fluxerr_opt2_1, dq_opt2_1, format=('F,F,F,I')

readcol, 'SN2021fxy_muv_20210401_1.dat', wave_uv2_1, flux_uv2_1, fluxerr_uv2_1, dq_uv2_1, format=('F,F,F,I')


;;;;;;;;;;;;;;


dqbits_uv2_1=intarr(16,n_elements(dq_uv2_1))
dqbits_opt2_1=intarr(16,n_elements(dq_opt2_1))

for n=0,n_elements(dq_opt2_1)-1 do dqbits_opt2_1[*,n]=binary(dq_opt2_1[n])
for n=0,n_elements(dq_uv2_1)-1  do dqbits_uv2_1[*,n]=binary(dq_uv2_1[n])

cgplot,   wave_opt1_1[where(dq_opt1_1 ge 0)],   flux_opt1_1[where(dq_opt1_1 ge 0)],   yrange=[0,2.5*10.0^(-14)], xrange=[2500,5300]
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=2
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[5,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[5,*] eq '1')]],   psym=16, color='red', symsize=1.5
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[6,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[6,*] eq '1')]],   psym=14, color='green'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[7,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[7,*] eq '1')]],   psym=15, color='orange'

cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] ne '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] ne '1')]], color='blue'

cgoplot,   wave_opt1_1[where(dqbits_opt1_1[0,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[0,*] eq '1')]],   psym=14, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[1,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[1,*] eq '1')]],   psym=15, color='red'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[2,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[2,*] eq '1')]],   psym=16, color='green'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[3,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[3,*] eq '1')]],   psym=14, color='orange'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[4,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[4,*] eq '1')]],   psym=15, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[8,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[8,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[9,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[9,*] eq '1')]],   psym=14, color='red'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[10,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[10,*] eq '1')]],   psym=14, color='green'
;cgoplot,   wave_opt1_1[where(dqbits_opt1_1[11,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=2
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[12,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[12,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[13,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[13,*] eq '1')]],   psym=14, color='black'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[14,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[14,*] eq '1')]],   psym=15, color='black'
cgoplot,   wave_opt1_1[where(dqbits_opt1_1[15,*] eq '1')],   flux_opt1_1[[where(dqbits_opt1_1[15,*] eq '1')]],   psym=16, color='black'



;;;;;;;;  this gives a good spectrum
cgplot,   wave_opt2_1[where(dqbits_opt2_1[11,*] ne '1' and dqbits_opt2_1[6,*] ne '1' )],   flux_opt2_1[[where(dqbits_opt2_1[11,*] ne '1' and dqbits_opt2_1[6,*] ne '1' )]], psym=14, color='blue'
cgoplot,   wave_opt2_1[where(dqbits_opt2_1[11,*] eq '1' or dqbits_opt2_1[6,*] eq '1' )],   flux_opt2_1[where(dqbits_opt2_1[11,*] eq '1' or dqbits_opt2_1[6,*] eq '1' )], psym=15, color='red'

cgoplot,   wave_opt2_1,   flux_opt2_1







cgplot,   wave_uv1_1,   flux_uv1_1, color='blue', yrange=[-1.0*10.0^(-15.0),1.2*10.0^(-14.0)]



print, wave_uv1_1[where(dqbits_uv1_1[5] eq '1')]

;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[0,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[0,*]  eq '1')]],   psym=14, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[1,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[1,*]  eq '1')]],   psym=15, color='red'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[2,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[2,*]  eq '1')]],   psym=16, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[3,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[3,*]  eq '1')]],   psym=14, color='orange'
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[4,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[4,*]  eq '1')]],   psym=15, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[5,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[5,*]  eq '1')]],   psym=16, color=red, symsize=2
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[6,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[6,*]  eq '1')]],   psym=14, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[7,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[7,*]  eq '1')]],   psym=15, color='orange'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[8,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[8,*]  eq '1')]],   psym=16, color='blue'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[9,*] eq '1')],    flux_uv1_1[[where(dqbits_uv1_1[9,*]  eq '1')]],   psym=14, color='red'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[10,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[10,*] eq '1')]],   psym=14, color='green'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[11,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[11,*] eq '1')]],   psym=15, color='orange', symsize=3
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[12,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[12,*] eq '1')]],   psym=16, color='blue'
cgoplot,   wave_uv1_1[where(dqbits_uv1_1[13,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[13,*] eq '1')]],   psym=14, color='black'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[14,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[14,*] eq '1')]],   psym=15, color='black'
;cgoplot,   wave_uv1_1[where(dqbits_uv1_1[15,*] eq '1')],   flux_uv1_1[[where(dqbits_uv1_1[15,*] eq '1')]],   psym=16, color='black'


;;;;;;;;;;;;;;

cgplot,  wave_uv2_1[where(dq_uv2_1 eq 0)], flux_uv2_1[where(dq_uv2_1 eq 0)], xrange=[1600,5800], yrange=[0,2.0*10.0^(-14)], psym=4
cgoplot, wave_opt2_1[where(dq_opt2_1 eq 0)], flux_opt2_1[where(dq_opt2_1 eq 0)], color='red', psym=4

cgplot,  wave_uv2_1[where(dq_uv2_1 eq 0)], flux_uv2_1[where(dq_uv2_1 eq 0)], xrange=[1600,3200], yrange=[0,2.0*10.0^(-15)]
cgoplot, wave_opt2_1[where(dq_opt2_1 eq 0)], flux_opt2_1[where(dq_opt2_1 eq 0)], color='red'

wave=[wave_opt2_1, wave_uv2_1]
flux=[flux_opt2_1, flux_uv2_1]
fluxerr=[fluxerr_opt2_1, fluxerr_uv2_1]
dq=[dq_opt2_1, dq_uv2_1]


newwave2=fltarr(800+1)
newflux2=make_array(800+1,/double, value=!Values.F_NAN)
newfluxerr2=make_array(800+1,/double, value=!Values.F_NAN)


for i=0,n_elements(newwave2)-1 do newwave2[i]=1600+i*5.0

weights=2D
for i=0,n_elements(newwave)-1 do begin

	match=where(abs(wave - newwave2[i]) le 5 and dq eq 0 and finite(flux) eq 1,count)
	if count gt 0 then	weights=1.0/(fluxerr[match]^2.0)
	if count gt 0 then	newflux2[i]=total( flux[match]*weights ) / total(weights)
	if count gt 0 then	newfluxerr2[i]=sqrt(1.0/total( weights) )

endfor

print, 'begin epoch 2'
for i=0,n_elements(newwave2)-1 do print, newwave2[i], newflux2[i], newfluxerr2[i]
print, 'end epoch 2'



;;;;;; 

readcol, '$SNFOLDER/UVSPECTRA/hst/SN2011fe/SN2011fe_sullivan/ptf11kly_20110913.obs.dat', fewave, feflux, feerr
readcol, '$SNFOLDER/github/SN2017erp/hstspec/SN2017erp_hst_20170629.dat', erpwave, erpflux, erperr               
readcol, 'SN2021fxy_hst_20210329.dat', fxywave, fxyflux, fxyerr                                                  
cgplot, fxywave, fxyflux        

cgoplot, fewave, feflux*3.0*10.0^(-2), color='blue'                                                           
cgoplot, fewave, feflux*2.5*10.0^(-2), color='blue'
cgoplot, erpwave, erpflux*0.8, color='red'         
cgoplot, erpwave, erpflux*0.7, color='red'
cgplot, fxywave, fxyflux                           
cgoplot, erpwave, erpflux*0.7, color='red'
cgoplot, fewave, feflux*2.5*10.0^(-2), color='blue'
cgplot, fxywave, fxyflux, xrange=[1600,3600]       
cgoplot, fewave, feflux*2.5*10.0^(-2), color='blue'
cgoplot, erpwave, erpflux*0.7, color='red'         
cgplot, fxywave, fxyflux, xrange=[1600,3000]       
cgoplot, erpwave, erpflux*0.7, color='red'  
cgplot, fxywave, fxyflux, xrange=[1600,2600]
cgoplot, erpwave, erpflux*0.7, color='red'  
cgoplot, fewave, feflux*2.5*10.0^(-3), color='blue'


;;;;; use error column/ data quality

readcol, 'SN2021fxy_nuv_20210408_1.dat', wave_opt4_1, flux_opt4_1, fluxerr_opt4_1, dq_opt4_1

readcol, 'SN2021fxy_muv_20210408_1.dat', wave_uv4_1, flux_uv4_1, fluxerr_uv4_1, dq_uv4_1

readcol, 'SN2021fxy_muv_20210408_2.dat', wave_uv4_2, flux_uv4_2, fluxerr_uv4_2, dq_uv4_2

readcol, 'SN2021fxy_muv_20210408_3.dat', wave_uv4_3, flux_uv4_3, fluxerr_uv4_3, dq_uv4_3




cgplot,  wave_uv4_1[where(dq_uv4_1 eq 0)], flux_uv1_1[where(dq_uv4_1 eq 0)], xrange=[1600,5800], yrange=[0,2.0*10.0^(-14)]
cgoplot, wave_uv4_2[where(dq_uv4_2 eq 0)], flux_uv1_2[where(dq_uv4_2 eq 0)], color='red'
cgoplot, wave_opt4_1[where(dq_opt4_1 eq 0)], flux_opt1_1[where(dq_opt4_1 eq 0)], color='blue'

cgplot,  wave_uv4_1[where(dq_uv4_1 eq 0)], flux_uv4_1[where(dq_uv4_1 eq 0)], xrange=[1600,3200], yrange=[0,4.0*10.0^(-16)]
cgoplot, wave_uv4_2[where(dq_uv4_2 eq 0)], flux_uv4_2[where(dq_uv4_2 eq 0)], color='red'
cgoplot, wave_uv4_3[where(dq_uv4_3 eq 0)], flux_uv4_3[where(dq_uv4_3 eq 0)], color='orange'
cgoplot, wave_opt4_1[where(dq_opt4_1 eq 0)], flux_opt4_1[where(dq_opt4_1 eq 0)], color='blue'

;1600 - 5680 combining both spectra


wave=[wave_opt4_1, wave_uv4_1, wave_uv4_2, wave_uv4_3]
flux=[flux_opt4_1, flux_uv4_1, flux_uv4_2, flux_uv4_3]
fluxerr=[fluxerr_opt4_1, fluxerr_uv4_1, fluxerr_uv4_2, fluxerr_uv4_3]
dq=[dq_opt4_1, dq_uv4_1, dq_uv4_2, dq_uv4_3]


;;;
newwave4=fltarr(800+1)
newflux4=make_array(800+1,/double, value=!Values.F_NAN)
newfluxerr4=make_array(800+1,/double, value=!Values.F_NAN)


for i=0,n_elements(newwave4)-1 do newwave4[i]=1600+i*5.0

weights=2D
for i=0,n_elements(newwave)-1 do begin

	match=where(abs(wave - newwave4[i]) le 5 and dq eq 0 and finite(flux) eq 1,count)
	if count gt 0 then	weights=1.0/(fluxerr[match]^2.0)
	if count gt 0 then	newflux4[i]=total( flux[match]*weights ) / total(weights)
	if count gt 0 then	newfluxerr4[i]=sqrt(1.0/total( weights) )

endfor

print, 'begin epoch 4'
for i=0,n_elements(newwave4)-1 do print, newwave4[i], newflux4[i], newfluxerr4[i]
print, 'end epoch 4'


cgplot, newwave, newflux              
cgoplot, newwave2, newflux2, color='red'
cgoplot, newwave4, newflux4, color='blue'




readcol, 'SN2021fxy_hst_20210401.dat', fxywave, fxyflux, fxyerr                                                  
cgplot, fxywave, fxyflux, color='blue', xrange=[3700,3800]

readcol, '~/Desktop/SN/github/UVSN/original_spectra/21fxy_210401_apo_full_trim_spex_smooth.ascii', optwave, optflux, opterr                                                  
cgoplot, optwave*(1.0+0.009393), optflux*(mean(fxyflux[where(fxywave gt  4000 and fxywave lt 5000)],/nan)/mean(optflux[where(optwave gt  4000 and optwave lt 5000)],/nan)), color='red'



readcol, 'SN2021fxy_nuv_20210401_1.dat', fxywave, fxyflux, fxyerr                                                  
cgoplot, fxywave, fxyflux, color='blue'


print, 'final stop'
stop
end