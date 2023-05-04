pro SN2016ccj_hst




;;;;; use error column/ data quality


readcol, 'SN2016ccj_opt_20160514_1.dat', wave_opt1_1, flux_opt1_1, fluxerr_opt1_1, dq_opt1_1

readcol, 'SN2016ccj_uv_20160514_1.dat', wave_uv1_1, flux_uv1_1, fluxerr_uv1_1, dq_uv1_1

readcol, 'SN2016ccj_uv_20160514_2.dat', wave_uv1_2, flux_uv1_2, fluxerr_uv1_2, dq_uv1_2




cgplot, wave_uv1_1[where(dq_uv1_1 eq 0)], flux_uv1_1[where(dq_uv1_1 eq 0)], xrange=[2800,3400], yrange=[0,2.0*10.0^(-15)]
cgoplot, wave_uv1_2[where(dq_uv1_2 eq 0)], flux_uv1_2[where(dq_uv1_2 eq 0)], color='red'
cgoplot, wave_opt1_1[where(dq_opt1_1 eq 0)], flux_opt1_1[where(dq_opt1_1 eq 0)], color='blue'



plot, wave_uv1_1[where(dq_uv1_1 eq 0)], flux_uv1_1[where(dq_uv1_1 eq 0)], xrange=[1600,8000], yrange=[0,2.0*10.0^(-15)]
oplot, wave_uv1_2[where(dq_uv1_2 eq 0)], flux_uv1_2[where(dq_uv1_2 eq 0)]
oplot, wave_opt1_1[where(dq_opt1_1 eq 0)], flux_opt1_1[where(dq_opt1_1 eq 0)]



plot, wave_uv1_2[where(dq_uv1_2 eq 0)], flux_uv1_2[where(dq_uv1_2 eq 0)], xrange=[1600,2100], yrange=[0,2.0*10.0^(-15)]



plot, wave_uv1_1[where(dq_uv1_1 eq 0)], flux_uv1_1[where(dq_uv1_1 eq 0)], xrange=[1600,4000], yrange=[0,2.0*10.0^(-15)]
oplot, wave_uv1_2[where(dq_uv1_2 eq 0)], flux_uv1_2[where(dq_uv1_2 eq 0)]
oplot, wave_opt1_1[where(dq_opt1_1 eq 0)], flux_opt1_1[where(dq_opt1_1 eq 0)]





readcol, 'SN2016ccj_opt_20160521_1.dat', wave_opt2_1, flux_opt2_1, fluxerr_opt2_1, dq_opt2_1

readcol, 'SN2016ccj_uv_20160521_1.dat', wave_uv2_1, flux_uv2_1, fluxerr_uv2_1, dq_uv2_1

readcol, 'SN2016ccj_uv_20160521_2.dat', wave_uv2_2, flux_uv2_2, fluxerr_uv2_2, dq_uv2_2


cgoplot, wave_opt2_1[where(dq_opt2_1 eq 0)], flux_opt2_1[where(dq_opt2_1 eq 0)], color='blue'
cgoplot, wave_uv2_2[where(dq_uv2_1 eq 0)], flux_uv2_1[where(dq_uv2_1 eq 0)] , color='blue'
cgoplot, wave_uv2_2[where(dq_uv2_2 eq 0)], flux_uv2_2[where(dq_uv2_2 eq 0)] , color='blue'


 mg2= [2796.0,2803] 
;; weak mg2 from MW, nothing obvious between z=0.02-0.03



readcol, 'SN2016ccj_opt_20160529_1.dat', wave_opt3_1, flux_opt3_1, fluxerr_opt3_1, dq_opt3_1

readcol, 'SN2016ccj_uv_20160529_1.dat', wave_uv3_1, flux_uv3_1, fluxerr_uv3_1, dq_uv3_1

readcol, 'SN2016ccj_uv_20160529_2.dat', wave_uv3_2, flux_uv3_2, fluxerr_uv3_2, dq_uv3_2

readcol, 'SN2016ccj_uv_20160529_3.dat', wave_uv3_3, flux_uv3_3, fluxerr_uv3_3, dq_uv3_3


cgoplot, wave_opt3_1[where(dq_opt3_1 eq 0)], flux_opt3_1[where(dq_opt3_1 eq 0)], color='red'
cgoplot, wave_uv3_1[where(dq_uv3_1 eq 0)], flux_uv3_1[where(dq_uv3_1 eq 0)] , color='red'
cgoplot, wave_uv3_2[where(dq_uv3_2 eq 0)], flux_uv3_2[where(dq_uv3_2 eq 0)] , color='red'
cgoplot, wave_uv3_3[where(dq_uv3_3 eq 0)], flux_uv3_3[where(dq_uv3_3 eq 0)] , color='red'

;1600 - 5680 combining both spectra

wave=[wave_opt3_1, wave_uv3_1, wave_uv3_2, wave_uv3_3]
flux=[flux_opt3_1, flux_uv3_1, flux_uv3_2, flux_uv3_3]
fluxerr=[fluxerr_opt3_1, fluxerr_uv3_1, fluxerr_uv3_2, fluxerr_uv3_3]
dq=[dq_opt3_1, dq_uv3_1, dq_uv3_2, dq_uv3_3]


wave=[wave_opt1_1, wave_uv1_1, wave_uv1_2]
flux=[flux_opt1_1, flux_uv1_1, flux_uv1_2]
fluxerr=[fluxerr_opt1_1, fluxerr_uv1_1, fluxerr_uv1_2]
dq=[dq_opt1_1, dq_uv1_1, dq_uv1_2]

wave=[wave_opt2_1, wave_uv2_1, wave_uv2_2]
flux=[flux_opt2_1, flux_uv2_1, flux_uv2_2]
fluxerr=[fluxerr_opt2_1, fluxerr_uv2_1, fluxerr_uv2_2]
dq=[dq_opt2_1, dq_uv2_1, dq_uv2_2]


;;;
newwave=fltarr(800+1)
newflux=make_array(800+1,/double, value=!Values.F_NAN)
newfluxerr=make_array(800+1,/double, value=!Values.F_NAN)


for i=0,n_elements(newwave)-1 do newwave[i]=1600+i*5.0

weights=2D
for i=0,n_elements(newwave)-1 do begin

	match=where(abs(wave - newwave[i]) le 5 and dq eq 0 and finite(flux) eq 1,count)
	if count gt 0 then	weights=1.0/(fluxerr[match]^2.0)
	if count gt 0 then	newflux[i]=total( flux[match]*weights ) / total(weights)
	if count gt 0 then	newfluxerr[i]=sqrt(1.0/total( weights) )

endfor


for i=0,n_elements(newwave)-1 do print, newwave[i], newflux[i], newfluxerr[i]




print, 'final stop'
stop
end