pro tweak

readcol,"SN2011fe_uv_20110913_1.dat",uvlambda,uvflux,/silent
readcol,"SN2011fe_fuv_20110913_1.dat",fuvlambda,fuvflux,/silent


plot, uvlambda, uvflux
oplot, fuvlambda, fuvflux

npoints=n_elements(uvlambda)
newuvflux=fltarr(npoints)
for n=0,npoints-1 do newuvflux[n]= (uvflux[n]*(1.0-(1.0-2.0e-15/uvflux[n])*(npoints-n)/npoints))


oplot, uvlambda, newuvflux

for n=0,npoints-1 do print, uvlambda[n], newuvflux[n]

stop
end
