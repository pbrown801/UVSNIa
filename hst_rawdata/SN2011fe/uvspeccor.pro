pro uvspeccor

;;; failed attempt to correct for the red leak

readcol,"SN2011fe_110913_full.dat",speclambda,specflux,/silent

readcol,"g230lb_sens.dat",g230lblambda,g230lbsens,/silent


interpspec=interpol(specflux,speclambda,g230lblambda)
countspec=interpspec*g230lbsens

range=where(g230lblambda ge 2000 and g230lblambda le 10000)
C0=5523.0*total(countspec[range]/g230lblambda[range]^3.0 *10.0)
print, C0


;  try redoing in other units


interpg230lbsens=interpol(g230lbsens,g230lblambda,speclambda)
countspec=interpg230lbsens*specflux

range=where(speclambda ge 2000 and speclambda le 10000)
C0=5523.0*total(countspec[range]/speclambda[range]^3.0 *1.00)
;  assuming the spacing is one pixel
print, C0


plot, g230lblambda, interpspec
oplot, speclambda, specflux


readcol,"SN2011fe_uv_20110913_net.dat",specnet_lambda,specnet_counts,/silent
readcol,"SN2011fe_uv_20110913_1.dat",specnet_lambda,specflux,/silent

netinterpg230lbsens=interpol(g230lbsens,g230lblambda,specnet_lambda)

specnet_flux=specnet_counts/netinterpg230lbsens
sens=specnet_counts/specflux


plot, specnet_lambda, specnet_counts
oplot, speclambda, countspec

C0=0.025
specnet_newcounts=fltarr(n_elements(specnet_lambda))
for x=0,n_elements(specnet_lambda)-1 do specnet_newcounts[x]=specnet_counts[x]-C0*(1+0.0024*x)

specnet_newflux=specnet_newcounts/sens


plot, specnet_lambda, specnet_counts, xrange=[1600,2600]
oplot, specnet_lambda, specnet_newcounts

plot, specnet_lambda, specnet_flux
oplot, specnet_lambda, specnet_newflux

stop
end