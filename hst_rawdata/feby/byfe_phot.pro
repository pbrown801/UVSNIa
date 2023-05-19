pro byfe_phot

restore, filename='SN2011fe_bolmags.sav'



snphot_array, '~/Desktop/SN/SwiftSNarchive/PhotArchive/SN2011by_uvotB13.dat', 'by.fits', dt=by
bybpeak=by.bbtime(where(min(by.bbmags) eq by.bbmags))


plot, SN2011fe_bolepochs, SN2011fe_bolmags[1,*]-SN2011fe_bolmags[2,*], psym=-4
bym2w1=where( finite(by.mag_array[2,*]) eq 1 and finite(by.mag_array[5,*]) eq 1 )
oplot, by.time_array[bym2w1]-bybpeak[0], by.mag_array[1,bym2w1]-by.mag_array[2,bym2w1], psym=-5

plot, SN2011fe_bolepochs, SN2011fe_bolmags[2,*]-SN2011fe_bolmags[5,*], psym=-4
byw1v=where( finite(by.mag_array[2,*]) eq 1 and finite(by.mag_array[5,*]) eq 1 )
oplot, by.time_array[byw1v]-bybpeak[0], by.mag_array[2,byw1v]-by.mag_array[5,byw1v], psym=-5

plot, SN2011fe_bolepochs, SN2011fe_bolmags[2,*]-SN2011fe_bolmags[4,*], psym=-4
byw1b=where( finite(by.mag_array[2,*]) eq 1 and finite(by.mag_array[4,*]) eq 1 )
oplot, by.time_array[byw1b]-bybpeak[0], by.mag_array[2,byw1b]-by.mag_array[4,byw1b], psym=-5



plot, SN2011fe_bolepochs, SN2011fe_bolmags[4,*]-SN2011fe_bolmags[5,*], psym=-4
bybv=where( finite(by.mag_array[4,*]) eq 1 and finite(by.mag_array[5,*]) eq 1 )
oplot, by.time_array[bybv]-bybpeak[0], by.mag_array[4,bybv]-by.mag_array[5,bybv], psym=-5


stop
end
