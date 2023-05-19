pro SN2017erp_epochs


readcol, 'SN2017erp_epochs.dat', tel, mjdstart, mjdend, lambdastart, lambdaend, format='A,F,F'

mjdmid=(mjdstart+mjdend)/2.0d
jdmid=2400000.5d+mjdmid

datecol=tel
for i=0,n_elements(datecol)-1 do datecol[i]=DATE_CONV( JDmid[i], 'S' )



specepochs=mjdmid-57934.9d


for i=0,n_elements(tel)-1 do print, specepochs[i], ' & ', datecol[i], ' & ',  tel[i], ' & ', lambdastart[i], ' & ', lambdaend[i], ' \\', format='(F5.1, A, A, A, A, A, A, A, I, A, I, A)'