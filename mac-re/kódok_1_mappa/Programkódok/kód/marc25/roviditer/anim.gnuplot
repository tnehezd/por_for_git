filename(n) = sprintf('dens_%i.dat', n) 
 set xlabel 'Distance [AU]' 
 set ylabel 'Surface density [M_Sun/AU/AU]' 
 set title 'Surface density profile' 
 plot filename(i) using 1:2 title sprintf('t=%i nap',i) with line  
 i=i+10 
 if (i < n) reread 
 