filename(n) = sprintf('dens_%i.dat', n) 
 set xlabel 'Distance [AU]' 
 set ylabel 'Temperature [K]' 
 set title 'Temperature profile' 
 plot filename(i) using 1:3 title sprintf('t=%i nap',i) with line  
 i=i+10 
 if (i < n) reread 
 