function T = minburntemp(phi)
%returns the minimum temperature for successful combustion in CELCIUS based
%on a given phi. Based on article by Lavoie, Martz, Wooldridge and Assanis
%in Combustion and Flame Vol 157, Issue 6, pps 1106-1110.

%only works on scalars currently.

if phi > 1
    T = 126.9;
elseif phi >= 0.35
    T = 1596.2 * phi^2 - 3116.3 * phi + 1622;
elseif phi >=.2
    T = 726.8;
elseif phi >= .175
    T = -4000 * phi + 1526.8;
else
    T = 826.8;
end    

end