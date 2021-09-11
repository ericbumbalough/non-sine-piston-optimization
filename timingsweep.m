function [timing,Wout] = timingsweep(start,stop,step,mass,Rtherm)
%sweeps CA0 from start to end and returns timing versus indicated work

timing = start:step:stop;

Wout = zeros(1,length(timing));

for n = 1:length(timing)
    Wout(n) = myfun(mass, 6.25, timing(n), Rtherm);
end

end
