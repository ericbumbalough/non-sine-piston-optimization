function [c] = nonlcon( in )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% mass = in(1);
% splinefit = in(2);
% CA0 = in(3);
% burndur = in(4);
phi = in(5);

c = [0,0,0];

[~, ~, Pmax, RI, minflameT] = engine(in);

c(1) = Pmax - 1.2e7;
c(2) = RI - 5e6;
c(3) = minburntemp(phi) - minflameT;

end

