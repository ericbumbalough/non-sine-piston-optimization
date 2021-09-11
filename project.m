function [ out ] = myfun( vTDC, CR, mass)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%% List of Inputs
%vTDC: specific volume at TDCf [m^3/kg]
%CR: compression ratio [1]
%mass: trapped mass [kg]

%% List of Outputs

%% Parameters
res = .1; %CA resolution [degree]

%% stuffs
vBDC = CR * vTDC;
CA = 0:res:180;


v=spline([0 180],[0 [vTDC vBDC] 0], CA);
plot(CA,v)
 

end
%%EOF