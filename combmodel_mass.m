function [ Wtot Pmax maxPRR TBDC] = myfun(mass, splinefit, CA0)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%% List of Inputs
%mass: trapped mass [kg]
%splinefit: instantanious compression ratios at CA in CAforce
%CA0: CA at start of burn [degree]

%% List of Outputs
%Wtot: indicated work without pumping loss [W]
%Pmax: maximum cylinder pressure [Pa]
%maxPRR: maximum pressure rise rate [Pa/degree]

%% Parameters
res = 1; %CA resolution [degree]
CR = 12.5; %compression ratio
displacement = 0.055; %engine displacement [m^3]
CAforce = [90]; %CA points that the spline fits at given vforce points
burndur = 30; %Wiebe burn duration [degrees]
lhv = 160000; %lower heating value of fuel [J]
Pint = 100000; %intake pressure [Pa]
gamma = 1.4; %ratio of specific heats. Assumed constant for non firing conditions.
R = 286.9; %ideal gas constant for air. [J/kg/C]
bore = 0.0085; %cylinder bore [m]
RPM = 2000;
Twall = 137; %this was from a Babajimapolis paper. find it. [C]
hmod = 1; %heat transfer multiplier

%% Geometry
VTDC = displacement /(CR - 1); % [m^3]
VBDC = VTDC * CR; % [m^3]
Vspline = splinefit * VTDC; %[m^3]

stoke = 4 * displacement / pi / bore^2; %[m]
S = stroke * 2 * RPM / 60; %mean piston velocity [m/s]
Apiston = pi * bore^2 / 4; %area of piston, assuming a flat piston [m^2]
Ahead = pi * bore^2 / 4; %area of head, assuming a flat piston [m^2]
dt = res / RPM / 30; %time step [s]
Area = Apiston + Ahead + V .* 4 / bore; %area used for heat transfer calculations [m^2]

CA = 0:res:180; %crank angle array
V=spline([0 CAforce 180],[0 [VTDC Vspline VBDC] 0], CA); %fit volume as a spline with 0 slope endpoints
%perhaps i could constrain the spline to have only 1 inflection point?
%constrain with a maximum slope
%constrain piston maximum acceleration
%constrain piston maximum velocity

CAmirror = zeros(1,20/res-1);
Vmirror = zeros(1,20/res-1);

for n=1:20/res-1 %mirror V and CA for bTDC.
    CAmirror(n) = -CA(20/res-n+1);
    Vmirror(n) = V(20/res-n+1);
end

CA = [CAmirror CA]; %mirror V and CA for bTDC
V = [Vmirror V];

%% Burn Curve
mfb =(1-exp(-5*((CA-CA0)/burndur).^3)); %Wiebe curve

i=1; %loop counter
while mfb(i) < 0 %overwrite negative Wiebe points
    mfb(i) = 0;
    i = i + 1;
end

Qfuel = lhv * diff([mfb 1]); %heat addition from fuel is the differential of mfb

%% Combustion Model Initial Conditions

P = zeros(1,length(CA));
P(1) = Pint * (VBDC / V(1))^gamma;%assuming isentropic compression with constant specific heat from Pint

T = zeros(1,length(CA));
T(1) = P(1) * V(1) / mass / R; %ideal gas law

W = zeros(1,length(CA)); %indicated work [J]
Wcomp = (P(1) * V(1)-Pint * VBDC)/(1-gamma); %isentropic, constant gamma compression work

%% Heat Transfer Model Initial Conditions
Pmotor = Pint * (VBDC ./V).^gamma; %constant specific heat isentropic motoring curve
gasV = zeros(1,length(CA));
h = zeros(1,length(CA)); %Convection coefficients [W/m^2/K]
Qloss = zeros(1,length(CA)); %heat loss [J]

%% Time Step Calculations
cp = zeros(1,length(CA)); %constant pressure specific heat [J/kg/C]

for a = 2:length(CA)

    if T(a-1) < 2500
        gammahold = (-9.967e-12)*(T(a-1)+273.15)^3+(6.207e-8)*(T(a-1)+273.15)^2-(1.436e-4)*(T(a-1)+273.15)+1.396; %from SAE 2004-01-2996
    else
        gammahold = 1.2692; %beyond range of correleation assumed constant
    end
    
    cp(a-1) = gammahold * R / (gammahold - 1);
    
    gasV(a) = 2.28*S+0.00324/6*displacement*T(1)/Pmotor(1)/V(1)*(P(A)-Pmotor(a)); %gas velocity for reduced woshni from SAE 2004-01-2996
    h(a) = hmod*(3.4*3.26*min((4*V(a)/(pi*bore^2)),bore/2)^-0.2*P(a-1)^0.8*T(a-1)^-0.73*gasV(a)^0.8); %reduced woshni from SAE 2004-01-2996
    Qloss(a) = Area(a) * h(a)* (T(a-1)-Twall) * dt; %[J]
    
    T(a) = (Qfuel(a) + Qloss(a) - P(a-1) * (V(a) - V(a-1))) / ((cp(a-1) - R) * mass) + T(a-1); %1st law at constant pressure
    P(a) = mass * R * T(a) / V(a); %ideal gas law
    
    W(a) =  (P(a) + P(a-1)) * (V(a) - V(a-1)) / 2; % PdV work based on average pressure. this is inconsistant with 1st law analysis used to find T
end

%% Total Calculations
Wtot = sum(W) + Wcomp;
Qlosstot = sum(Qloss);
Pmax = max(P);

end
%%EOF