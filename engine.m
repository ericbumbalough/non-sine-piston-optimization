function [negative_eff Wtot Pmax RI minflameT TBDC Qpistontot amax have sdist] = engine(in)

mass = in(1);
splinefit = in(2);
CA0 = in(3);
burndur = in(4);
phi = in(5);

%% List of Inputs
%Pint: intake pressure [Pa]
%Tint: intake temperature [C]
%splinefit: instantanious compression ratios at CA in CAforce
%CA0: CA at start of burn [degree]
%burndur: burn length parameter of Weibe function
%phi: normalized air fuel ratio

%% List of Outputs
%Wtot: indicated work without pumping loss [W]
%Pmax: maximum cylinder pressure [Pa]
%RI: ringing intensity. the sound power in the engine cylinder per area. [W/m^2]
%minflameT: minimum temperature during flame propogation. used to verify combustion is feasible [C]
%TBDC: the temperature at bottom dead center. assumed to be the exhaust temperature [C]
%mass: trapped mass [kg]
%Qpiston: heat loss throught piston [J]
%amax: maximum acceleration [m/s^2]
%have: CA averaged convection coefficient [W/m^2/K]
%snet: net stroke. the distance traveled by the piston [m]

%% Parameters
CR = 12.5; %compression ratio
displacement = 0.00055; %engine displacement [m^3]
CAforce = 90; %CA points that the spline fits at given vforce points
gamma = 1.4; %ratio of specific heats. Assumed constant for non firing conditions.
R = 286.9; %ideal gas constant for air. [J/kg/C]
bore = 0.085; %cylinder bore [m]
RPM = 2000;
Twall = 137; %this was from a Babajimapolis paper. find it. [C]
res = .1; %Crank angle resolution. 0.1 is >99% accurate compared to 0.00001 and runs in about .01 seconds
hmod = 4; %heat transfer modifier determined experimentally
beta = 2; %experimental scaling factor for RI. correlates maximum pressure rise rate to amplitude of pressure waves [s]

%% Geometry
VTDC = displacement /(CR - 1); % [m^3]
VBDC = VTDC * CR; % [m^3]
Vspline = splinefit * VTDC; %[m^3]

stroke = 4 * displacement / pi / bore^2; %[m]
S = stroke * 2 * RPM / 60; %mean piston velocity [m/s]
Apiston = pi * bore^2 / 4; %area of piston, assuming a flat piston [m^2]
Ahead = pi * bore^2 / 4; %area of head, assuming a flat head [m^2]
dt = res / RPM / 6; %time step [s]

CA = 0:res:180; %crank angle array %%%this will cause a mistake if res is not a factor of 180
V=spline([0 CAforce 180],[0 [VTDC Vspline VBDC] 0], CA); %fit volume as a spline with 0 slope endpoints
%perhaps i could constrain the spline to have only 1 inflection point?
%constrain with a maximum slope
%constrain piston maximum acceleration
%constrain piston maximum velocity

roundres = ceil(20/res)-1; %find a integer number for the length of the vectors 

CAmirror = zeros(1,roundres);
Vmirror = zeros(1,roundres);

for n=1:roundres %mirror V and CA for bTDC.
    CAmirror(n) = -CA(roundres-n+2);
    Vmirror(n) = V(roundres-n+2);
end

CA = [CAmirror CA]; %mirror V and CA for bTDC
V = [Vmirror V];

Area = Apiston + Ahead + V .* 4 / bore; %area used for heat transfer calculations [m^2]

%% Combustion Model Initial Conditions

Tint = (mass * R / 100000 / VBDC)^(gamma-1) * 300^gamma - 273.15; %[C] intake temp is an isntropic ideal gas process from ambient state
Pint = mass * R * (Tint+273.15) / VBDC;

P = zeros(1,length(CA));
P(1) = Pint * (VBDC / V(1))^gamma; %assuming isentropic compression with constant specific heat from Pint

T = zeros(1,length(CA));
T(1) = P(1) * V(1) / mass / R - 273.15; %ideal gas law [C]

W = zeros(1,length(CA)); %indicated work [J]
Wcomp = (P(1) * V(1)-Pint * VBDC)/(1-gamma); %isentropic, constant gamma compression work

%% Burn Curve
mfb =(1-exp(-5*((CA-CA0)/burndur).^3)); %Wiebe curve

i=1; %loop counter
while (i <= length(mfb)) && (mfb(i) < 0) %overwrite negative Wiebe points
    mfb(i) = 0;
    i = i + 1;
end

%%% 43.6 MJ/kg and 14.7 as stoiciometric AF should be made into parameters
lhv = phi * mass * 43e6 / 14.7;

Qfuel = lhv * diff([mfb 1]); %heat addition from fuel is the differential of mfb

%% Heat Transfer Model Initial Conditions
Pmotor = Pint * (VBDC ./V).^gamma; %constant specific heat isentropic motoring curve
gasV = zeros(1,length(CA)); %approximate gas velocity
h = zeros(1,length(CA)); %Convection coefficients [W/m^2/K]
Qloss = zeros(1,length(CA)); %heat loss [J]
Qpiston = zeros(1,length(CA)); %heat loss through piston [J]

%% Time Step Calculations
cp = zeros(1,length(CA)); %constant pressure specific heat [J/kg/C]
gammamin = 1.4; %holds minimum value of gamma for use in RI calculation

for a = 2:length(CA)

    if T(a-1) < 2500
        gammahold = (-9.967e-12)*(T(a-1)+273.15)^3+(6.207e-8)*(T(a-1)+273.15)^2-(1.436e-4)*(T(a-1)+273.15)+1.396; %from SAE 2004-01-2996
    else
        gammahold = 1.2626; %beyond range of correleation assumed constant
    end
    
    if gammahold < gammamin %hold lowest value of gamma for RI calculation
        gammamin = gammahold;
    end
    
    cp(a-1) = gammahold * R / (gammahold - 1);
    
    if P(a-1) > Pmotor(a) %do not account for additional heat transfer caused by burning prior to ignition
        gasV(a) = 2.28*S+0.00324/6*displacement*(Tint+273.15)/Pint/VBDC*(P(a-1)-Pmotor(a)); %gas velocity for reduced woshni from SAE 2004-01-2996
    else
        gasV(a) = 2.28*S;
    end
    
    h(a) = hmod*(3.4*3.26*min((4*V(a)/(pi*bore^2)),bore/2)^-0.2*(P(a-1)/10^3)^0.8*(T(a-1)+273.15)^-0.73*gasV(a)^0.8); %reduced woshni from SAE 2004-01-2996
    Qloss(a) = -Area(a) * h(a)* (T(a-1)-Twall) * dt; %[J]
    Qpiston(a) = -Apiston * h(a)* (T(a-1)-Twall) * dt; %[J]
    
    T(a) = (Qfuel(a) + Qloss(a) - P(a-1) * (V(a) - V(a-1))) / ((cp(a-1)) * mass) + T(a-1); %1st law at constant pressure
    P(a) = mass * R * (T(a) + 273.15) / V(a); %ideal gas law
    
    W(a) =  (P(a) + P(a-1)) * (V(a) - V(a-1)) / 2; % PdV work based on average pressure. this is inconsistant with 1st law analysis used to find T
end

%% Total Calculations
Wtot = sum(W) + Wcomp;
negative_eff = -Wtot / lhv;
Qpistontot = sum(Qpiston);
Pmax = max(P);

dPdtmax = max(diff(P))/dt/1e6; %kpa/msec
RI = beta^2 * dPdtmax^2 / 2 / gammamin / (max(P)/1000) * sqrt(gammamin * R * (max(T)+273.15)); %ringing intensity [W/m^2]
TBDC = T(end); %[C] CELCIUS
x = V ./ (bore^2 * pi / 4); %piston position [m]
amax = max(diff(x,2)) / dt^2; %maximum piston acceleration [m/s^2]
have = mean(h); %crank angle averaged convection coefficient [W/m^2/C]
sdist = stroke + x(1) - VTDC / (bore^2 * pi / 4); %total distance traveled by piston

%% Flame Temperature Calculations
% this section is to avoid finding an empty matix

flamestart = find(CA>=CA0,1);
flameend = find(CA>=CA0+burndur,1);

if max(CA) < CA0 + burndur
    flameend = length(T); %if burn is still happening at BDC search till the end of T
end

minflameT = min(T(flamestart:flameend)); %find the minimum temp in the burning range [C]

if max(CA) < CA0
    minflameT = 0; %set burn to fail if CA0 is outside of model range
end

end
%%EOF