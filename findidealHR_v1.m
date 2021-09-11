function out=findidealHR(CAstep, Pstep, bore, stroke, r, Apiston, Ahead, RPM, hmod)

%Find the ideal HR curve to maximize the balance between early
%combustion/heat loss and late combustion/low expansion work

% version one finds the best us of a pressure change. However this is not
% necessarily the best use of a energy change because the pressure change
% can vary with a given energy input because of pressure, gamma and volume.

%%%%%%%%%%%%%%%%%%%%%%%
%ver 2: sweep instant energy release along CA
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%assumbtions:
%motor trace is isentropic
%semi constant gamma
%%%%%%%%%%%%%%%%%%%%%%%

%step = step size [deg CA]
%Pstep = Pressure change step [kpa]
%Qin = energy addition from fuel [J]
%bore and stroke [m]
%Apiston & Ahead area of piston and head [m^2]
%r = compression ratio
%hmod = h modifier to force heat transfer lower

%% Prelim calcs
CA = linspace(0,180,180/CAstep+1); %only expansion stroke
Vswept = pi*bore^2/4*stroke; %[m^3]
Vclear = Vswept/(r-1); %clearance vol [m^3]
mass = Vswept*1.2; %trapped mass [kg]
S = stroke * 2 * RPM / 60; %mean piston velocity [m/s]
dt = CAstep / RPM / 30;

%% Determine Motoring Trace
V = Vclear + (Vswept - Vswept*cos(CA*pi/180))/2; %m^3
Pmotor(1) = 101.3 * (r)^1.4; %[kpa]
Pmotor = Pmotor(1)*(V(1)./V).^1.4; %Motor trace [kPa]
Tr = Pmotor(1)*V(1)/mass/0.287; %reference temp in K

%% Brute force ideal P
% 1) try adding Pstep to each CA point
% 2) keep the point from above that produces the most work
% 3) repeat, keep occasional traces for different fueling rates

guessP = Pmotor; %stores pressure trace in loop %[kpa]
holdP = Pmotor; %stores best pressure trace between loops %[kpa]
bestP = Pmotor; %stores best pressure trace at the highest level, used to initialize next loop %[kpa]
out=Pmotor;
guessT = guessP.*V/mass/0.287; %initialize to motored [K]

Twall = 410; %[K]
Area = Apiston + Ahead + V.*4/bore;

gasV = zeros(1,length(Pmotor)); %initialize with zeros
h = zeros(1,length(Pmotor));
Qloss = zeros(1,length(Pmotor)); %causes no heat loss on first iteration, could move calcs to front of loop and adjust indicies accordingly


for a = 1:10 %arbitrary number of pressure search iterations
    
    bestwork = -inf; %reset work for new sweep
    
    for b = 1:length(CA)-1 %search along CA to find best place for pressure increase
                 
        guessP = bestP;
            
        for c =2:length(CA) %carry out calc of work for each pressure increase location
         
            if c > b %after pressure increase
                
                if guessT(c-1) < 2500
                    gamma = (-9.967e-12)*guessT(c-1)^3+(6.207e-8)*guessT(c-1)^2-(1.436e-4)*guessT(c-1)+1.396; %from SAE 2004-01-2996
                else
                    gamma = 1.2692; %beyond range of correleation assumed constant
                end
                
                guessP(c) = (-Qloss(c-1) * (gamma-1)/1000+V(c)*guessP(c-1))/(V(c)+gamma*(V(c)-V(c-1)));%[kpa]
                
            elseif c == b %at pressure increase
                guessP(c) = guessP(c) + Pstep; %[kpa]
            %else: it keeps the pre existing pressure trace from bestP    
            end 
               
            %these could be reduced to overwriting scalars to save memory
            guessT(c) = guessP(c)*V(c)/mass/0.287; %temp in K
            gasV(c) = 2.28*S+0.00324/6*Vswept*Tr/Pmotor(1)/V(1)*(guessP(c)-Pmotor(c)); %gas velocity for reduced woshni from SAE 2004-01-2996
            h(c) = hmod*3.4*3.26*min((4*V(c)/(pi*bore^2)),bore/2)^-0.2*guessP(c)^0.8*guessT(c)^-0.73*gasV(c)^0.8; %reduced woshni from SAE 2004-01-2996
            Qloss(c) = Area(c) * h(c)* (guessT(c)-Twall) * dt;
            
        end
        
        plot(CA,guessP)
        
        work = guessP .* diff([0,V]);
        
        if sum(work) > bestwork
            holdP = guessP; %keep pressure trace incase no other pressure increase locations in this loop outperform it %[kpa]
            bestwork = sum(work);
        end
        
    end
    
    bestP = holdP; %store best trace to start next iteration with %[kpa]
    
    out = [out;bestP]; %store all pressure traces. first line is motored
    
end
%% EOF
