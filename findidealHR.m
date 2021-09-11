function out=findidealHR(CAstep, Qstep, bore, stroke, r, Apiston, Ahead, RPM, hmod)

%Find the ideal HR curve to maximize the balance between early
%combustion/heat loss and late combustion/low expansion work

%%%%%%%%%%%%%%%%%%%%%%%
%ver 2: sweep instant energy release along CA
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%assumbtions:
%motor trace is isentropic
%semi constant gamma
%what does my selection of reference temperature and pressure do to the
%   woshni correlation
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
% 1) try adding Qstep to each CA point
% 2) keep the point from above that produces the most work
% 3) repeat, keep occasional traces for different fueling rates

guessP = Pmotor; 
guessQ = zeros(1,length(Pmotor));%stores HR trace in loop %[J]
holdQ = guessQ; %stores best HR trace between loops %[J]
bestQ = guessQ; %stores best HR trace at the highest level, used to initialize next loop %[J]
out=guessQ;


Twall = 410; %[K]
Area = Apiston + Ahead + V.*4/bore;




for a = 1:10 %arbitrary number of pressure search iterations
    
    bestwork = -inf; %reset work for new sweep
    
    for b = 1:length(CA)-1 %search along CA to find best place for pressure increase
                 
        guessQ = bestQ;
        gasV = zeros(1,length(Pmotor)); %initialize with zeros
        h = zeros(1,length(Pmotor));
        Qloss = zeros(1,length(Pmotor)); %causes no heat loss on first iteration, could move calcs to front of loop and adjust indicies accordingly
        work = zeros(1,length(Pmotor));
        guessT = guessP.*V/mass/0.287; %initialize to motored [K

        for c =2:length(CA) %carry out calc of work for each pressure increase location
         
            if c > b %after pressure increase
                
                if guessT(c-1) < 2500
                    gamma = (-9.967e-12)*guessT(c-1)^3+(6.207e-8)*guessT(c-1)^2-(1.436e-4)*guessT(c-1)+1.396; %from SAE 2004-01-2996
                else
                    gamma = 1.2692; %beyond range of correleation assumed constant
                end
                
                guessP(c) = (guessQ(c-1)-Qloss(c-1) * (gamma-1)/1000+V(c)*guessP(c-1))/(V(c)+gamma*(V(c)-V(c-1)));%[kpa]
                
            elseif c == b %at Energy increase
                guessQ(c) = guessQ(c) + Qstep; %[kpa]
            %else: it keeps the pre existing pressure trace from bestP    
            end 
               
            %these could be reduced to overwriting scalars to save memory
            guessT(c) = guessP(c)*V(c)/mass/0.287; %temp in K
            gasV(c) = 2.28*S+0.00324/6*Vswept*Tr/Pmotor(1)/V(1)*(guessP(c)-Pmotor(c)); %gas velocity for reduced woshni from SAE 2004-01-2996
            h(c) = hmod*3.4*3.26*min((4*V(c)/(pi*bore^2)),bore/2)^-0.2*guessP(c)^0.8*guessT(c)^-0.73*gasV(c)^0.8; %reduced woshni from SAE 2004-01-2996
            Qloss(c) = Area(c) * h(c)* (guessT(c)-Twall) * dt; %[J]
            work(c) = guessP(c)*(V(c)-V(c-1));
            
        end
        
        %% fails to forget previous cycles. when walking along HR placement
        %% it keeps the previous data for everything except heat release.
        %% also units must be wrong because inputting 5 causes giant
        %% pressure changes.
        
        
        

        subplot(3,3,1),plot(CA,guessP)
        subplot(3,3,2),plot(CA,guessT)
        subplot(3,3,3),plot(CA,gasV)
        subplot(3,3,4),plot(CA,h)
        subplot(3,3,5),plot(CA,Qloss)
        subplot(3,3,6),plot(CA,guessQ) 
        subplot(3,3,7),plot(CA,Pmotor,'r',CA,guessP,'b')
        subplot(3,3,8),plot(CA,work)
       % subplot(2,2,3),plot(CA,guessQ)%         
%         work = guessP .* diff([0,V]);
                
        if sum(work) > bestwork
            holdQ = guessQ; %keep pressure trace incase no other pressure increase locations in this loop outperform it %[kpa]
            bestwork = sum(work);
        end
        
    end
    
    bestQ = holdQ; %store best trace to start next iteration with %[kpa]
    
    out = [out;bestQ]; %store all pressure traces. first line is motored
    
end
%% EOF
