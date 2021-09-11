function out=idealHR_v1(step, Qin, bore, stroke, r, Apiston, Ahead, RPM)

%Find the ideal HR curve to maximize the balance between early
%combustion/heat loss and late combustion/low expansion work

%step = step size [deg CA]
%Qin = energy addition from fuel [J]
%bore and stroke [m]
%Apiston & Ahead area of piston and head [m^2]
%r = compression ratio

%% Prelim calcs
CA = linspace(0,180,180/step+1); %only expansion stroke
Vswept = pi*bore^2/4*stroke; %[m^3]
Vclear = Vswept/(r-1); %clearance vol [m^3]
mass = Vswept*1.2; %trapped mass [kg]
S = stroke * 2 * RPM / 60; %mean piston velocity [m/s]
dt = step / RPM / 30;


%% Determine Motoring Trace
V = Vclear + (Vswept - Vswept*cos(CA*pi/180))/2; %m^3
Pmotor(1) = 101.3 * (r)^1.4;
Pmotor = Pmotor(1)*(V(1)./V).^1.4; %Motor trace [kPa]
Tr = Pmotor(1)*V(1)/mass/287; %reference temp in K

%% Brute force ideal P

Twall = 410; %[K]
idealP = zeros(1,length(Pmotor));

for n = 1:length(CA)-1
    
Area = Apiston + Ahead + V(n)*4/bore;

guessP = linspace(Pmotor(n), 10*Pmotor(n),length(Pmotor));
guessT = guessP.*V(n)/mass/0.287; %temp in K
gasV = 2.28*S+0.00324/6*Vswept*Tr/Pmotor(1)/V(1)*(guessP-Pmotor(n));
h = 3.4*3.26*min((4*V/(pi*bore^2)),bore/2).^-0.2.*guessP.^0.8.*guessT.^-0.73.*gasV.^0.8;
work = guessP .* (V(n+1)-V(n));
Qloss = Area .* h .* (guessT-Twall) * dt;
net = work - Qloss;

[C,I]=max(net);
idealP(n) = guessP(I);

end

plot(CA,idealP,'r',CA,Pmotor)

out.CA=CA;
out.idealP=idealP;
out.Pmotor=Pmotor;
out.Qloss=Qloss;
out.work=work;
out.net=net;
out.guessP=guessP;

% subplot(3,3,1),plot(CA,guessP)
% subplot(3,3,2),plot(CA,guessT)
% subplot(3,3,3),plot(CA,gasV)
% subplot(3,3,4),plot(CA,h)
% subplot(3,3,5),plot(CA,work)
% subplot(3,3,6),plot(CA,Qloss)
% subplot(3,3,7),plot(CA,net)
% subplot(3,3,8),plot(CA,guessP)
% subplot(3,3,9),plot(CA,guessP)



%% Determine Ideal Variation from Motor Trace Pressure
%for n=1:length(CA)-1 % the -1 is to be able to calc the last dV
% n=1;
% Twall = 410; %[K]
% guessP = Pmotor(1);
% Tr = Pmotor(1)*V(1)/mass/287; %reference temp in K
% 
% prevnet = 0; %net work from previous iteration
% net = 1;

% %% Search for Max P
% %work up to max P ideal by steps
% 
% 
% while net > prevnet
%     
%     prevnet=net;
%     prevguess = guessP;
%     
%     guessP = prevguess + 2*Pmotor(n);
%     guessT = guessP*V(n)/mass/0.287; %temp in K
%     gasV = 2.28*S+0.00324*Vswept*Tr/Pmotor(1)/V(1)*(guessP-Pmotor(n));
%     h = 3.26*bore^-0.2*guessP^0.8*guessT^-0.55*gasV^0.8;
%     work = guessP * (V(n+1)-V(n));
%     Qloss = Area * h * (guessT-Twall) * dt;
%     net = work - Qloss;
%      
% end
% 
% %% Refine Max P
% %Search for max P ideal thats between steps used in "Search for Max P"
% %above
% 
% maxP = guessP
% minP = prevguess
% 
% while abs(1-maxP/minP) > 0.05
%     
%     prevnet=net;
%     
%     guessP = (maxP+minP)/2;
%     guessT = guessP*V(n)/mass/0.287; %temp in K
%     gasV = 2.28*S+0.00324*Vswept*Tr/Pmotor(1)/V(1)*(guessP-Pmotor(n));
%     h = 3.26*bore^-0.2*guessP^0.8*guessT^-0.55*gasV^0.8;
%     work = guessP * (V(n+1)-V(n));
%     Qloss = Area * h * (guessT-Twall) * dt;
%     net = work - Qloss;
%     
%     if net > prevent
%         
%     end
%     
% end
%     
% while abs(1-prevguess/guessP) < 0.05
%     guessP = (maxP + minP)/2;
%     guessT = guessP*V(n)/mass/287; %temp in K
%     gasV = 2.28*S+0.00324*Vswept*Tr/Pmotor(1)/V(1)*(Pguess-Pmotor(n));
%     h = 3.26*bore^-0.2*guessP^0.8*guessT^-0.55*gasV^0.8;
%     work = guessP * (V(n+1)-V(n));
%     Qloss = Area * h * (guessT-Twall) * dt;
%     net = work - Qloss;
%     if net > prevnet
%         prevnet = net;
% 
%     end
% end

    
    
%end

%% EOF