function out=sweep(fuelstart,fuelstop)

%sweep fueling rates to determine best timing for inst combustion

Qin = linspace(fuelstart,fuelstop,20);
LocPeakP = zeros(1,length(Qin));
work = zeros(1,length(Qin));
loss = zeros(1,length(Qin));
T = zeros(181,length(Qin));

for n=1:length(Qin)
    
    store = idealHR(1,Qin(n),.086,.0946,12.5,.007,.01,2000);
    LocPeakP(n)=store.LocPeakP;
    work(n)=store.bestwork;
    loss(n)=store.bestloss;

    
end

out.loc=LocPeakP;
out.work=work;
out.Qin=Qin;
out.loss=loss;


