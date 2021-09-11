%combustion optimization script,
%runs DIRECT and fmincon on a latin hypercube of starting points
%x = [mass, splinefit, CA0, burndur, phi]

%% DIRECT
vlb = [0.0007, 3.6, -20, 000, 0.175]; % lower bounds
vub = [0.0015, 9.7, 180, 200, 1.000]; % upper boudns

nc=3;% number of non-linear constraints

I=[ ];% set discrete variable

c_L=-inf*ones(nc,1);% non-linear constraint lower bounds (-infinity for negative null form)
c_U=zeros(nc,1);% non-linear constraint upper bounds (zero for negative null form)

GLOBAL.MaxIter=0;% set max. iterations to zero
GLOBAL.MaxEval=8000;% DIRECT will terminate after 1000 function evals

A=[ ]; b_L=[ ]; b_U=[ ];% place linear constraints, if any--otherwise, these remain empty
xopt = gclsolve('engine', 'nonlcon', vlb, vub, A, b_L, b_U, c_L, c_U, I, GLOBAL)

%combustion optimization script
%x = [mass, splinefit, CA0, burndur, phi]

%% fmincon hyper cube optimization

[mass_0 splinefit_0 CA0_0 burndur_0 Pmax_0 RI_0 flameT_0 TBDc_0...
    Qpistontot_0 amax_0 have_0 sdist_0] = get_scaling(  );

lb = [.0007/mass_0, 3.6/splinefit_0, -20/CA0_0, 0, 0.175]; %scaled lower bounds
ub = [.0015/mass_0, 9.7/splinefit_0, 180/CA0_0, 200/burndur_0, 1]; %scaled upper boudns
options = optimset('Algorithm','sqp','MaxFunEvals',2000);

%x0 = [.0003/mass_0,3.6/splinefit_0,-15/CA0_0,20/burndur_0,1;...
%      .0003/mass_0,3.6/splinefit_0,36/CA0_0,20/burndur_0,1];

points = 128;
X=lhsdesign(points,5,'criterion','maximin'); %determine test points from a hypercube

x0 = zeros(points,5);
xstar = zeros(points,5);
fstar = zeros(1,points);
flag = zeros(1,points);

for i = 1:points

    x0(i,:) = (ub - lb) .* X(i,:) + lb; %transform hypercube variables to scaled problem variables

    [xstar(i,:) fstar(i), flag(i), output(i)] = fmincon(@s_engine, x0(i,:), [], [], [], [], lb, ub, @s_nonlcon, options);
    
    if flag(i) < 0 %if optimizer failed
        fstar(i) = 0;
    end
    i
end

