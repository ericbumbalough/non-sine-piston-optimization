%combustion optimization script
%x = [mass, splinefit, CA0, burndur, phi]

clc

[mass_0 splinefit_0 CA0_0 burndur_0 Pmax_0 RI_0 flameT_0 TBDc_0...
    Qpistontot_0 amax_0 have_0 sdist_0] = get_scaling(  );

lb = [0/mass_0, 3.6/splinefit_0, -20/CA0_0, 0, 0.175]; %scaled lower bounds
ub = [.0015/mass_0, 9.7/splinefit_0, 180/CA0_0, 200/burndur_0, 1]; %scaled upper boudns
options = optimset('Display','iter','Algorithm','sqp');

x0 = xopt.x_k' ./ [mass_0 splinefit_0 CA0_0 burndur_0 1]; %use optimum values from DIRECT as start to fmincon and scale it

[xstar fstar] = fmincon(@s_engine, x0, [], [], [], [], lb, ub, @s_nonlcon, options)