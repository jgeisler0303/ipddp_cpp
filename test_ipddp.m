
alg.tol=1e-7;
alg.maxiter=1000;
alg.mu=0; % 0 for automatic selection

alg.infeas=0;
%alg.infeas=1;

alg.show_progress= 1;

%% InvPend
[funcs, fp, bp] = dynamics_invpend(1);
[fp_opt, bp, trace, time] = ipddp(fp, bp, funcs, alg);
mex_ipddp_problem('InvPend', @dynamics_invpend)
[success, x, u]= ipddp_InvPend_mex(fp.x(:, 1), fp.u, alg);

tiledlayout(size(fp.x, 1), 1)
k= 1:fp.horizon+1;
for i= 1:size(fp.x, 1)
    nexttile
    plot(k, fp_opt.x(i, :), k, x(i, :))
end

%% Car
[funcs, fp, bp] = dynamics_car(1);
[fp_opt, bp, trace, time] = ipddp(fp, bp, funcs, alg);
mex_ipddp_problem('Car', @dynamics_car)
[success, x, u]= ipddp_Car_mex(fp.x(:, 1), fp.u, alg);

tiledlayout(size(fp.x, 1), 1)
k= 1:fp.horizon+1;
for i= 1:size(fp.x, 1)
    nexttile
    plot(k, fp_opt.x(i, :), k, x(i, :))
end

%% ConCar
[funcs, fp, bp] = dynamics_concar(1);
[fp_opt, bp, trace, time] = ipddp(fp, bp, funcs, alg);
mex_ipddp_problem('ConCar', @dynamics_concar)
[success, x, u]= ipddp_ConCar_mex(fp.x(:, 1), fp.u, alg);

tiledlayout(size(fp.x, 1), 1)
k= 1:fp.horizon+1;
for i= 1:size(fp.x, 1)
    nexttile
    plot(k, fp_opt.x(i, :), k, x(i, :))
end

%% Arm
[funcs, fp, bp] = dynamics_arm(1);
[fp_opt, bp, trace, time] = ipddp(fp, bp, funcs, alg);
mex_ipddp_problem('Arm', @dynamics_arm)
[success, x, u]= ipddp_Arm_mex(fp.x(:, 1), fp.u, alg);

tiledlayout(size(fp.x, 1), 1)
k= 1:fp.horizon+1;
for i= 1:size(fp.x, 1)
    nexttile
    plot(k, fp_opt.x(i, :), k, x(i, :))
end
