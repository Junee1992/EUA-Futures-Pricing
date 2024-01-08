%% Parameter Estimation
function [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett, vtt, fitresult_linear, fitresult_season, att_1, Ptt, Ptt_1] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names, LT, correlation, par_init)

[nobsn, ncontracts] = size(y);
if sum(contains(par_names, "phi")) > 0
    n_lag = sum(contains(par_names, "phi"))/ncontracts;
    serial = "yes";
else
    n_lag = 0; serial = "no";
end

% Parameter constraints
[A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag);
for i = find(contains(par_names, "phi"))'
    if par_init(i) == 0
        lb(i) = 0; ub(i) = 0;
    end
end

if detrend_price == "yes"
    % Detrend the price
    [y_detrend, trend, fitresult_linear] = linear_detrend(y);
    seasonality = "sinusoid_b";
    for i = 1:10
        [season, p, fitresult_season, gof] = deseasonalise(y_detrend, i, seasonality);
        rmse_seasonality(i) = gof.rmse;
    end
    no_season = find(rmse_seasonality == min(rmse_seasonality));
    [season, p, fitresult_season, gof] = deseasonalise(y_detrend, no_season, seasonality);
    y_deseason = y_detrend - season;
else
    trend = []; season = [];
    y_deseason = y;
    fitresult_linear = []; fitresult_season = [];
end

global save_ett save_ytt save_att save_vtt save_att_1 save_Ptt_1 save_Ptt save_vtt_star

% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 300, 'MaxFunEvals', 50000, 'algorithm', 'sqp');
if LT == "GBM"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v3_arp, par_init, [], [], [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
elseif LT == "OU"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v3_arp, par_init, A, b, [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
end

par_optim = par_optimised;
log_L_optim = log_L;
att = save_att;
ytt = save_ytt;
ett = save_ett;
vtt = save_vtt;
att_1 = save_att_1;
Ptt = save_Ptt;
Ptt_1 = save_Ptt_1;

end