%% Parameter Estimation
function [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett] = param_estim_arp(y, ttm, deltat, detrend_price, n_par, par_names, LT, correlation, par_init)

[nobsn, ncontracts] = size(y);
if sum(contains(par_names, "phi")) > 0
    n_lag = sum(contains(par_names, "phi"))/ncontracts;
    serial = "yes";
else
    n_lag = 0; serial = "no";
end

% Parameter constraints
[A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag);

if detrend_price == "yes"
    % Detrend the price
    [y_detrend, trend] = linear_detrend(y);
    seasonality = "sinusoid_b"
    for i = 1:10
        [season, p, fitresult, gof] = deseasonalise(y_detrend, i, seasonality);
        rmse_seasonality(i) = gof.rmse;
    end
    no_season = find(rmse_seasonality == min(rmse_seasonality));
    [season, p, fitresult, gof] = deseasonalise(y_detrend, no_season, seasonality);
    y_deseason = y_detrend - season;
else
    trend = []; season = [];
    y_deseason = y;
end

% Grid search
% init = grid_search(lb, ub, 3, 6, n_par, n_lag, ncontracts, LT, correlation);
% log_L = kf_v2_arp(init(1,:), par_names, y_deseason, deltat, ttm, LT, correlation, serial);
% % parpool(8);
% parfor i = 1:size(init,1)
%     log_L = kf_v2_arp(init(i,:), par_names, y_deseason, deltat, ttm, LT, correlation, serial);
%     logL(i) = log_L;
% end
% par_init = init(find(logL == max(logL)),:);

global save_ett save_ytt save_att save_vt

% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 300, 'MaxFunEvals', 50000);
% kf_v1_arp(par_init, par_names, y, deltat, ttm, LT, correlation, serial)
if LT == "GBM"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v1_arp, par_init, [], [], [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
elseif LT == "OU"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v2_arp, par_init, A, b, [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
end

par_optim = par_optimised;
log_L_optim = log_L;
att = save_att;
ytt = save_ytt;
ett = save_ett;
vtt = save_vt;


end