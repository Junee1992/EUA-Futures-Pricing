%% Parameter Estimation
function [par_optim, log_L_optim, par_init, trend, season, att, ytt, ett] = param_estim(y, ttm, deltat, detrend_price, n_par, ncontracts, par_names, LT, correlation, max_lags)

% Initial estimation - assume no AR process
n_lag = 0;
serial = "no";

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
init = grid_search(lb, ub, 3, 6, n_par, n_lag, ncontracts, LT, correlation);
log_L = kf_v1(init(1,:), par_names, y_deseason, deltat, ttm, LT, correlation, serial);
% parpool(8);
parfor i = 1:size(init,1)
    log_L = kf_v1(init(i,:), par_names, y_deseason, deltat, ttm, LT, correlation, serial);
    logL(i) = log_L;
end
par_init = init(find(logL == max(logL)),:);

global save_ett save_ytt save_att save_vt

% Setting options for optimisation
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 300, 'MaxFunEvals', 50000);

if LT == "GBM"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v1, par_init, [], [], [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
elseif LT == "OU"
    [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v1, par_init, A, b, [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
end

% Checking AR
if max_lags > 0
    for j = 1:ncontracts
        ARmodel = estimate(arima('Constant', 0, 'ARLags', [], 'Distribution', 'Gaussian'), save_ett(:,j), 'display', 'off');
        ARmodelAIC(1, j) = ARmodel.summarize.AIC;
        for l = 1:max_lags
            ARmodel = estimate(arima('Constant', 0, 'ARLags', 1:l, 'Distribution', 'Gaussian'), save_ett(:,j), 'display', 'off');
            AR_coef{l+1,j} = ARmodel.AR;
            ARmodelAIC(l+1,j) = ARmodel.summarize.AIC;
        end
    end
    contractNames = arrayfun(@(x) sprintf('C%d', x), 1:ncontracts, 'UniformOutput', false);
    AICsummary = array2table([[0:max_lags]', ARmodelAIC], 'VariableNames', ['Lags' contractNames]);
    for j = 1:ncontracts
        lags_temp(j) = find(ARmodelAIC(:,j) == min(ARmodelAIC(:,j)))-1;
    end

    n_lag = max(lags_temp); % number of lags for serial correlation of measurement errors
    m_init = [];
    for j = 1:ncontracts
        m_init(j,:) = [cell2mat(AR_coef{n_lag+1,j})];
    end
    par_optimised = [par_optimised reshape(m_init', 1, size(m_init,1)*size(m_init,2))];
else
    par_optimised = par_optimised;
end

par_optim = par_optimised;
log_L_optim = log_L;
att = save_att;
ytt = save_ytt;
ett = save_ett;
vtt = save_vt;

% if n_lag > 0
%     serial = "yes";
%     par_names = define_parameters(LT, ncontracts, correlation, n_lag)';
%     n_par = size(par_names, 1);
%     [A, b, lb, ub] = param_constr(n_par, ncontracts, LT, correlation, n_lag);
%     m_init = [];
%     for j = 1:ncontracts
%         m_init = [m_init cell2mat(AR_coef{n_lag+1,j})];
%     end
%     par_init_ar = [par_optimised m_init];
%     global save_vt
%     [par_optimised, log_L, exitflag, output, lambda, grad, hessian] = fmincon(@kf_v1, par_init_ar, A, b, [], [], lb, ub, [], options, par_names, y_deseason, deltat, ttm, LT, correlation, serial);
% 
% end

